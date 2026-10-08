#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Parameterize G6P at multiple pH values for Amber/GAFF2 and generate ACPYPE outputs
INPUT_SDF="G6P.sdf"

# pH values you want to model
PH_VALUES=(2.0 3.0 5.0 7.0)

for PH in "${PH_VALUES[@]}"; do
    # Integer pH for filenames: 2.0 -> 2, 7.0 -> 7
    PH_INT="${PH%.*}"
    BASE_NAME="g6p_pH${PH_INT}"

    echo "=== Working on pH ${PH} (PH_INT=${PH_INT}) -> base ${BASE_NAME} ==="

    echo "[1] Generate 3D structure and add hydrogens (OpenBabel)"
    # No -p; avoid buggy pH-dependent protonation for this polyphosphate
    obabel "$INPUT_SDF" -O "${BASE_NAME}_3D_h.mol2" --gen3D -h \
        || { echo "obabel 3D/H generation failed for pH ${PH}"; exit 1; }

    echo "[2] (Optional) check and handle EP (extra point) atoms"
    if grep -q " EP " "${BASE_NAME}_3D_h.mol2"; then
        echo ">>> Found 'EP' or extra-point entries in ${BASE_NAME}_3D_h.mol2. These must be removed or handled."
        echo ">>> Saving original and creating a cleaned version without EP."
        cp "${BASE_NAME}_3D_h.mol2" "${BASE_NAME}_3D_h_withEP.mol2"
        awk '!/ EP /{print}' "${BASE_NAME}_3D_h.mol2" > "${BASE_NAME}_noEP.mol2"
        ANTECH_INPUT="${BASE_NAME}_noEP.mol2"
    else
        ANTECH_INPUT="${BASE_NAME}_3D_h.mol2"
    fi

    echo "[3] Decide net charge (-nc) based on pH"
    # Approximate scheme for G6P:
    #   pH 2 ~ 0 (more protonated)
    #   pH 3 ~ 0 (still mostly protonated)
    #   pH 5 ~ -2 (deprotonated phosphate)
    #   pH 7 ~ -2
    NC=0
    case "${PH_INT}" in
        2) NC=0 ;;
        3) NC=0 ;;
        5) NC=-2 ;;
        7) NC=-2 ;;
        *) echo "No net charge rule defined for pH=${PH}"; exit 1 ;;
    esac
    echo "    -> Using net charge NC=${NC} for pH ${PH}"

    echo "[4] Run antechamber to assign GAFF2 atom types and AM1-BCC charges (G6P, pH=${PH})"
    # -nc ${NC}: enforce pH-dependent net charge
    # -rn G6P: residue name G6P (matches PDB ligand)
    antechamber -i "${ANTECH_INPUT}" -fi mol2 \
                -o "${BASE_NAME}_ante.mol2" -fo mol2 \
                -c bcc -s 2 -at gaff2 -nc "${NC}" -rn G6P \
                || { echo "antechamber failed for pH ${PH} — inspect output"; exit 1; }

    echo "[INFO] Check above antechamber output for 'net charge' confirmation."

    echo "[5] Check the antechamber mol2 for unknown atom types / EP"
    if grep -q " EP " "${BASE_NAME}_ante.mol2"; then
        echo "ERROR: antechamber output still contains EP. Stop and inspect the ligand file for pH ${PH}." >&2
        exit 2
    fi

    echo "[6] Generate frcmod for missing parameters (parmchk2)"
    parmchk2 -i "${BASE_NAME}_ante.mol2" -f mol2 -o "${BASE_NAME}.frcmod" -a Y \
        || { echo "parmchk2 failed for pH ${PH}"; exit 1; }

    echo "[7] (Optional) ligand geometry optimization with sqm"
    if command -v sqm >/dev/null 2>&1; then
        echo "NOTE: sqm usually uses an input (.in) file. Skipping automatic sqm optimization for now."
        # If you want proper sqm optimization later, you’d write a .in file and call:
        # sqm -i ${BASE_NAME}.in -o ${BASE_NAME}.out
    else
        echo "sqm not found; skipping ligand QM geometry optimization."
    fi

    echo "[8] Build AMBER prmtop/inpcrd using tleap (GAFF2) for pH ${PH}"
    tleap -f - <<EOF > "${BASE_NAME}_tleap.log" 2>&1
source leaprc.gaff2
loadamberparams ${BASE_NAME}.frcmod
lig = loadmol2 ${BASE_NAME}_ante.mol2
saveamberparm lig ${BASE_NAME}.prmtop ${BASE_NAME}.inpcrd
saveoff lig ${BASE_NAME}.lib
quit
EOF

    # Only stop on REAL Leap issues, not just the summary "Errors = 0; Warnings = 0"
    if grep -Ei "Error:|WARNING:|Unknown|missing parameter" "${BASE_NAME}_tleap.log" >/dev/null; then
        echo "tleap produced warnings/errors for pH ${PH}; inspect ${BASE_NAME}_tleap.log"
        grep -nEi "Error:|WARNING:|Unknown|missing parameter" "${BASE_NAME}_tleap.log" || true
        echo "Stopping due to tleap warnings. Fix ligand params/atom types first."
        exit 3
    fi

    echo "[9] Run ACPYPE on the ANTECHAMBER-produced mol2 (for GROMACS/CHARMM, etc.)"
    # This will create a directory like: g6p_pH2.acpype/, g6p_pH3.acpype/, ...
    acpype -i "${BASE_NAME}_ante.mol2" -b "${BASE_NAME}" -o all \
        || { echo "acpype failed for pH ${PH}"; exit 1; }

    echo "[DONE] Generated files for pH=${PH}:"
    ls -1 "${BASE_NAME}."* "${BASE_NAME}_ante.mol2" "${BASE_NAME}.frcmod" "${BASE_NAME}_tleap.log" || true
    echo "ACPYPE directory contents (${BASE_NAME}.acpype):"
    ls -1 "${BASE_NAME}.acpype" || true
    echo "-----------------------------------------------"
done

echo "All done for G6P at pH 2, 3, 5, 7."

