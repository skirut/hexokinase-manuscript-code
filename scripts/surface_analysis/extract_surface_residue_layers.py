"""Surface residue composition extraction from MolLayers output.

Purpose
-------
Matches MolLayers ``.RTP``/surface PDB outputs and counts protein and
surface residues by physicochemical class. The results are combined
into a single Excel workbook.

Inputs
------
- MolLayers output under the configured
  ``results/surface_analysis/mol_layers`` directory.
- Corresponding ``_surf1.pdb`` surface files.

Outputs
-------
- ``results/tables/combined_surface_analysis.xlsx``

Run from
--------
Run from the repository root with ``HEXOKINASE_DATA_ROOT`` configured.

Usage
-----
``python scripts/surface_analysis/extract_surface_residue_layers.py``

Dependencies
------------
Python, pandas, and ``scripts.utils.paths``.
"""

import os
import glob
import pandas as pd

from scripts.utils.paths import MOL_LAYERS_DIR, TABLES_DIR

# Residue categories
NONPOLAR = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY"}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
POSITIVE = {"LYS", "ARG", "HIS"}
NEGATIVE = {"ASP", "GLU"}

# Paths
rtp_dir = MOL_LAYERS_DIR
pdb_dir = MOL_LAYERS_DIR
output_file = TABLES_DIR / "combined_surface_analysis.xlsx"

# Store results
all_data = []

# --- Helper Function ---
def count_atoms_on_surface(pdb_file, surface_res_ids):
    count = 0
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                res_name = line[17:20].strip()
                res_no = line[22:26].strip()
                chain = line[21].strip()
                res_id = f"{res_name}_{res_no}_{chain}"
                if res_id in surface_res_ids:
                    count += 1
    return count

# Match each .RTP with its corresponding surf1.pdb file
for rtp_file in glob.glob(os.path.join(rtp_dir, "*.RTP")):
    pdb_id = os.path.basename(rtp_file).split(".")[0]
    pdb_file = os.path.join(pdb_dir, f"{pdb_id}_surf1.pdb")
    if not os.path.isfile(pdb_file):
        print(f"Skipping {pdb_id}: PDB file not found.")
        continue

    # --- Parse RTP File ---
    total_residues = 0
    surface_residues = set()
    polar_surf = set()
    nonpolar_surf = set()
    pos_surf = set()
    neg_surf = set()

    with open(rtp_file) as f:
        for line in f:
            if line.startswith("#") or len(line.strip().split()) < 4:
                continue
            res_name, res_no, chain, is_surface = line.strip().split()[:4]
            res_id = f"{res_name}_{res_no}_{chain}"
            total_residues += 1

            if is_surface == "1":
                surface_residues.add(res_id)
                if res_name in NONPOLAR:
                    nonpolar_surf.add(res_id)
                elif res_name in POLAR:
                    polar_surf.add(res_id)
                elif res_name in POSITIVE:
                    pos_surf.add(res_id)
                elif res_name in NEGATIVE:
                    neg_surf.add(res_id)

    # --- Parse PDB File ---
    atom_count = 0
    residue_ids = set()
    nonpolar_all = set()
    polar_all = set()
    pos_all = set()
    neg_all = set()

    with open(pdb_file) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_count += 1
            res_name = line[17:20].strip()
            res_no = line[22:26].strip()
            chain = line[21].strip()
            res_id = f"{res_name}_{res_no}_{chain}"
            residue_ids.add(res_id)

            if res_name in NONPOLAR:
                nonpolar_all.add(res_id)
            elif res_name in POLAR:
                polar_all.add(res_id)
            elif res_name in POSITIVE:
                pos_all.add(res_id)
            elif res_name in NEGATIVE:
                neg_all.add(res_id)

    # --- Count Surface Atoms ---
    surface_atom_count = count_atoms_on_surface(pdb_file, surface_residues)

    # --- Combine Data ---
    row = [
        pdb_id,
        atom_count,
        surface_atom_count,
        len(residue_ids),
        total_residues,
        len(surface_residues),
        len(nonpolar_all), len(nonpolar_surf),
        len(polar_all), len(polar_surf),
        len(pos_all), len(pos_surf),
        len(neg_all), len(neg_surf)
    ]
    all_data.append(row)

# --- Save to Excel ---
columns = [
    "Isoform (PDB ID)",
    "Total No. of Atoms in Protein",
    "Total No. of Atoms on Surface",
    "Total No. of Residues in Protein (PDB)",
    "Total No. of Residues in Protein (RTP)",
    "Total No. of Residues on Surface (RTP)",
    "No. of Nonpolar Residues in Protein",
    "No. of Nonpolar Residues on Surface",
    "No. of Polar Residues in Protein",
    "No. of Polar Residues on Surface",
    "No. of Positively Charged Residues in Protein",
    "No. of Positively Charged Residues on Surface",
    "No. of Negatively Charged Residues in Protein",
    "No. of Negatively Charged Residues on Surface"
]

df = pd.DataFrame(all_data, columns=columns)
df.to_excel(output_file, index=False)

print(f"\nAnalysis complete. Excel saved to: {output_file}")


