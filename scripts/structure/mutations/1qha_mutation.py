"""MODELLER mutation modeling for HK1 G87E.

Purpose
-------
Loads the HK1 starting structure, completes the PDB with MODELLER,
mutates residue 87 on chain A from its starting residue to GLU, performs
the local optimization steps in the script, and writes the mutated
structure.

Inputs
------
- ``1qha.pdb``
- Chain A residue 87 in the PDB numbering used by the script.

Outputs
-------
- ``1qha_G87E.pdb``

Run from
--------
Run from the directory containing ``1qha.pdb``.

Usage
-----
``python 1qha_mutation.py``

Dependencies
------------
MODELLER.
"""

# mutate_1qha_G87E.py
from modeller import Environ, Selection, log
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics

# --- settings you may edit ---
input_pdb  = '1qha.pdb'   # your starting structure
output_pdb = '1qha_G87E.pdb'
chain_id   = 'A'          # change if your residue is on a different chain
resnum     = '87'         # PDB numbering (string allows insertion codes like '87A')
target_aa  = 'GLU'        # target residue type (MODELLER 3-letter code)
# -----------------------------

log.none()  # or log.minimal() for more verbosity

env = Environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# Read and complete the PDB (adds missing atoms consistent with topology)
mdl = complete_pdb(env, input_pdb)

# Basic sanity checks
if chain_id not in [c.name for c in mdl.chains]:
    raise ValueError(f"Chain '{chain_id}' not found. Chains present: {[c.name for c in mdl.chains]}")

chain = mdl.chains[chain_id]
if resnum not in chain.residues:
    raise ValueError(f"Residue {resnum} not found on chain {chain_id}. "
                     f"Available example indices: {[r.index for r in chain.residues[:5]]} "
                     f"and names: {[r.name for r in chain.residues[:5]]}")

# Perform mutation
res = chain.residues[resnum]
sel = Selection(res)
sel.mutate(residue_type=target_aa)

# Rebuild restraints and do a light relaxation to resolve clashes
mdl.restraints.make(mdl, restraint_type='stereo', spline_on_site=False)

# Conjugate gradients minimization
cg = conjugate_gradients(output='REPORT')
cg.optimize(mdl, max_iterations=300, min_atom_shift=0.01)

# Short MD (optional but helpful)
md = molecular_dynamics(output='REPORT')
md.optimize(mdl, temperature=300, max_steps=500)

# Final CG to settle
cg.optimize(mdl, max_iterations=200, min_atom_shift=0.01)

# Save the mutated model
mdl.write(file=output_pdb)
print(f"Done. Wrote mutated structure to: {output_pdb}")


