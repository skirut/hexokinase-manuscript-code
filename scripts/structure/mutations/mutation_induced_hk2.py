"""MODELLER mutation modeling for HK2 variants.

Purpose
-------
Applies the HK2 mutation set encoded in the script, handles the residue
numbering offset between the model and HK2 numbering, and writes
corrected mutant structures.

Inputs
------
- ``2nzt.pdb``
- Chain A and residue numbering defined in the ``mutations`` list.

Mutations encoded in this script
--------------------------------
T680M, D861Y, T863A, and D84Y.

Outputs
-------
- ``hk2_T680M.pdb``
- ``hk2_D861Y.pdb``
- ``hk2_T863A.pdb``
- ``hk2_D84Y.pdb``

Run from
--------
Run from the directory containing ``2nzt.pdb``.

Usage
-----
``python mutation_induced_hk2.py``

Dependencies
------------
MODELLER and Python.
"""

from modeller import Environ, Selection
from modeller.scripts import complete_pdb
import os

env = Environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

input_pdb = '2nzt.pdb'

# Mutations: (original_residue_number_in_HK2, chain_ID, target_residue_name, label_suffix)
mutations = [
    (680, 'A', 'MET', 'T680M'),
    (861, 'A', 'TYR', 'D861Y'),
    (863, 'A', 'ALA', 'T863A'),
    (84,  'A', 'TYR', 'D84Y'),
]

start_residue_number = 15  # First residue in HK2 (PDB) numbering
offset = start_residue_number - 1

def renumber_pdb_residues(pdb_file, offset, output_file):
    with open(pdb_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith(('ATOM', 'HETATM')):
                res_num = int(line[22:26])
                new_res_num = res_num + offset
                line = f"{line[:22]}{new_res_num:4d}{line[26:]}"
            fout.write(line)

for hk2_resi, chain, new_res, suffix in mutations:
    modeller_resi = hk2_resi - offset
    print(f"\nMutating HK2 residue {hk2_resi} (model residue {modeller_resi}) in chain {chain} to {new_res}...")

    mdl = complete_pdb(env, input_pdb)

    try:
        res = mdl.chains[chain].residues[str(modeller_resi)]
    except KeyError:
        print(f"Residue {modeller_resi} in chain {chain} not found in structure.")
        continue

    sel = Selection(res)
    sel.mutate(residue_type=new_res)

    # Verification
    mutated_res = mdl.chains[chain].residues[str(modeller_resi)]
    actual_res = mutated_res.name

    if actual_res.upper() == new_res.upper():
        print(f"Mutation successful: Residue {hk2_resi} (model resi {modeller_resi}) is now {actual_res}")
    else:
        print(f"Mutation failed: Expected {new_res}, but found {actual_res}")

    # Save model
    temp_pdb = f"temp_{suffix}.pdb"
    final_pdb = f"hk2_{suffix}.pdb"
    mdl.write(file=temp_pdb)

    # Fix residue numbering back to HK2
    renumber_pdb_residues(temp_pdb, offset, final_pdb)
    os.remove(temp_pdb)

    print(f"Saved corrected structure with HK2 numbering: {final_pdb}")
