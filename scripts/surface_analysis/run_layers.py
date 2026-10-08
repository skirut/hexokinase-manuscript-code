#!/usr/bin/env python3
"""MolLayers surface-layer calculation for modeled Hexokinase structures.

Purpose
-------
Runs the local MolLayers engine on every PDB file in the configured
modeled-structure directory with all-layer peeling enabled.

Inputs
------
- PDB files under ``$HEXOKINASE_DATA_ROOT/modelled_structures``
- The local MolLayers Python package/source used by the study.

Outputs
-------
- One ``*_layers.txt`` file per input PDB under
  ``results/surface_analysis/mol_layers/``.

Run from
--------
Run from the repository root with the local MolLayers package available
to Python and ``HEXOKINASE_DATA_ROOT`` configured.

Usage
-----
``python scripts/surface_analysis/run_layers.py``

Dependencies
------------
Python, the local MolLayers package, and ``scripts.utils.paths``.
"""


import os
from MolLayers import engine

from scripts.utils.paths import MODELLED_DIR, MOL_LAYERS_DIR

SOURCE_PATH = MODELLED_DIR
OUTPUT_PATH = MOL_LAYERS_DIR

OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

for pdb_file in os.listdir(SOURCE_PATH):
    if pdb_file.endswith(".pdb"):
        pdb_path = SOURCE_PATH / pdb_file

        data = engine.Layers(
            str(pdb_path),
            peel_all_layers=True
        )

        output_file = OUTPUT_PATH / pdb_file.replace(
            ".pdb",
            "_layers.txt"
        )

        with open(output_file, "w") as f:
            f.write(str(data))

print(
    "MolLayers analysis complete for all PDB files. "
    f"Results saved in: {OUTPUT_PATH}"
)

