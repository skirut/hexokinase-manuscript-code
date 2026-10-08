#!/usr/bin/env python3
"""Centralized data and results paths for the Hexokinase repository.

Purpose
-------
Loads ``config/systems.yaml`` and resolves the study data directories
from the ``HEXOKINASE_DATA_ROOT`` environment variable. It also defines
the repository's standardized results directories and creates them when
requested.

Inputs
------
- ``config/systems.yaml``
- Environment variable ``HEXOKINASE_DATA_ROOT``. The current module
  requires this variable to be set.

Outputs
-------
Provides reusable Path objects including ``APO_DIR``, ``BOUND_DIR``,
``MODELLED_DIR``, ``LIGAND_DIR``, and the analysis-specific results
directories.

Run from
--------
Imported by analysis scripts from the repository root or any context
where the repository package is importable.

Usage
-----
Set the data root before running analysis scripts, for example:
``export HEXOKINASE_DATA_ROOT=/path/to/hexokinase_data``

Dependencies
------------
Python standard library and PyYAML.
"""


from pathlib import Path
import os
import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = REPO_ROOT / "config" / "systems.yaml"


def load_config():
    with open(CONFIG_FILE, "r") as f:
        return yaml.safe_load(f)


CONFIG = load_config()

DATA_ROOT_ENV = os.environ.get("HEXOKINASE_DATA_ROOT")

if not DATA_ROOT_ENV:
    raise RuntimeError(
        "HEXOKINASE_DATA_ROOT is not set. "
        "Set it to the directory containing structures/, G6P_structure/, "
        "modelled_structures/, and ligand_regulatory_site/."
    )

DATA_ROOT = Path(DATA_ROOT_ENV).expanduser().resolve()
APO_DIR = DATA_ROOT / CONFIG["paths"]["apo"]
BOUND_DIR = DATA_ROOT / CONFIG["paths"]["bound"]
MODELLED_DIR = DATA_ROOT / CONFIG["paths"]["modelled_structures"]
LIGAND_DIR = DATA_ROOT / CONFIG["paths"]["ligand_regulatory_site"]

RESULTS_ROOT = REPO_ROOT / CONFIG["results_root"]

FIGURES_DIR = RESULTS_ROOT / "figures"
TABLES_DIR = RESULTS_ROOT / "tables"

LIGAND_INTERACTIONS_DIR = RESULTS_ROOT / "ligand_interactions"
RMSD_RMSF_RG_DIR = RESULTS_ROOT / "rmsd_rmsf_rg"
DSSP_DIR = RESULTS_ROOT / "dssp"
BINDING_ENERGY_DIR = RESULTS_ROOT / "binding_energy"
DIM_REDUCTION_DIR = RESULTS_ROOT / "dimensionality_reduction"
SURFACE_ANALYSIS_DIR = RESULTS_ROOT / "surface_analysis"
MOL_LAYERS_DIR = SURFACE_ANALYSIS_DIR / "mol_layers"


def ensure_results_dirs():
    for directory in [
        RESULTS_ROOT,
        FIGURES_DIR,
        TABLES_DIR,
        LIGAND_INTERACTIONS_DIR,
        RMSD_RMSF_RG_DIR,
        DSSP_DIR,
        BINDING_ENERGY_DIR,
        DIM_REDUCTION_DIR,
        SURFACE_ANALYSIS_DIR,
        MOL_LAYERS_DIR,
    ]:
        directory.mkdir(parents=True, exist_ok=True)

