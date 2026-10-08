# External Software

The Hexokinase study combines Python-based trajectory analysis with
specialized molecular modeling, molecular dynamics, and structural
analysis software.

These tools are not bundled with this repository.

## Molecular dynamics and system preparation

| Software | Version used in the working environment | Role |
|---|---:|---|
| HTMD | 2.0.3 | Molecular dynamics system preparation and workflow |
| HTMD-deps | 2.3.12 | HTMD dependency package |
| ACEMD | 3.7.6 | Molecular dynamics simulation engine |
| AmberTools | 22.0 | Molecular mechanics parameterization and system preparation |
| MoleculeKit | 1.8.36 | Molecular structure handling |
| Packmol | 20.010 | Molecular system construction |

## Structural analysis

| Software | Version / status | Role |
|---|---:|---|
| NACCESS | 2.1.1 | Solvent-accessible surface analysis |
| MolLayers | Local source package | Molecular surface/layer analysis |

The MolLayers implementation used in the study was maintained as a
local Python source package and was not installed as a standard
PyPI/Conda package.

## Homology modelling

| Software | Version / status | Role |
|---|---:|---|
| MODELLER | 10.6 available in local installation | Homology modelling and mutation modelling |

MODELLER is a separately distributed software package and is not
included in this repository.

## Reproducibility note

The exact installation procedure and licensing requirements of
specialized third-party software may differ from those of the Python
packages listed in `requirements.txt` and `environment.yml`.

The repository therefore separates the portable Python analysis
environment from specialized external software.
