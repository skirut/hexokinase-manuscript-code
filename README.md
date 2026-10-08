# Hexokinase Manuscript – Computational Analysis Code

## Associated Manuscript

**Understanding the structural and functional dynamics of Hexokinase Isoforms in Cancer Metabolism using Molecular Dynamics simulation studies**

**Authors:** Kiruthiga Shankar Kumar¹, Naga Bhushana Rao Karampudi¹  
¹ Department of Biological Sciences, School of Engineering and Sciences, SRM University-AP, Amaravati, Andhra Pradesh, India

Computational analysis and molecular dynamics workflows associated with the study:

**"Understanding the structural and functional dynamics of Hexokinase Isoforms in Cancer Metabolism using Molecular Dynamics simulation studies"**

## About This Repository

This repository contains the source code, configuration files, documentation, and computational workflow associated with the above manuscript. It is intended to provide a version-controlled record of the computational procedures used in the study and to support reproducibility of the reported analyses.

The repository does not contain large molecular-dynamics trajectories, topology files, coordinate files, restart files, or generated computational outputs. These files are intentionally excluded because of their size and because the repository is intended to preserve the analysis workflow and source code.

## Overview

This repository contains the computational scripts and workflow documentation used to investigate the structural, dynamic, evolutionary, and ligand-interaction properties of the four major human hexokinase isoforms:

- HK1
- HK2
- HK3
- HK4

The study combines structural modeling, mutation analysis, surface and residue-layer analysis, molecular dynamics simulation, trajectory analysis, ligand-interaction analysis, binding-energy analysis, dimensionality reduction, sequence/evolutionary analysis, and visualization.

The repository is intentionally limited to **source code, configuration files, documentation, and software-environment specifications**. Large molecular dynamics trajectories, topology files, coordinate files, simulation restarts, generated analysis outputs, and other computational data are not included.

---

## Structures Studied

| Isoform | PDB ID |
|---|---|
| HK1 | 1QHA |
| HK2 | 2NZT |
| HK3 | 3HM8 |
| HK4 | 1V4S |

The structures listed above were used as the starting structural references for the corresponding hexokinase isoforms.

---

## Mutation Analysis

The study includes cancer-associated and structurally investigated mutations in the hexokinase isoforms:

| Isoform | Mutations |
|---|---|
| HK1 | G87E |
| HK2 | D84Y, D861Y, T680M |
| HK3 | T546M, D867N, D867Y |
| HK4 | D205N |

Mutation-modeling scripts are provided under:

```text
scripts/structure/mutations/
```

---

## Computational Workflow

The repository follows the major computational stages used in the study:

1. Structure preparation and homology modeling
2. Mutation modeling and structural preparation
3. Structural validation
4. Sequence and evolutionary analysis
5. Surface and residue-layer analysis
6. G6P ligand preparation and parameterization
7. Molecular dynamics system preparation
8. Molecular dynamics simulation
9. RMSD, RMSF, and radius-of-gyration analysis
10. Secondary-structure analysis
11. Protein-G6P hydrogen-bond and contact analysis
12. Center-of-mass analysis
13. Binding-energy analysis
14. t-SNE and UMAP dimensionality-reduction analysis
15. Visualization and interpretation of computational results

The scripts are organized according to these workflow stages.

---

## Repository Structure

```text
hexokinase_simulation/
│
├── config/
│   └── systems.yaml
│
├── docs/
│   ├── script_guide.md
│   └── workflow.md
│
├── environment/
│   ├── environment.yml
│   ├── external_software.md
│   └── requirements.txt
│
├── LICENSE
├── README.md
│
└── scripts/
    │
    ├── __init__.py
    │
    ├── binding_energy/
    │   ├── binding_energy.py
    │   └── energy.py
    │
    ├── dimensionality_reduction/
    │   ├── global_tsne_umap_integrated_annotated_consistent.py
    │   └── plot.py
    │
    ├── dssp/
    │   └── dssp_hk.py
    │
    ├── ligand_interactions/
    │   ├── com.py
    │   └── hbond_com_contacts.py
    │
    ├── ligand_preparation/
    │   └── g6p_param.sh
    │
    ├── md_analysis/
    │   ├── rmsd_rmsf_rg.py
    │   └── vmd/
    │       ├── load_molecule.tcl
    │       ├── rmsd.tcl
    │       ├── rmsf.tcl
    │       └── rmsg.tcl
    │
    ├── sequence/
    │   └── visualize_tree.py
    │
    ├── simulation/
    │   ├── extended_simulation.py
    │   ├── hexokinase_g6p.py
    │   ├── hk_sim.py
    │   └── mutants_sim.py
    │
    ├── structure/
    │   ├── modeller/
    │   │   ├── 1qha_modeller.py
    │   │   ├── 1v4s_modeller.py
    │   │   ├── 2nzt_modeller.py
    │   │   └── 3hm8_modeller.py
    │   │
    │   └── mutations/
    │       ├── 1qha_mutation.py
    │       ├── mutation_g6p.py
    │       └── mutation_induced_hk2.py
    │
    ├── surface_analysis/
    │   ├── extract_surface_residue_layers.py
    │   ├── hk_layers_analysis.py
    │   └── run_layers.py
    │
    └── utils/
        ├── __init__.py
        └── paths.py
```

---

## Script Categories

### Structure Modeling

Located in:

```text
scripts/structure/modeller/
```

These scripts perform structural preparation and homology modeling for the hexokinase isoforms.

Included structures:

- HK1
- HK2
- HK3
- HK4

### Mutation Modeling

Located in:

```text
scripts/structure/mutations/
```

These scripts generate and prepare the investigated mutant structures and mutation-specific systems.

### Surface and Residue-Layer Analysis

Located in:

```text
scripts/surface_analysis/
```

This workflow uses molecular surface and residue-layer information to characterize structural exposure and residue distributions.

The workflow includes:

- MolLayers-based layer generation
- Surface residue extraction
- Surface analysis
- Summary analysis of layer-dependent residue properties

MolLayers is an external/local software component and is not bundled with this repository.

### G6P Ligand Preparation

Located in:

```text
scripts/ligand_preparation/
```

The G6P parameterization workflow uses AmberTools-based molecular mechanics parameterization.

The workflow includes ligand preparation and generation of parameters suitable for molecular dynamics simulations.

### Molecular Dynamics Simulation

Located in:

```text
scripts/simulation/
```

The simulation workflow includes:

- apo hexokinase system preparation/simulation
- G6P-bound system preparation/simulation
- mutant simulations
- selected extended simulations

The `extended_simulation.py` script was developed for extending selected production simulations and retains the environment-specific paths used in the original computational workflow.

Large simulation trajectories and simulation output files are intentionally excluded from this repository.

### Trajectory Analysis

Located in:

```text
scripts/md_analysis/
```

The Python workflow includes:

- RMSD
- RMSF
- radius of gyration

VMD scripts are provided under:

```text
scripts/md_analysis/vmd/
```

for trajectory and structural analysis workflows performed using VMD.

### Secondary-Structure Analysis

Located in:

```text
scripts/dssp/
```

This workflow performs secondary-structure analysis from molecular dynamics trajectories.

### Ligand-Interaction Analysis

Located in:

```text
scripts/ligand_interactions/
```

The scripts analyze interactions between hexokinase and G6P, including:

- hydrogen bonds
- contacts
- center-of-mass relationships

### Binding-Energy Analysis

Located in:

```text
scripts/binding_energy/
```

These scripts process and analyze molecular-mechanics-based energy terms and binding-energy calculations generated during the study.

### Dimensionality Reduction

Located in:

```text
scripts/dimensionality_reduction/
```

The workflow uses:

- t-SNE
- UMAP

to integrate and visualize multidimensional molecular-dynamics-derived descriptors.

### Sequence and Evolutionary Analysis

Located in:

```text
scripts/sequence/
```

The repository contains scripts used for visualization of sequence/evolutionary analysis results, including phylogenetic-tree visualization.

---

## Data and File Organization

The repository does **not** contain the molecular dynamics datasets used in the study.

The external data directory is specified through the environment variable:

```bash
HEXOKINASE_DATA_ROOT
```

For example:

```bash
export HEXOKINASE_DATA_ROOT=/media/simulations/hexokinase
```

The analysis scripts expect the external study data to contain the relevant directories, including:

```text
structures/
G6P_structure/
modelled_structures/
ligand_regulatory_site/
```

These directories remain outside the Git repository because they contain structural files, simulation systems, trajectories, and/or generated computational data.

The repository therefore separates:

```text
GitHub repository
    ├── source code
    ├── configuration
    ├── documentation
    └── environment specifications

External study data
    ├── structures
    ├── modelled structures
    ├── G6P systems
    ├── simulation trajectories
    └── generated computational outputs
```

---

## Configuration

The main configuration file is:

```text
config/systems.yaml
```

The repository uses the environment variable:

```bash
HEXOKINASE_DATA_ROOT
```

to locate the external computational data.

Before running analysis scripts, set:

```bash
export HEXOKINASE_DATA_ROOT=/path/to/hexokinase_data
```

For the original working environment, the data root was:

```bash
export HEXOKINASE_DATA_ROOT=/media/simulations/hexokinase
```

For other users, replace this with the location of their own study-data directory.

The path configuration is handled centrally by:

```text
scripts/utils/paths.py
```

This avoids requiring analysis scripts to contain the locations of the external data directories.

---

## Software Environment

The repository separates Python dependencies from specialized molecular-dynamics and structural-biology software.

### Python environment

The Python analysis environment is specified in:

```text
environment/requirements.txt
```

and:

```text
environment/environment.yml
```

The environment includes packages used for:

- molecular-dynamics trajectory analysis
- numerical analysis
- statistical analysis
- dimensionality reduction
- plotting and visualization
- sequence analysis
- data processing

### External software

Specialized software used in the workflow includes:

- ACEMD
- HTMD
- AmberTools
- MoleculeKit
- MODELLER
- NACCESS
- MolLayers
- VMD

Installation and version information for these tools is documented in:

```text
environment/external_software.md
```

These specialized programs are not distributed with this repository.

---

## Reproducibility

The repository is designed to preserve the computational workflow while keeping large simulation data outside version control.

To reproduce the analysis:

1. Clone the repository.
2. Install the required Python environment.
3. Install the required external molecular-dynamics and structural-analysis software.
4. Obtain or prepare the corresponding input structures and simulation data.
5. Set `HEXOKINASE_DATA_ROOT` to the external data directory.
6. Review the workflow documentation.
7. Run the relevant scripts in the documented workflow order.

Detailed workflow information is available in:

```text
docs/workflow.md
```

Script-specific information is available in:

```text
docs/script_guide.md
```

---

## Important Reproducibility Note

Some scripts are intended to be run within specific software environments or from specific working directories because they were developed as part of the original computational workflow.

In particular, simulation and specialized structural-analysis scripts may depend on:

- ACEMD/HTMD installation
- MODELLER installation
- MolLayers availability
- VMD
- external molecular structures
- simulation topology and trajectory files
- environment-specific executable locations

Therefore, this repository should be regarded as the **version-controlled source code and workflow record for the study**, rather than a self-contained simulation dataset.

---

## Excluded Files

Large or generated computational files are intentionally excluded from Git.

Examples include:

```text
*.xtc
*.dcd
*.trr
*.coor
*.vel
*.chk
*.restart
*.prmtop
*.psf
*.tpr
*.gro
*.top
```

as well as generated analysis outputs and temporary files.

This keeps the repository lightweight and focused on the reproducible computational workflow.

---

## Documentation

Additional documentation is provided in:

```text
docs/
├── script_guide.md
└── workflow.md
```

### `script_guide.md`

Provides script-level information including:

- purpose
- inputs
- outputs
- dependencies
- execution location
- usage

### `workflow.md`

Describes the overall computational workflow and the relationship between the different analysis stages.

---

## Citation

If you use the computational workflows or scripts from this repository, please cite the associated research publication:

**Kiruthiga Shankar Kumar and Naga Bhushana Rao Karampudi.**

*Understanding the structural and functional dynamics of Hexokinase Isoforms in Cancer Metabolism using Molecular Dynamics simulation studies.*

Publication and DOI information will be added here when available.

---

## License

This repository is distributed under the license specified in:

```text
LICENSE
```

Third-party software used by the workflow is subject to its own licensing terms.

---

## Contact

For questions regarding the computational workflow, analysis scripts, or scientific methodology, please refer to the associated research publication and project documentation.
