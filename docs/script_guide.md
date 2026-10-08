# Hexokinase Script Guide

This guide documents the retained scripts in the repository. The
documentation is derived from the script source and records the
intended inputs, outputs, execution context, and external dependencies.
It does not include molecular-dynamics trajectories or generated result
files.

## Repository-wide requirement

Analysis scripts that import `scripts.utils.paths` require:

```bash
export HEXOKINASE_DATA_ROOT=/path/to/hexokinase_data
```

The data root should contain the directories configured in
`config/systems.yaml`.

## Workflow order

1. **Structure modeling**
   - `scripts/structure/modeller/`
2. **Mutation modeling**
   - `scripts/structure/mutations/`
3. **Surface/layer analysis**
   - `scripts/surface_analysis/run_layers.py`
   - `scripts/surface_analysis/extract_surface_residue_layers.py`
   - `scripts/surface_analysis/hk_layers_analysis.py`
4. **G6P/system preparation and MD**
   - `scripts/ligand_preparation/`
   - `scripts/simulation/`
   - `scripts/structure/mutations/mutation_g6p.py`
5. **Trajectory analysis**
   - `scripts/md_analysis/rmsd_rmsf_rg.py`
   - `scripts/md_analysis/vmd/`
   - `scripts/dssp/dssp_hk.py`
6. **Ligand-interaction analysis**
   - `scripts/ligand_interactions/`
7. **Binding-energy analysis**
   - `scripts/binding_energy/`
8. **Dimensionality reduction**
   - `scripts/dimensionality_reduction/`
9. **Manuscript visualization**
   - `scripts/visualization/`
10. **Sequence/phylogenetic visualization**
   - `scripts/sequence/visualize_tree.py`

## Important execution distinction

Some scripts are repository-root workflows and use
`scripts.utils.paths`; others intentionally operate from a directory
containing their input files (for example the MODELLER scripts and VMD
Tcl scripts). Each script now states its own execution context in its
module documentation.
