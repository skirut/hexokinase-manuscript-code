# VMD visualization setup for a Hexokinase trajectory.
#
# Purpose
# -------
# Loads ``structure.pdb`` and ``output.xtc`` into VMD and creates a
# cartoon representation colored by structure while excluding water.
#
# Inputs
# ------
# - ``structure.pdb``
# - ``output.xtc``
# Both are expected in the current working directory.
#
# Outputs
# -------
# No file is written by the script. It configures the active VMD molecule
# and representation.
#
# Run from
# --------
# Run inside VMD from the directory containing the input PDB and
# trajectory.
#
# Usage
# -----
# ``source scripts/md_analysis/vmd/load_molecule.tcl``
#
# Dependencies
# ------------
# VMD with the molecular structure and trajectory files available.

mol new structure.pdb
mol delrep 0 top
mol representation NewCartoon
mol color Structure
mol selection {not water}
mol material Opaque
mol addrep top
mol addfile output.xtc

