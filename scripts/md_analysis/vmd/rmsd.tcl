# VMD RMSD analysis script.
#
# Purpose
# -------
# Loads the current Hexokinase structure/trajectory and performs the RMSD
# analysis implemented in the Tcl script.
#
# Inputs
# ------
# - ``structure.pdb``
# - ``output.xtc``
# - The selection and frame settings defined in the script.
#
# Outputs
# -------
# Writes the RMSD data file specified by the Tcl commands in the script.
#
# Run from
# --------
# Run inside VMD from the directory containing the structure and
# trajectory files.
#
# Usage
# -----
# ``source scripts/md_analysis/vmd/rmsd.tcl``
#
# Dependencies
# ------------
# VMD Tcl scripting and the VMD ``measure``/atom-selection commands.

# Load the structure
mol new structure.pdb

# Delete the default representation
mol delrep 0 top

# Create a new representation
mol representation NewCartoon
mol color Structure
mol selection {not water}
mol material Opaque
mol addrep top

# Load the trajectory
mol addfile output.xtc

# Wait until the trajectory is fully loaded
while {[molinfo top get numframes] < 1001} {
    after 100
}

# Set up atom selections
set frame0 [atomselect top "protein and backbone and noh" frame 0]
set sel [atomselect top "protein and backbone and noh"]
set all [atomselect top all]

# Get the number of frames
set nf [molinfo top get numframes]

# Open output file for writing RMSD values
set outfile [open "rmsd_analysis.dat" w]

# Loop through each frame
for {set i 0} {$i < $nf} {incr i} {
    $sel frame $i
    $all frame $i
    $all move [measure fit $sel $frame0]
    set rmsd [measure rmsd $sel $frame0]
    puts $outfile "$i $rmsd"
}

# Close the output file
close $outfile

# Print completion message
puts "RMSD analysis completed and written to rmsd_analysis.dat"


