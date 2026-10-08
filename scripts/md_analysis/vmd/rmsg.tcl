# VMD radius-of-gyration analysis script.
#
# Purpose
# -------
# Calculates the radius of gyration of the protein selection for every
# trajectory frame.
#
# Inputs
# ------
# - The molecule and trajectory already loaded in VMD.
# - Protein atom selection defined in the script.
#
# Outputs
# -------
# - ``rmsg.dat`` in the current working directory.
#
# Run from
# --------
# Run inside VMD after loading the desired structure and trajectory.
#
# Usage
# -----
# ``source scripts/md_analysis/vmd/rmsg.tcl``
#
# Dependencies
# ------------
# VMD Tcl scripting and the VMD ``measure rgyr`` command.

set sel "protein"
set outfile [open rmsg.dat w]
set nf [molinfo top get numframes]
set select [atomselect top "$sel"]
for {set i 0} {$i < $nf} {incr i} {
    $select frame $i
    puts $outfile "[measure rgyr $select]"
}
close $outfile
puts "All Done!"


