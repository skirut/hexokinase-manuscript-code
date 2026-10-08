# VMD per-residue RMSF analysis script.
#
# Purpose
# -------
# Calculates residue-level RMSF for the protein selection over the
# configured trajectory frames and writes the resulting values to a data
# file.
#
# Inputs
# ------
# - The molecule already loaded in VMD.
# - Protein atom selection and frame range defined in the script.
#
# Outputs
# -------
# - ``rmsf.dat`` in the current working directory.
#
# Run from
# --------
# Run inside VMD after loading the desired structure and trajectory.
#
# Usage
# -----
# ``source scripts/md_analysis/vmd/rmsf.tcl``
#
# Dependencies
# ------------
# VMD Tcl scripting and the VMD atom-selection/measure commands.

set start_frame first
set end_frame last
set sel "protein"
set outfile [open rmsf.dat w] 
set select [atomselect top "$sel"]    
set nf [molinfo top get numframes] 
set reslist [lsort -integer -unique [$select get resid]]
if {$start_frame == "first"} {
        set start_frame 1
}
if {$end_frame == "last"} {
        set end_frame [expr $nf - 1]
}
foreach res $reslist {
    set sel_tmp [atomselect top "$sel and resid $res"]
    set atomlist [$sel_tmp get index]
    set natoms [expr [llength $atomlist]]
    set rmsf_list [measure rmsf $sel_tmp first $start_frame last $end_frame]
    set rmsf_tmp 0
    for { set j 0 } { $j <= [expr $natoms - 1] } { incr j } { 
        set rmsf_tmp [expr $rmsf_tmp + [expr [lindex $rmsf_list $j]] ** 2 ]
    }
    set rmsf_tmp [expr ($rmsf_tmp / $natoms) ** 0.5 ]
    puts -nonewline $outfile "$res "
    puts $outfile "$rmsf_tmp"
    unset rmsf_list
}
close $outfile
puts "All Done!"

mol new structure.pdb
mol delrep 0 top
mol representation NewCartoon
mol color Structure
mol selection {not water}
mol material Opaque
mol addrep top
mol addfile output.xtc
set start_frame first
set end_frame last
set sel "protein"
set outfile [open rmsf.dat w]
set select [atomselect top "$sel"]
set nf [molinfo top get numframes]
set reslist [lsort -integer -unique [$select get resid]]
if {$start_frame == "first"} {
    set start_frame 0
}
if {$end_frame == "last"} {
    set end_frame [expr $nf - 1]
}
foreach res $reslist {
    set sel_tmp [atomselect top "$sel and resid $res"]
    set atomlist [$sel_tmp get index]
    set natoms [llength $atomlist]
    set rmsf_list [measure rmsf $sel_tmp first $start_frame last $end_frame]
    set rmsf_tmp 0
    for { set j 0 } { $j < $natoms } { incr j } {
        set rmsf_tmp [expr $rmsf_tmp + [lindex $rmsf_list $j] ** 2]
    }
    set rmsf_tmp [expr sqrt($rmsf_tmp / $natoms)]
    puts $outfile "$res $rmsf_tmp"
    
    $sel_tmp delete
    unset rmsf_list
}
$select delete
close $outfile
puts "All Done!"


