# name of an input dcd file  
set dcd_name  [lindex $argv 0]
set out_name  [lindex $argv 1]
# name of ligand residues in PDB file
set ligand  UNL
#  residue number of Lys that is placed is the ligand exis passway
#set residue  [lindex $argv 1]
#puts "========================= $number"
set name $dcd_name
mol load  parm7 ./ref.prmtop dcd $name
set mol1 [molinfo top get id]
put "END INI"

package require pbctools
pbc set [pbc get -all] -all
pbc wrap -center com -centersel protein -compound res -all

# use use frame 0 for the refereince
sest reflig [atomselect $mol1 "resname $ligand and noh"  frame 0]
set reference [atomselect $mol1 "water within 5 of resname $ligand" frame 0]
#compute number of atoms in the ligand 
set atomnumber  0
  foreach a [$reference get name]  {
  set atomnumber [expr $atomnumber+ 1]
 }
puts "Atoms in ligand water: $atomnumber"
# the frame being compared
set sel [atomselect $mol1 "water within 5 of resname $ligand"]
set sel1 [atomselect $mol1 "all"]
set n [molinfo $mol1  get numframes]
#output files
set outfile2 [open $out_name-ligWater_rmsd.dat w]

# loop over frames
for { set frame 0 } {$frame < $n } { incr frame } {
 # get correct frame
 $sel frame $frame
 $ligatoms frame $frame
 $sel1 frame $frame
 #compute the transformation
 set trans_mat [measure fit $sel $reference]
 #do the  alignment
 $sel1 move $trans_mat
 # compute RMSD
# RMSD of ligand 
 set rmsd2 [measure rmsd $sel $reference]
 puts $outfile2 "$frame $rmsd2"
 }
close $outfile2

set scale [expr $n-1]
puts $scale
set sel1 [atomselect $mol1 "water within 5 of resname $ligand" frame $scale]
set sel2 [atomselect $mol1 "water within 5 of resname $ligand" frame 0]
#animate write pdb $outpdbfile $sel1 beg $scale end $scale
$sel1 writepdb $out_name-lastWater.pdb
$sel2 writepdb $out_name-firstWater.pdb
close $outpdbfile
quit

