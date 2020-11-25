# name of an input dcd file  
set dcd_name  [lindex $argv 0]
set out_name  [lindex $argv 1]
set Hdist [lindex $argv 2]
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
set reference [atomselect $mol1 "protein and noh" frame 0]
set reflig [atomselect $mol1 "water within 5 of resname $ligand"  frame 0]

#compute number of atoms in the ligand 
set atomnumber  0
  foreach a [$reflig get name]  {
  set atomnumber [expr $atomnumber+ 1]
 }
puts "Atoms in ligand: $atomnumber"

# the frame being compared
set selProt [atomselect $mol1 "protein or resname $ligand"]
#would like to have a logic check for all other, but didnt work
#set selAll [atomselect $mol1 "all"]
set selSolv [atomselect $mol1 "water"]
set n [molinfo $mol1  get numframes]
# all ligand atoms
#set ligatoms [atomselect $mol1 "resname $ligand and noh"]
set ligatoms [atomselect $mol1 "water within 5 of resname $ligand" frame 0]
#output files
set outfile [open $out_name-ligWater_Hbond_$Hdist.dat w]

# loop over frames
for { set frame 0 } {$frame < $n } { incr frame } {
 # get correct frame
 $selProt frame $frame
 $ligatoms frame $frame
 $selSolv frame $frame


 set listProt [measure hbonds $Hdist 30 $selProt $ligatoms]
 set atomsProt [lindex $listProt 0]
 set NumProt [llength $atomsProt]
 set ProtID {}
 set ResNum {}
 foreach i [lindex $atomsProt] {
      	lappend ProtID [[atomselect top "index $i" ] get resname ] 
	lappend ResNum [[atomselect top "index $i" ] get resid ] 
 }

 set listSolv [measure hbonds $Hdist 30 $ligatoms $selSolv]                                                                                         
 set atomsSolv [lindex $listSolv 0]                                                                                                             
 set NumSolv [llength $atomsSolv]
 

 puts $outfile "$frame	 $NumProt $NumSolv	$ProtID $ResNum"
}

close $outfile


set scale [expr $n-1]
puts $scale

#set sel1 [atomselect $mol1 "all" frame $scale]
#set sel2 [atomselect $mol1 "all" frame 0]

#animate write pdb $outpdbfile $sel1 beg $scale end $scale
quit


