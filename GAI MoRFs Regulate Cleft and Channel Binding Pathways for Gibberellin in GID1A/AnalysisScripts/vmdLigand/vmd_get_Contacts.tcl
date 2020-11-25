# name of an input dcd file  
set dcd_name  [lindex $argv 0]
set out_name  [lindex $argv 1]
set Hdist [lindex $argv 2]
set step [lindex $argv 3]
# name of ligand residues in PDB file
set ligand  UNL

set name $dcd_name
mol load  parm7 ./ref.prmtop dcd $name
set mol1 [molinfo top get id]
put "END INI"

package require pbctools
pbc set [pbc get -all] -all
pbc wrap -center com -centersel protein -compound res -all

# use use frame 0 for the refereince
set reference [atomselect $mol1 "protein and noh" frame 0]
set reflig [atomselect $mol1 "resname $ligand and noh"  frame 0]

#compute number of atoms in the ligand 
set atomnumber  0
  foreach a [$reflig get name]  {
  set atomnumber [expr $atomnumber+ 1]
 }
puts "Atoms in ligand: $atomnumber"

# the frame being compared
set selProt [atomselect $mol1 "hydrophobic"]
#would like to have a logic check for all other, but didnt work
set n [molinfo $mol1  get numframes]
# all ligand atoms
set ligatoms [atomselect $mol1 "resname $ligand"]

#output files
set outfile [open $out_name-lig_H0bond_$Hdist.dat w]

# loop over frames
for { set frame 0 } {$frame < $n } { incr frame $step} {
 # get correct frame
 $selProt frame $frame
 $ligatoms frame $frame

 set listProt [measure contacts $Hdist $selProt $ligatoms]
 set atomsProt [lindex $listProt 0]
 set NumProt [llength $atomsProt]
 set ProtID {}
 set ResNum {}
 foreach i [lindex $atomsProt] {
      	lappend ProtID [[atomselect top "index $i" ] get resname ] 
	lappend ResNum [[atomselect top "index $i" ] get resid ] 
 }

 puts $outfile "$frame	 $NumProt $ProtID $ResNum"
}
# end for loop
close $outfile

set scale [expr $n-1]
puts $scale

quit
