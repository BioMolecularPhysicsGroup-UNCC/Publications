#!/bin/sh

#module load ambertools
#module load anaconda3

cp ../cleaned_DELLA_GID.pdb .
cp ../waters2.pdb .
cp ../tleap_all_in .
cp ../amber_prep.sh .
cp ../run_AmberMin.pbs .
cp ../equilRAMD.pbs .

name="GA4"
echo "${name}"
grep ATOM cleaned_DELLA_GID.pdb > prot_lig.pdb
echo "TER" >> prot_lig.pdb
grep ATOM "${name}"minus.pdb >> prot_lig.pdb
echo "TER" >> prot_lig.pdb
grep HOH  waters2.pdb  >>  prot_lig.pdb
sed 's/GA3/'"${name}"'/g' tleap_all_in > tleap_all_in_2
#tleap -f tleap_all_in_2
#qsub run_AmberMin.pbs

#old stuff from ligand parameterization
#cp ../GA9minus/tleap_ligand_in .
#cp ../GA9minus/acpype.py .
#name="GA8"
#echo "${name}"
#antechamber -i "${name}"minus.pdb -fi pdb -o "${name}"minus.mol2 -fo mol2 -c bcc -s 2 -nc -1
#parmchk2 -i "${name}"minus.mol2 -f mol2 -o "${name}"minus.frcmod
#sed 's/GA9/'"${name}"'/g' tleap_ligand_in > tleap_ligand_in_2
#tleap -f tleap_ligand_in_2
#python3 acpype.py -i "${name}minus.mol2"
