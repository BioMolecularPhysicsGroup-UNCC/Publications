#!/bin/bash
maxprot=5365 
minlig=1
maxlig=47
molname=GID-GA4minus
force=20

topfile=../ref.prmtop
rstfile=../ref-min.crd
bincoor=../md_restart.coor
binvel=../md_restart.vel
xscfile=../md_restart.xsc

## Generate a random integer for the initial direction of the force
## The integer is between 0 - 32767

number=$1

outdir=${molname}_${number}_ramd_force_${force}_out
if [ -d "$outdir" ]; then rm -rf $outdir; mkdir $outdir; else mkdir $outdir; fi

cd $outdir
echo -e "#!/bin/bash\nmodule load slurm\nscancel $SLURM_JOB_ID" > qdel_jobid.sh
chmod +x qdel_jobid.sh
cat << EOF > ${molname}_${number}_ramd_force_${force}_rMin_025.namdin

# *** AMBER force field ********************************************************
amber                                         on
parmfile                                   $topfile
ambercoor                                  $rstfile
bincoordinates                             $bincoor
readexclusions                               yes
exclude                                scaled1-4
1-4scaling                                     0.83333333   #=1/1.2
scnb                                           2

#*** approximations for nonbonded interactions***********************************
cutoff                                        10
switching                                     on
switchdist                                    8.0
pairlistdist                                  36
outputpairlists                             1000
stepspercycle                                 20
nonbondedfreq                                  1
fullelectfrequency                             2
pairlistsPerCycle   2
twoawayx            yes
ldbUnloadOne      yes
noPatchesOnOne      yes

#***timestep*********************************************************************
minimization                                 off
numsteps                                    10000000
timestep                                       2

#***SHAKE use
rigidbonds                                   all
rigidTolerance                             1e-08

#***Basic Dynamics***************************************************************
zeroMomentum                                  no

#***temperature control**********************************************************
binvelocities                            $binvel
langevin                                      on
langevintemp                                 300.0
langevindamping                                1.0
langevinhydrogen                             off

#***initial contraints***********************************************************
constraints                                  off

#***initial fixed atoms***********************************************************
fixedAtoms                                   off

#***pressure control*************************************************************
useGroupPressure                             yes
useFlexibleCell                               no
useConstantArea                               no
LangevinPiston                                on
LangevinPistonTarget                           1.01325
LangevinPistonPeriod                        1000
LangevinPistonDecay                         1000
LangevinPistonTemp                           300
#SurfaceTensionTarget                          60

#***PME and PBC******************************************************************
PME                                           on
PMETolerance                               1e-06
PMEGridSpacing                                 1
extendedSystem                             $xscfile
wrapAll                                      off
XSTfile                               ${molname}_eq04_${number}_ramd_1.xst
XSTfreq                                     1000

#***Interactive Molecular Dynamics***********************************************
IMDon                                        off

#***output***********************************************************************
outputname                             ${molname}_${number}_ramd_0${force}
outputenergies                              500
outputtiming                                1000
restartname                            ${molname}_${number}_ramd_0${force}.rst
restartfreq                                 2000
dcdfile                                ${molname}_${number}_ramd_0${force}.dcd
dcdfreq                                     500
veldcdfile                             ${molname}_${number}_ramd_0${force}.vcd
veldcdfreq                                  2000
binaryoutput                                 off
binaryrestart                                 on

#*** Random Acceleration Molecular Dynamics *************************************
source /nesi/project/nesi02703/Software/Scripts/ramd-5.0.5/ramd-5.0.5/scripts/ramd-5.tcl
ramd namdVersion 2.12
ramd ramdfilename "someRamd.log"
#source /nesi/project/nesi02703/Software/Scripts/ramd/scripts/ramd-4.1.tcl
#source /apps/pkg/namd-2.11/rhel7_u2-x86_64-mpi/gnu/namd/lib/ramd/scripts/ramd-4.1.tcl
#*** sources the wrapper script ramd-4.1.tcl;
#*** please change the directory '../scripts/' to '$dir' ( the correct path );
#*** directory '$dir' should contain the scripts: ramd-4.1.tcl, ramd-4.1_script.tcl, and vectors.tcl

ramd debugLevel                       0
#*** activates verbose output if set to something else than 0

ramd ramdSteps                       50
#*** specifies the number of steps in 1 ramd stint;
#*** defaults to 50

ramd forceRAMD                          ${force}
#*** specifies the acceleration to be applied;
#*** defaults to 0.25 kcal/mol*A*amu

ramd rMinRamd                         0.025
#*** specifies the minimum distance to be travelled by the ligand in 1 ramd stint;
#*** defaults to 0.01 Angstr

ramd forceOutFreq                    10
#*** every 'forceOutFreq' steps detailed output of forces will be written;
#*** defaults to 0 (no detailed output)

ramd maxDist                        30
#*** specifies the distance between the COMs of the ligand and the protein when the simulation is stopped
#*** defaults to 50 Angstr

ramd firstProtAtom                    1
#*** specifies the index of the first protein atom
#*** defaults to 1 (assumes first atom in the system corresponds to first protein atom

ramd lastProtAtom                  ${maxprot}
#*** specifies the index of the last protein atom
#*** required; simulation exits if this parameter is not set

ramd firstRamdAtom                 ${minlig}
#*** specifies the index of the first ligand atom
#*** required; simulation exits if this parameter is not set

ramd lastRamdAtom                  ${maxlig}
#*** specifies the index of the last ligand atom
#*** required; simulation exits if this parameter is not set

ramd ramdSeed                     ${number}
#*** specifies the seed for the random number generator (for the generation of acceleration directions)
#*** defaults to 14253
#*** please change if you wish to run different trajectories
EOF
####################################################################################

srun namd2 ${molname}_${number}_ramd_force_${force}_rMin_025.namdin +isomalloc_sync  >  ${molname}_${number}_ramd_force_${force}_rMin_025.out 

#cm-launcher
#time mpirun -np $NSLOTS  $NAMD_HOME/namd2 \
#./qdel_jobid.sh
#exit
