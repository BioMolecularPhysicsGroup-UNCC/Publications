#!/bin/bash


echo "
# Minimization performed for relaxing the solute : restraint_wt=500
 &cntrl
  imin=1, ntx=1, irest=0, ntrx=1, ntxo=1,
  ntpr=25, ntwr=500, ntwx=0, ntwv=0, ntwe=0,
  ntf=1, ntb=1, dielc=1.0, cut=10.0,
  nsnb=25, igb=0,
  ntr=1, restraint_wt=500, restraintmask='(:1-236 & !@H= & !@Na= & !@Cl=)',
  maxcyc=1500, ntmin=1, ncyc=500, dx0=0.01, drms=0.0001,
  ntc=1,
 /
" >> ref_min-500
sander -O -i ref_min-500 -p ref.prmtop -c ref.inpcrd -ref ref.inpcrd -o ref-min-500.ou -r ref-min-500.crd -inf ref-min-500.log

echo "
# Minimization performed for relaxing the solute : restraint_wt=1
 &cntrl
  imin=1, ntx=1, irest=0, ntrx=1, ntxo=1,
  ntpr=25, ntwr=500, ntwx=0, ntwv=0, ntwe=0,
  ntf=1, ntb=1, dielc=1.0, cut=10.0,
  nsnb=25, igb=0,
  ntr=1, restraint_wt=100, restraintmask='(:1-236 & !@H= & !@Na= & !@Cl=)',
  maxcyc=1500, ntmin=1, ncyc=100, dx0=0.01, drms=0.0001,
  ntc=1,
 /
" >> ref_min-100
sander -O -i ref_min-100 -p ref.prmtop -c ref-min-500.crd -ref ref-min-500.crd -o ref-min-100.ou -r ref-min-100.crd -inf ref-min-100.log

echo "
# Minimization performed for relaxing the solute : restraint_wt=1
 &cntrl
  imin=1, ntx=1, irest=0, ntrx=1, ntxo=1,
  ntpr=25, ntwr=500, ntwx=0, ntwv=0, ntwe=0,
  ntf=1, ntb=1, dielc=1.0, cut=10.0,
  nsnb=25, igb=0,
  ntr=1, restraint_wt=5, restraintmask='(:1-236 & !@H= & !@Na= & !@Cl=)',
  maxcyc=1500, ntmin=1, ncyc=500, dx0=0.01, drms=5,
  ntc=1,
 /
" >> ref_min-5
sander -O -i ref_min-5 -p ref.prmtop -c ref-min-100.crd -ref ref-min-100.crd -o ref-min-5.ou -r ref-min-5.crd -inf ref-min-5.log

echo "
# Minimization performed for relaxing the solute  : no positional restraints
&cntrl
  imin=1, ntx=1, irest=0, ntrx=1, ntxo=1,
  ntpr=25, ntwr=500, ntwx=0, ntwv=0, ntwe=0,
  ntf=1, ntb=1, dielc=1.0, cut=10.0,
  nsnb=25, igb=0,
  ntr=0,
  maxcyc=1500, ntmin=1, ncyc=500, dx0=0.01, drms=0.0001,
  ntc=1, 
 /
" >> ref_min
sander -O -i ref_min -p ref.prmtop -c ref-min-5.crd -ref ref-min-5.crd -o ref-min.ou -r ref-min.crd -inf ref-min.log

