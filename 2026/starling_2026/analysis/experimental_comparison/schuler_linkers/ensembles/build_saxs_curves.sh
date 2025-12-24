#!/bin/zsh

for i in dCh_minus  dTRBP  sNh_plus  sNrich
do
    mkdir $i
    autoSCTR --pdb ${i}_STARLING.pdb --xtc ${i}_STARLING.xtc -n 400

    mv average_* $i/
    mv logfile_autoSCRT.txt $i/
    
done
