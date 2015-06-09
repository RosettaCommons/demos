#!/bin/bash

PATH_TO_APP=../../doug_dock_design_min_mod2_cal_cal.macosgccrelease
PATH_TO_DB=../../minirosetta_database

for i in pos_*; do
    cd $i
    echo "Working on $i"

    $PATH_TO_APP -database $PATH_TO_DB -s ../../inputs/1NX1_clean_repack_min_all.pdb -resfile ../resfile_$i -nstruct 255 -inner_num 45 -pert_num 25 -ia_ener 100 -use_input_sc -pdb_gz >& $i.log

    echo "Finished $i"
    cd ../
done