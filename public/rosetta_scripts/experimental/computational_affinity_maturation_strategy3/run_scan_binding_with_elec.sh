#!/bin/bash
if [ "$#" -eq 0 ];then
   echo "run.sh pdb_file"
   exit
fi

path_name=$(dirname $1)
pdb_name=$(basename $1 .pdb)
first_res=`grep "ATOM.................B" $1|head -1|awk '{print $6}'`
last_res=`grep "ATOM.................B" $1|tail -1|awk '{print $6}'`
DB="~/minirosetta_database"

rosetta_scripts.linuxiccrelease -database ${DB}  -s ${path_name}/${pdb_name}.pdb  -parser:protocol ppk_minbb_use_correct.xml -ex1 -ex2 -use_input_sc -correct -overwrite

/work/yfsong/SVN/fresh/mini/bin/PB_potential.linuxgccrelease -database ${DB}  -s ${path_name}/${pdb_name}_0001.pdb -correct -parse_charge -nooutput 

/work/yfsong/dist/pymol/freemol/bin/apbs.exe ${path_name}/${pdb_name}_0001_0001.in

rosetta_scripts.linuxiccrelease -database ${DB}  -s ${path_name}/${pdb_name}_0001.pdb  -parser:protocol ppk_minbb_add_elec.xml -ex1 -ex2 -use_input_sc -correct -parse_charge -PB_potential_file ${path_name}/${pdb_name}_0001_0001.dx -PB_score_residue_range ${first_res} ${last_res} -score:weights score12_w_corrections.wts -score:patch add_PB_elec.wts_patch -overwrite

PB_potential.linuxgccrelease -database ${DB}  -s ${path_name}/${pdb_name}_0001_0001.pdb -correct -parse_charge -nooutput

apbs.exe ${path_name}/${pdb_name}_0001_0001_0001.in

rosetta_scripts.linuxiccrelease -database ${DB}  -s ${path_name}/${pdb_name}_0001_0001.pdb -parser:protocol scan_binding_with_elec.xml   -ex1 -ex2 -use_input_sc -nooutput -correct -parse_charge -PB_potential_file ${path_name}/${pdb_name}_0001_0001_0001.dx -PB_score_residue_range ${first_res} ${last_res} -score:weights score12_w_corrections.wts -score:patch add_PB_elec.wts_patch -overwrite
