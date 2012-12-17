#!/bin/bash
if [ "$#" -eq 0 ];then
   echo "run.sh pdb_file"
   exit
fi

rosetta_home=${HOME}/rosetta-trunk
mode=release
path_name=$(dirname $1)
pdb_name=$(basename $1 .pdb)
first_res=`grep "ATOM.................B" $1|head -1|awk '{print $6}'`
last_res=`grep "ATOM.................B" $1|tail -1|awk '{print $6}'`
db=${rosetta_home}/rosetta_database
rosetta_scripts=${rosetta_home}/rosetta_source/bin/rosetta_scripts.linuxgcc${mode}

${rosetta_scripts} -database ${db}  -s ${path_name}/${pdb_name}.pdb  -parser:protocol ppk_minbb_use_correct.xml -ex1 -ex2 -use_input_sc -correct -overwrite

${rosetta_scripts} -dadtabase ${db}  -s ${path_name}/${pdb_name}_0001.pdb  -parser:protocol ppk_minbb_add_elec.xml -ex1 -ex2 -use_input_sc -correct -parse_charge -overwrite

${rosetta_scripts} -database ${db}  -s ${path_name}/${pdb_name}_0001_0001.pdb -parser:protocol scan_binding_with_elec.xml   -ex1 -ex2 -use_input_sc -nooutput -correct -parse_charge -overwrite
