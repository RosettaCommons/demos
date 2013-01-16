#!/bin/bash
#if [ "$#" -eq 0 ];then
#   echo "run.sh pdb_file"
#   exit
#fi

rosetta_home=${HOME}/rosetta-trunk
mode=debug
#pdb_path=$1
pdb_path=~/pdbs/HB80.3_relax.pdb
echo input structure=${pdb_path}
path_name=`dirname ${pdb_path}`
echo dir=${path_name}
pdb_name=`basename ${pdb_path} .pdb`
echo filename=${pdb_name}
db=${rosetta_home}/rosetta_database
echo db=${db}
rosetta_scripts=${rosetta_home}/rosetta_source/bin/rosetta_scripts.linuxgcc${mode}
echo executable=${rosetta_scripts}

echo I am here

#${rosetta_scripts} -database ${db}  -s ${path_name}/${pdb_name}.pdb  -parser:protocol ppk_minbb_use_correct.xml -ex1 -ex2 -use_input_sc -correct -overwrite

echo here I am

${rosetta_scripts} -database ${db}  -s ${pdb_name}_0001.pdb  -parser:protocol ppk_minbb_add_elec.xml -ex1 -ex2 -use_input_sc -correct -parse_charge -overwrite

echo here we go
#${rosetta_scripts} -database ${db}  -s ${path_name}/${pdb_name}_0001_0001.pdb -parser:protocol scan_binding_with_elec.xml   -ex1 -ex2 -use_input_sc -nooutput -correct -parse_charge -overwrite
