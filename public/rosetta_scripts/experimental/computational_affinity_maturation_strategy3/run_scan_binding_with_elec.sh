#!/bin/bash
#if [ "$#" -eq 0 ];then
#   echo "run.sh pdb_file"
#   exit
#fi

rosetta_home=${HOME}/dvlp/rosetta-trunk
mode=release
#pdb_path=$1
pdb_path=~/pdbs/HB80.3_relax.pdb
echo input structure=${pdb_path}
path_name=`dirname ${pdb_path}`
echo dir=${path_name}
pdb_name=`basename ${pdb_path} .pdb`
echo filename=${pdb_name}
db=${rosetta_home}/rosetta_database
echo db=${db}
rosetta_scripts=${rosetta_home}/rosetta_source/bin/rosetta_scripts.apbs.linuxgcc${mode}
echo executable=${rosetta_scripts}

export LD_LIBRARY_PATH=${rosetta_home}/rosetta_source/external/apbs/apbs-1.4-rosetta/lib:$LD_LIBRARY_PATH

#${rosetta_scripts} -database ${db}  -s ${path_name}/${pdb_name}.pdb  -parser:protocol ppk_minbb_use_correct.xml -ex1 -ex2 -use_input_sc -correct -overwrite

#${rosetta_scripts} -database ${db}  -s ${pdb_name}_0001.pdb  -parser:protocol ppk_minbb_add_elec.xml -ex1 -ex2 -use_input_sc -correct -parse_charge -overwrite

${rosetta_scripts} -database ${db}  -s ${path_name}/${pdb_name}.pdb -parser:protocol scan_binding_with_elec.xml   -ex1 -ex2 -use_input_sc -nooutput -correct -parse_charge -overwrite
