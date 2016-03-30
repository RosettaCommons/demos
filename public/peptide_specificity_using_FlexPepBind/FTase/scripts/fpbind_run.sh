#! /bin/bash

ROSETTA_BIN="rosetta/main/source/bin"
ROSETTA_DB="rosetta/main/database"

	for i in $(cat input_files/peptide.list)
	do

		if [ -d $i/ ] ;then
		rm -r $i/
		fi
		
		mkdir $i/
		cd  $i/

		#creating the resfile for the peptide
		awk '{print $1, $2, $3}' ../input_files/resfile >tmp1
		awk '{print $5, $6, $7, $8, $9, $10, $11, $12}' ../input_files/resfile >tmp3
		echo "  $i" | fold -w1 >tmp2
		paste tmp1 tmp2 tmp3 >resfile.$i
		rm tmp1 tmp2 tmp3 

		#generating initial structure from the template complex
		$ROSETTA_BIN/fixbb.linuxgccrelease -database $ROSETTA_DB -s ../input_files/template.pdb -resfile resfile.$i -ex1 -ex2 -ex3 -ex4 -use_input_sc -ndruns 100 -scorefile design.score.sc -extra_res_fa ../input_files/fpp.params -nstruct 1 >design.log
		mv template_0001.pdb $i.pdb
		
    #running minimization 
		mkdir minimization;
		cd minimization;
		ln -s ../../input_files/minimization_flags .
    ln -s ../../input_files/constraints.cst .
		ln -s ../$i.pdb ./start.ppk.pdb
    $ROSETTA_BIN/FlexPepDocking.linuxgccrelease -database $ROSETTA_DB  @minimization_flags -extra_res_fa ../../input_files/fpp.params > minimization.log
	
		cd ../..
	done;

exit 0;
