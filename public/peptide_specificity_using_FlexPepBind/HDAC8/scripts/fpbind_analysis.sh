#! /bin/bash

	for i in $(cat input_files/peptide.list)
	do

		cd $i/Minimization/

		# @extracting different scoring terms
		for j in total_score I_sc pep_sc_noref reweighted_sc
			do ../../scripts/printScoreFile_byHeader.pl min.score.sc $j | awk '{print $2}' >${j}_1
		done
		cd ../..

	done


	if [ -d score_analysis ]; then
		rm -r score_analysis;
	fi
	mkdir score_analysis;

	for i in $(cat input_files/peptide.list)
	do
		for j in total_score I_sc pep_sc_noref reweighted_sc;do
			sed 1d $i/Minimization/${j}_1 | sort -nk 1 | head -1 >> score_analysis/${j}
		done	
 	done



exit 0;

