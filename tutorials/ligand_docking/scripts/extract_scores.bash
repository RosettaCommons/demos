#!/bin/bash

#extract_scores.bash score.sc

awk 'BEGIN{OFS=",";}{
	
		if( NR !=1){
		if( $2 == "total_score" ){
			for(x=1; x<= NF; ++x){
				if($x=="ligand_rms_no_super_X"){
					rms=x
				}
				if($x == "interface_delta_X"){
					score=x
				}
				if($x == "total_score"){
					total=x
				}
				if($x == "description"){
					file=x
				}
		} 
		}else {
		print $file, $total, $score, $rms
			}
		}
	}' $1

