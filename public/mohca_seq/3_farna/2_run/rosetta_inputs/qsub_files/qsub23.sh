#!/bin/bash
#PBS -N _Users_clarencecheng_Rosetta_demos_public_mohca_seq_3_farna_2_run_rosetta_inputs_out_23
#PBS -o /dev/null
#PBS -e /dev/null
#PBS -m n
#PBS -M nobody@stanford.edu
#PBS -l walltime=16:00:00

cd /Users/clarencecheng/Rosetta/demos/public/mohca_seq/3_farna/2_run/rosetta_inputs

/Users/clarencecheng/Rosetta//main/source/bin/rna_denovo -nstruct 500 -params_file glycine_riboswitch.params -fasta glycine_riboswitch.fasta -out:file:silent glycine_riboswitch.out -include_neighbor_base_stacks -minimize_rna false -native glycine_riboswitch_3p49_native_RNA.pdb -in:file:silent helix0.out helix1.out helix2.out helix3.out helix4.out helix5.out helix6.out helix7.out -input_res 2-9 65-72 16-21 26-31 33-35 54-56 39-42 48-51 81-85 155-159 92-97 101-106 108-110 145-147 114-117 139-142 -cst_file glycine_riboswitch_constraints -staged_constraints -cycles 20000 -ignore_zero_occupancy false -output_res_num 1-159 > /dev/null 2> /dev/null 
