#!/usr/bin/python
from os import system
import string

command= '~/SWA_RNA_python/SWA_dagman_python/SWA_DAG/SWA_rna_build_dagman.py '

command+= '-s  template.pdb -fasta fasta '

command+= '-native native.pdb '

command+= '-clusterer_num_pose_kept 1000 '

command+= '-force_field_file rna/rna_loop_hires_04092010.wts '

command+= '-sampler_extra_anti_chi_rotamer false '

command+= '-sampler_extra_syn_chi_rotamer false '

command+= '-clusterer_quick_alignment true '

command+= '-clusterer_optimize_memory_usage true '

command+= '-clusterer_keep_pose_in_memory true '

command+= '-allow_combine_DS_regions False '

command+= '-dinucleotide_at_single_element_cc pure_append_prepend '

command+= '-allow_bulge True '

command+= '-input_res  1 2 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 '

command+= '-rmsd_res  3 4 5 6 7 8 '

command+= '-jump_point_pair_list  1-47 '

command+= '-fixed_res  1 2 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 '

command+= '-alignment_res_list  1-47 '

command+= '-optimize_long_loop_mode True '

command+= '-OLLM_chain_closure_only True '

command+= '-rna_torsion_potential_folder ps_03242010/ '

command+= '-sample_virt_ribose_in_sep_DAG True '

command+= '>LOG_SWA_rna_build_dagman.out '

system(command) 

