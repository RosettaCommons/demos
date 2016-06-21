bcl svn revision: 4434

The three cone model parameters can be calculated from multiple mutants each having multiple relaxed structures. The command line given below will calculate the statistics for the pdbs of the two mutants provided in the inputs directory.
bin/bcl.exe CalculateConeModelParameters -list_of_pdb_lists input/pdb_lists.ls -aaclass AAComplete -prefix calc_params_ -sl_cb_sl_max_angle 0 20 9 sl_cb_sl_max_angle.histogram -sleffective_cb_ca_angle 0 20 9 sleffective_cb_ca_angle.histogram -sleffective_cb_distance 0 1 10 sleffective_cb_distance.histogram -use_degrees -parameter_table_filename out_summary.table > & calc_params.log &
The outputted files are histograms of the statistics for the three parameters
as well as a summary file (out_summary.table). 
