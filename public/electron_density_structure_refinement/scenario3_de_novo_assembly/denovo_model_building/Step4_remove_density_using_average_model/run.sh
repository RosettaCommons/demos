# chewing back 4 residues at gap regions
ln -s ../Step3_Simulated_annealing_Monte_Carlo_sampling/average_model/average.pdb
../../scripts/chewing_back_gap_regions.py -p average.pdb -f ../../input_files/trpv1.fasta  > average_to_clean_density_map.pdb 

# using this to maskout density where there are residues being assigned
../../rosetta/densityMap_eraser.static.linuxgccrelease  -pdb_fn average_to_clean_density_map.pdb -mapfile ../../input_files/transmem.mrc  -database ../../rosetta/rosetta_database -radius 2 -outmap_fn round2_r2.mrc 

# setup jobs
../../scripts/print_out_remaining_rsds_given_averagemodel.py --pdb average.pdb -f ../../input_files/trpv1.fasta  > residues_to_search_for_2nd_round.txt

