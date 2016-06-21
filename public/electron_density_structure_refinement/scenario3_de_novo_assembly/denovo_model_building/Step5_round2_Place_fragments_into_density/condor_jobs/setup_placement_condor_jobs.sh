ln -s ../../Step4_remove_density_using_average_model/residues_to_search_for_2nd_round.txt .
../../../scripts/setup_placements_condor_jobs.py \
    --target_list  residues_to_search_for_2nd_round.txt \
    --mapfile     ../../Step4_remove_density_using_average_model/round2_r2.mrc \
    --fasta       ../../../input_files/trpv1.fasta \
    --fragfile    ../../../input_files/t001_.200.9mers  \
    --condor_job_fn placement_condor_job \
    --submit_script_fn submit_placement_condor_jobs.sh
