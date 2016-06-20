cd condor_jobs/
    cat ./*/clustering_selected.sc | grep -v "^#" > candidate_fragment_placements.sc
cd ../

../../scripts/prepare_round2.py \
    --round1_mc_picked_frags_dir ../Step3_Simulated_annealing_Monte_Carlo_sampling/5percent_selected_models/ \
    --round1_clustering_selected_results_file ../Step1_Place_fragments_into_density/condor_jobs/candidate_fragment_placements.sc \
    --round1_frags_dir ../Step1_Place_fragments_into_density/candidate_fragment_placements \
    --round2_clustering_selected_results_file condor_jobs/candidate_fragment_placements.sc \
    --round2_frags_dir candidate_fragment_placements  \
    --total_rsd 315 

