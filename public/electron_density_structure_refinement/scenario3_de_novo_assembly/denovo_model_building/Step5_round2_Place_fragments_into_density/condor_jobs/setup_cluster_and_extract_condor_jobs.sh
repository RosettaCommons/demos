../../../scripts/setup_condor_jobs.py \
    --target_list  residues_to_search_for_2nd_round.txt \
    --fragfile ../../../input_files/t001_.200.9mers \
    --run_script cluster_and_extract.sh \
    --condor_job_fn cluster_and_extract_condor_job \
    --submit_script_fn submit_cluster_and_extract_condor_jobs.sh

