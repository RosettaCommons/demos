../../../scripts/setup_placements_condor_jobs.py \
    --mapfile     ../../../input_files/transmem.mrc \
    --fasta       ../../../input_files/trpv1.fasta \
    --fragfile    ../../../input_files/t001_.200.9mers  \
    --condor_job_fn placement_condor_job \
    --submit_script_fn submit_placement_condor_jobs.sh
