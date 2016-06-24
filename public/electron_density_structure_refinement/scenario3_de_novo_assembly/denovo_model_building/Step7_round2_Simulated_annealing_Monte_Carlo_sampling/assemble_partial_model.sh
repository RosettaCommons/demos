rt=$PWD
fragdir="../../Step5_round2_Place_fragments_into_density/combined_round1_round2_frags/"
scripts_dir="../../scripts"

echo the candidate_fragment_placements directory is set at $fragdir
echo the scripts dir is set at $scripts_dir

echo gathering results and pick 5 percent low-scoring SAMC trajectories 
$scripts_dir/get_toppercent_samc_results.py

echo clustering results to get rid of redundancy assignment from these low-scoring trajectories
$scripts_dir/cluster_mc_results_based_on_totalscore.py -f 5percent_selected_models.txt  -o  | bash

cd 5percent_selected_models
    echo finding consensus fragment assignment among these low-scoring SAMC trajectories
    ../$scripts_dir/find_consensus_placements_among_SAMC_trajectories.py -f run*txt -d $fragdir | tee consensus_placements.txt

    echo copying fragments from the consensus_placements.txt 
    for frag in `cat consensus_placements.txt   | grep -v "#" | awk '{print $NF}'`; do
        cp $fragdir/$frag .;
        echo $frag was copied
    done
cd $rt

mkdir -p average_model/frags_to_match_fasta_numbering/
cd average_model/frags_to_match_fasta_numbering/
        echo map fragments to fasta numbering
        ../../../../scripts/map_frags_to_seq_numbering.py \
            -f ../../../../input_files/trpv1.fasta \
            -i ../../5percent_selected_models/after_rotation_frags.9.*pdb
        cd ../
    echo making average model from consensus fragments from SAMC
    ../../../scripts/make_average_model_no_superpose.pl frags_to_match_fasta_numbering/after_rotation_frags.9.*.*.*.????.renum.pdb
    if [ ! -f  average.pdb ]; then
        echo "ERROR: average.pdb is not being made. Please contact: wangyr@u.washington.edu to report bug."
        exit
    else
        echo "Congrats!!!!  partial model is made. Now go ahead to complete it using RosettaEM (Step4)."
    fi
cd $rt

        
    


