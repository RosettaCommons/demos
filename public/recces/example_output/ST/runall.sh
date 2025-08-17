    recces_turner -score:weights stepwise/rna/turner -seq1 gu -seq2 ac -n_cycle 300000 -temps -1 -out_prefix kT_inf -save_score_terms &
    recces_turner -score:weights stepwise/rna/turner -seq1 gu -seq2 ac -n_cycle 9000000 -temps 0.8 1 1.4 1.8 3 7 30 -st_weights 0 7.33 14.6 17.32 18.87 18.34 17.09 -out_prefix ST -save_score_terms &

