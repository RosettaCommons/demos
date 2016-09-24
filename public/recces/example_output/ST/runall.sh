    recces_turner -score:weights stepwise/rna/turner -seq1 gu -seq2 ac -rna:farna:thermal_sampling:n_cycle 300000 -rna:farna:thermal_sampling:temps -1 -rna:farna:thermal_sampling:out_prefix kT_inf -rna:farna:thermal_sampling:save_score_terms &
    recces_turner -score:weights stepwise/rna/turner -seq1 gu -seq2 ac -rna:farna:thermal_sampling:n_cycle 9000000 -rna:farna:thermal_sampling:temps 0.8 1 1.4 1.8 3 7 30 -st_weights 0 7.33 14.6 17.32 18.87 18.34 17.09 -rna:farna:thermal_sampling:out_prefix ST -rna:farna:thermal_sampling:save_score_terms &

