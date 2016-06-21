# Analysis
This tutorial will go through some examples of what to do with the output that Rosetta produces. In general the output is one or several of the following files:

* pdb file
* silent file
* score file
* log file

###1. PDB file

The models that Rosetta produces are in pdb file format. If you open the file, you will of course find all atom positions. Additionally, Rosetta will append all scoring terms for every residue at the end of the pdb file.

* Score the provided pdb file and open the scored structure in your favourite text editor.

        $> ../../../main/source/bin/score.jd2.linuxclangrelease -s 1ubq.pdb -out:output -out:file:scorefile score.sc
    
* Now, search for the word "pose". This will bring you to the end of the coordinate section. It should look similar to this:

        # All scores below are weighted scores, not raw scores.
        #BEGIN_POSE_ENERGIES_TABLE S_0004_design_0062.pdb
        label fa_atr fa_rep fa_sol fa_intra_rep fa_elec pro_close hbond_sr_bb hbond_lr_bb hbond_bb_sc hbond_sc dslf_fa13 rama omega fa_dun p_aa_pp ref total
        weights 0.8 0.44 0.75 0.004 0.7 1 1.17 1.17 1.17 1.1 1 0.2 0.5 0.56 0.32 1 NA
        pose -513.648 51.38 327.036 1.19669 -61.0293 2.86533 -26.313 -42.5744 -9.85907 -17.8907 -11.8656 -9.94805 10.566 150.226 -17.7945 -29.4897 -197.143
        THR:NtermProteinFull_1 -2.22882 0.1511 2.60394 0.00934 -0.65245 0 0 0 -0.29914 -0.28075 0 0 0.10172 0.02367 0 0.16454 -0.40685 
        VAL_2 -2.91409 0.33743 1.30547 0.01945 -0.14243 0 0 0 0 0 0 -0.10417 0.10036 0.61358 0.17219 0.74484 0.13264
        
    
 
 |  line starts with     | What does this mean? |
 |--------|--------------------------------------------|  
 | label | names of the included energy terms (scores) |
 | weights | every scoring term has a weighting factor |
 | pose    | every scoring term as the sum of all residues |
 | THR:NtermP[...] | scores for residue #1 (N-terminus, threonine)   |
 | VAL_2 | scores for residue #2 (here a Valine)
 
* Extract the total score for each individual residue:  

        $ egrep '[^A-Z]*_[0-9]+\s' 1ubq_0001.pdb | cut -d'_' -f2 | awk '{print $1 "\t" $NF}' > total_score_per_res.dat   
 This selects all lines that contains the pattern "Letters_[some-number] ", then gets rid of everything before the underscore and finally prints out the first field (the number after _), a tab, and the very last field (here: the total_score)

* Now, you can plot the residue position vs. total score, e.g. as a histogram. You should notice, that there are regions with relatively bad energies (around positions 35-40). 

### The score file