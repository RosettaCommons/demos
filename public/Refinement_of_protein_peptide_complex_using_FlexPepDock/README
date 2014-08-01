
General purpose
---------------
The FlexPepDock Refinement protocol is designed to create high-resolution models of complexes between flexible peptides and globular proteins, with side chains of binding motifs modeled at nearly atomic accuracy. The Refinement protocol is intended for cases where an approximate, coarse-grain model of the interaction is available. The protocol iteratively optimizes the peptide backbone and its rigid-body orientation relative to the receptor protein, in addition to on-the-fly side-chain optimization. 

The protocol is able to account for a considerable diversity of peptide conformations within a given binding site. However, it is not intended for cases where the peptide backbone is unknown (the ab-initio protocol is best optimized for these cases - see demo in Demos/public/abinitio_fold_and_dock_of_peptides_using_FlexPepDock/ ).

Protocol workflow
-----------------

1. prepacking step

This preliminary step involves the packing of the side-chains in each monomer to remove internal clashes that are not related to inter-molecular interactions. The prepacking guarantees a uniform conformational background in non-interface regions, prior to molecular docking.
The resulting files will be: a start.ppk.pdb structure which will be used as a starting structure in the refinement, a ppk.score.sc score file and a prepack.log file of the run. 

2. docking step

This is the main part of the protocol. In this step, the peptide backbone and its rigid-body orientation are optimized relative to the receptor protein using the Monte-Carlo with Minimization approach, in addition to on-the-fly side-chain optimization. An optional low-resolution (centroid) pre-optimization may improve performance further.

3. Model selection

This step is not performed in the demo since there are several possible ways of processing results, so we only provide here recommendations. 
In any case, desicion should be made after looking at the RMS vs total score plot, but there are additional scoring terms that can be considered when analyzing results: I_sc (score of interface residues only), pep_sc (score of the peptide residues and the interface residues), reweighted_sc (calculated sum of I_sc, pep_sc and total_score). 

You can either decide to concentrate on the top models accoring to one of these different terms, or if you feel that you detect in the plots several energetic minima, you can cluster the results and look at the top result from the largest clusters. For this purpose a dedicated clustering script is provided under scripts/clustering/cluster.sh (usage: cluster.sh <pdb name> <number of top structures> <radius> <score file name> <native file name> <list of structures to be clustered> <stype=column number of score according to which to select the top x structures>) 
