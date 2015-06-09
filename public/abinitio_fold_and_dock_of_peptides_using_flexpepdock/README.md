AbInitio fold-and-dock of peptides using FlexPepDock
====================================================

This demo will illustrate how to run abinitio docking of peptide on receptor.

Purpose
-------
A wide range of regulatory processes in the cell are mediated by flexible peptides that fold upon binding to globular proteins. The FlexPepDock ab-initio protocols is designed to create high-resolution models of complexes between flexible peptides and globular proteins. Both protocols were benchmarked over a large dataset of peptide-protein interactions, including challenging cases such as docking to unbound (free-form) receptor models (see References).

Algorithm
---------
The input to the ab-initio protocol is: (1) A model of the peptide-protein complex in PDB format,starting from arbitrary (e.g., extended) peptide backbone conformation. It is required that the peptide is initially positioned in some proximity to the true binding pocket, but the exact starting orientation may vary. A preiminary step for the ab-initio protocol is the generation of fragment libraries for the peptide sequence, with 3-mer, 5-mer and 9-mer fragments. Another preliminary step is pre-packing. The first step in the main part of the protocol involves a Monte-Carlo simulation for de-novo folding and docking of the peptide over the protein surface in low-resolution (centroid) mode, using a combination of fragment insertions, random backbone perturbations and rigid-body transformation moves. In the second step, the resulting low-resolution model is refined with FlexPepDock Refinement. As in the independent refinement protocol, the output models are then ranked by the used based on their energy score, or also subjected to clustering for improved performance. 

Refinement vs. ab-initio protocol
---------------------------------
The Refinement protocol is intended for cases where an approximate, coarse-grain model of the interaction is available. The protocol iteratively optimizes the peptide backbone and its rigid-body orientation relative to the receptor protein, in addition to on-the-fly side-chain optimization ( look at /demos/Refinement_of_protein_peptide_complex_using_FlexPepDock/ to learn how to run refinement of protein-peptide complexes).
The ab-initio protocol extends the refinement protocol considerably, and is intended for cases where no information is available about the peptide backbone conformation. It simultaneously folds and docks the peptide over the receptor surface, starting from any arbitrary (e.g., extended) backbone conformation. It is assumed that the peptide is initially positioned close to the correct binding site, but the protocol is robust to the exact starting orientation. The resulting low-resolution models are refined using FlexPepDock Refinement protocol.

Running the FlexPepDock ab-initio protocol
------------------------------------------
1. Create your initial complex structure ( input/1AWR.ex.pdb : extended structure, input/1AWR.pdb : native ).
2. Pre-pack your initial complex.        ( use input/prepack_flags, input/1AWR.ex.ppk.pdb : prepacked extended structure, input for low resolution modeling )
3. Prepare 3-mer, 5-mer and 9-mer fragment files for the peptide using the fragment picker, as in any other Rosetta application (fragment libraries are not required for the receptor). ( scripts to create frags provided in scripts/frags/ )
4. Assuming the receptor chain precedes the peptide chain, offset the indexing of the fragment file to account for it. ( use scripts/frags/shift.sh )

More tips
---------
1. Always pre-pack:
Unless you know what you are doing, always pre-pack the input structure (using the pre-packing mode), before running the peptide docking protocol. Our docking protocol focuses on the interface between the peptide and the receptor. However, we rank the structures based on their overall energy. Therefore, it is important to create a uniform energetic background in non-interface regions. The main cause for irrelevant energetic differences between models is usage of sub-optimal side-chain rotamers in these regions. Therefore, pre-packing the side-chains of each monomer before docking is highly recommended, and may significantly improve the eventual model ranking.
2. Model Selection:
In order to get good results, it is recommended to generate a large number of models ( 50,000 or more for ab-initio). The selection of models should be made based on their score. While selection of the single top-scoring model may suffice in some cases, it is recommended to inspect the top-5 or top-10 scoring models. Our tests indicate that for ab-initio peptide docking, using the score file column labeled "reweighted_sc" may be better than using the default score (score-12). Clustering may also help (see example runs above and protocol capture of the ab-initio protocol). Use clustering script provided in scripts directory.
3. The unbound rotamers flag:
In many cases, the unbound receptor (or peptide) may contain side-chain conformations that are more similar to the final bound structure than those in the rotamer library. In order to save this useful information, it is possible to specify a structure whose side-chain conformations will be appended to the rotamer library during prepacking or docking, and may improve the chances of getting a low-scoring near-native result. This option was originally developed for the RosettaDock protocol.
4. Extra rotamer flags:
It is highly recommended to use the Rosetta extra rotamer flags that increase the number of rotamers used for prepacking (we used the -ex1 and -ex2aro flags in our own runs, but feel free to experiment with other flags if you think you know what you are doing. Otherwise, stick to -ex1 and -ex2aro).

When you should / should not use FlexPepDock
--------------------------------------------
1. For blind docking: 
This protocol is not intended for fully blind docking. The ab-initio protocol assumes that the peptide is located at the vicinity of the binding site, but does not assume anything about the initial peptide backbone conformation. In general, the approximate binding site can be estimated from available computational methods for interface prediction, or from experimental procedures such as site-directed mutagenesis or known homologues. It may be useful to use a constraint file to force the peptide to reach the vicinity of a known binding site or to force specific interactions.
2. Secondary structure assignment: 
Formally, the protocol requires initial secondary structure assignment. It may benefit implicitly from accurate secondary structure prediction when building the fragment libraries, but this is not necessary. 
3. Receptor model: 
This protocol allows full receptor side-chain flexibility, and was shown to perform quite well when docking to unbound receptors or to alternative conformations. However, it is assumed that the receptor backbone does not change too much at the interface, as we do not yet model receptor backbone flexibility. We expect receptor backbone flexibility would be added in future extensions to the protocol.
4. Peptide length:
Our benchmarks consist of peptides of length 5-15, and our protocols performed well on this benchmark regardless of peptide length. We experimented sporadically also with larger peptides (up to 30 residues), but we do not have elaborate benchmark results for these.

Expected Outputs
----------------

The output of a FlexPepDock run is a score file (score.sc by default) and k model structures (as specified by the -nstruct flag and the other common Rosetta input and output flags). The score of each model is the second column of the score file. Model selection should be made based on either the score or I_sc or reweighted-score columns (which exhibited superior performance in the ab-initio benchmarks).

Post Processing
---------------
Except for model selection by total score or I_sc or reweighted score, and possibly clustering, no special post-processing steps are needed. 
All the scripts required are provided in the scripts/ directory

How to run
----------
1. Prepacking:<br/>`$ROSETTA_BIN/FlexPepDocking.linuxgccrelease @input/prepack_flags -database <PATH TO ROSETTA DATABASE>`
2. Low-resolution modeling:<br/>`$ROSETTA_BIN/FlexPepDocking.linuxgccrelease @input/abinitio_flags -database <PATH TO ROSETTA DATABASE>`
3. Selection of top scoring models and clustering

Further information
-------------------
please visit https://www.rosettacommons.org/manuals/archive/rosetta3.4_user_guide/d7/d14/_flex_pep_dock.html


