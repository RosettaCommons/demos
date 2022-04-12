# ConfChangeMover for soluble proteins

KEYWORDS: STRUCTURE_PREDICTION GENERAL

Author: Davide Sala

Citation: REF
Authors: 

--------------------------------------------------------------------------------------------------
## Modeling of Adenosylcobinamide kinase apo-to-holo conformational change with distance restraints
In this demo the transition of a small helix and its two adjacent loops (from green to red in the picture) will be modeled with the ConfChangeMover. The small helix in the target structure is longer than in the input conformation and is shifted in the sequence composition, making the modeling quite complex. Its Cα RMSD between native conformations is 10.6 Å. The total Cα RMSD of the region to model (helix + loops) is 9.9 Å. To intensify the sampling of loops, in stage2 three different loops conformations will be created for each helix transition. 


## Preparation of the input files
The apo (PDB ID 1CBU) and holo (PDB ID 1C9K) structures can be downloaded from ProteinDataBank. All the residues except range 1-180 of chain A are removed. In the target holo conformation all the missing residues can be modeled by uploading sequence and template in the user-friendly portal http/swissmodel.expasy.org (click on “start modeling” followed by “User template: on the right). 
Fragments can be collected through htt/old.robetta.org. The website requires the protein fasta sequence to generate 3-mer and 9-mer fragments.  
Four Cα distances involving two residues and the termini of the helix have been derived from the target structure to be used as restraints. 
All the following input files have been placed inside the “input_files” folder.  “1cbu.pdb” and “1c9k.pdb” for input apo and target holo structures, respectively. “aat000_03_05.200_v1_3.txt” and “aat000_09_05.200_v1_3.txt” for 3-mer and 9-mer fragments, respectively. The fasta sequence file “1c9k.fasta”. Restraints are in “1c9k.cst”.


## RosettaScripts XML file
RosettaScripts XML file used for modeling has been deposited in the “scripts” folder. To increase loops variability the frequencies of using template segments and of targeting exclusively gaps have been kept low to 0.1. The “MultiplePoseMover” is required to perform actions on multiple stage2 models out of a single stage1 sampling.
Besides ConfChangeMover, a couple of metrics are calculated: 1) the RMSD from the native structure of the helix (_H suffix) and 2) the RMSD of the whole modeled region (_HL suffix). 

## Running the command to perform the modeling
By executing the command below, 9 models are generated in approximatively 20 minutes on a single core. User may need to change the executable depending on the OS used, MAC/linux/Windows. Rosetta executables are located in “path/to/Rosetta/main/source/bin”. 

$> rosetta_scripts.default.linuxgccrelease -parser:protocol scripts/ccm.xml @input_files/flags.options

## Analyzing results
Depending on the number of models generated user can find the corresponding number of PDB files and a score file. The score file contains the RMSD metrics for each model. Pre-generated output files can be found in the “output_files” folder. 
--------------------------------------------------------------------------------------------------
