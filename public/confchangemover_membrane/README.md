Author: Davide Sala

KEYWORDS: STRUCTURE_PREDICTION GENERAL

--------------------------------------------------------------------------------------------------
"Modeling of Rhodopsin active-to-inactive conformational change with simulated DEER distances"

Here, Rhodopsin transition from active to inactive state will be modeled. Such transition involves the three helices TM5-TM6-TM7 that have a Cα RMSD of 4.2 Å, the inclusion of loops does not change RMSD. The cytoplasmic ends of TM5 and TM6 differ in the two protein states. The longer helices in the active state are due to crystal contacts. Thus, in that region the option to deactivate stage2 dihedral constraints will be used and fragments insertion will be performed. 

"Preparation of the input files"
The active (PDB ID 2X72) and inactive (PDB ID 1GZM) structures can be downloaded from ProteinDataBank. All the residues except range 3-326 of chain A are removed. In both the structures all the missing residues can be modeled by uploading sequence and template in the user-friendly portal https://swissmodel.expasy.org (click on “start modeling” followed by “User template: on the right). 
Fragments can be collected through http://old.robetta.org. The website requires the protein fasta sequence to generate 3-mer and 9-mer fragments. 
For membrane proteins, Rosetta requires a trans-membrane span topology file in which protein regions embedded in membrane are specified. There are three ways to generate a span file. 1) PDB file with membrane coordinates can be downloaded from the OPM webserver (http://www.opm.phar.umich.edu) and converted through the “mp_span_from_pdb” executable. 2) trans-membrane residues can be predicted with the OCTOPUS webserver at https://octopus.cbr.su.se and the file converted to the span format with the “octopus2span.pl” script located in the “scripts” folder. 3) trans-membrane residues can be predicted with the TOPCONS webserver at https://topcons.cbr.su.se and the file converted to the span format with the “topcons2span.pl” script located in the “scripts” folder.
Ten DEER distance distributions involving TM5-TM6-TM7 helices have been predicted with the DEER Spin-Pair Distributor in the CHARMM-GUI webserver at https://charmm-gui.org. The median value has been used in the restraints file.
All the following input files have been placed inside the “input_files” folder. “2x72.pdb” and “1gzm.pdb” for input active and target inactive structures, respectively. “aat000_03_05.200_v1_3.txt” and “aat000_09_05.200_v1_3.txt” for 3-mer and 9-mer fragments, respectively. The fasta sequence file “rhodopsin.fasta”.  The span file “rhodopsin.span”. The restraints file “1gzm_simulated.cst”. 

"RosettaScripts XML file" 
RosettaScripts XML file used for modeling has been deposited in the “scripts” folder. The “modify_segments” option is used to include residue 306 in the segment 9 that otherwise would end at residue 305. Because residue 306 is restrained, its inclusion in a segment makes sure that such restraints will be fully exploited in stage 1. The “stage2_residues_no_dihedral_csts” is used to remove dihedral constraints of regions selected through “residue_selectors”. “stage1_multi_sse_freq” defines the frequency of performing multiple randomly chosen SSEs movements. 
Besides ConfChangeMover, a couple of metrics are calculated: 1) the RMSD from the native structure of TM5-TM6-TM7 (_H suffix) and 2) the RMSD of the while modeled region (_HL suffix). 

"Running the command to perform the modeling"
First, copy the scoring function stage2_membrane.wts in the current folder. By executing the command below, 5 models are generated in approximatively 45 minutes on a single core. User may need to change the executable depending on the OS used, MAC/linux/Windows. Rosetta executables are located in “path/to/Rosetta/main/source/bin”. 

rosetta_scripts.linuxgccrelease -in:file:s input_files/2x72.pdb -in:file:native input_files/1gzm.pdb -in:file:spanfile input_files/rhodopsin.span -parser:script_vars frag3=input_files/aat000_03_05.200_v1_3.txt frag9=input_files/aat000_09_05.200_v1_3.txt -epr_deer:input_files input_files/1gzm_simulated.cst -parser:protocol scripts/ccm.xml  -out:nstruct 5 -out:path:pdb ./  -mute core.optimization.Minimizer core.scoring.MembranePotential core.pose.util

"Analyzing results in output"
Depending on the number of models generated user can find the corresponding number of PDB files and a score file. The score file contains the RMSD metrics for each model. Pre-generated output files can be found in the “output_files” folder. 
--------------------------------------------------------------------------------------------------

"This is for testing purpose only"
$> rosetta_scripts.default.linuxgccrelease -in:file:s input_files/1cbu.pdb -in:file:native input_files/1c9k.pdb -parser:script_vars frag3=input_files/aat000_03_05.200_v1_3.txt frag9=input_files/aat000_09_05.200_v1_3.txt cstfile=input_files/1c9k.cst -parser:protocol scripts/ccm.xml  -out:nstruct 1 -out:path:pdb ./  -mute core.optimization.Minimizer 
