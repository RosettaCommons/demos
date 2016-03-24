AbInitio fold-and-dock of peptides using FlexPepDock
----------------------------------------------------
This demo will illustrate how to run FlexPepDock ab-initio folding and docking of peptide on receptor. The FlexPepDock ab-initio protocols is designed to create high-resolution models of complexes between flexible peptides and globular proteins. Given the approaximate the location of the peptide binding site the ab-initoio procol samples both rigid-body and torsional space extensively. No prior knowledge about the peptide backbone is necessary as the protocol uses fragments to sample peptide backbone conformational space rigorously.

Protocol Overview
-----------------
The input to the ab-initio protocol is a model of the peptide-protein complex in PDB format,starting from arbitrary (e.g., extended) peptide backbone conformation. It is required that the peptide is initially positioned in some proximity to the true binding pocket, but the exact starting orientation may vary. A preiminary step for the ab-initio protocol is the generation of fragment libraries for the peptide sequence, with 3-mer, 5-mer and 9-mer fragments. Another preliminary step is pre-packing. The first step in the main part of the protocol involves a Monte-Carlo simulation for de-novo folding and docking of the peptide over the protein surface in low-resolution (centroid) mode, using a combination of fragment insertions, random backbone perturbations and rigid-body transformation moves. In the second step, the resulting low-resolution model is refined with FlexPepDock Refinement. As in the independent refinement protocol, the output models are then ranked by the used based on their energy score, or also subjected to clustering for improved performance.

Refinement vs. ab-initio protocol
---------------------------------
The Refinement protocol is intended for cases where an approximate, coarse-grain model of the interaction is available which is close to the correct solution both in Cartesian and dihedral (phi, psi) space. The protocol iteratively optimizes the peptide backbone and its rigid-body orientation relative to the receptor protein including on-the-fly side-chain optimization ( look at /demos/refinement_of_protein_peptide_complex_using_FlexPepDock/ to learn how to run refinement of protein-peptide complexes).
The ab-initio protocol extends the refinement protocol considerably, and is intended for cases where no information is available about the peptide backbone conformation. It simultaneously folds and docks the peptide over the receptor surface, starting from any arbitrary (e.g., extended) backbone conformation. It is assumed that the peptide is initially positioned close to the correct binding site, but the protocol is robust to the exact starting orientation. The resulting low-resolution models are refined using FlexPepDock Refinement protocol.

Running the FlexPepDock ab-initio protocol
------------------------------------------
1. Create an initial complex structure: An initial model can be built by placing the peptide in close promity to the binding site in an arbitary conformation. In this we have have provided 2A3I.ex.pdb in which the peptide is in extended conformation. This will serve as the starting strucure and our goal will be to model it to obtain a near-native model with a helical peptide conformation. The native structure is 2A3I.pdb. Both 2A3I.ex.pdb and 2A3I.pdb are located in the input directory.

2. Prepacking the input model: This step involves the packing of the side-chains in each monomer to remove internal clashes that are not related to inter-molecular interactions. The prepacking guarantees a uniform conformational background in non-interface regions, prior to refinement. The prepack_flags file contains the flags for running the prepacking job. The run_prepack script will run prepacking of the input structure 2A3I.ex.pdb located in the input directory.

You need to change the paths of the Rosetta executables and database directories in the run_prepack script (also for run_refine; see below).

  ROSETTA_BIN="rosetta/main/source/bin"
  ROSETTA_DB="rosetta/main/database/"

After changing the paths run the run_prepack script as:
   $./run_prepack

The output will be a prepacked structure, 2A3I.ex.ppk.pdb located in the input directory; a scorefile named ppk.score.sc and a log file named prepack.log file located in the output directory. This prepacked structure will be used as the input for the refinement step.

3. Create 3mer, 5mer & 9mer fragment: The scripts necessary for creating fragments are provided in the fragment_picking directory.
    a. Go to fragment_picking directory.
    b. Save the peptide sequence in the xxxxx.fasta file.
    c. Run make_fragments.pl script to generate PSIPred secondary structure and PSI-Blast sequence profile. You need to chnage the paths in the upper section of the make_fragments.pl file.
       Run as $perl make_fragments.pl -verbose -id xxxxx xxxxx.fasta
       This will create xxxxx.psipred_ss2, xxxxx.checkpoint along with other files.
    d. Run the executable fragment_picker.linuxgccrelease to create the frags. The flags are provided in the flags file and fragment scoring weights are provided in the psi_L1.cfg file.
    Run as $ROSETTA_BIN/fragment_picker.linuxgccrelease -database $ROSETTA_DB @flags >log
    e. Change the fragment numbering using shift.sh script.
    Run as $bash shift.sh frags.500.3mer X >frags.3mers.offset ; where X is the number of residues in the receptor. Do the same for 5mer and 9mer frags
    The offset fragment file will be used as input to the FlexPepDock abinitio protocol. Put then in input/frags directory.

4. Ab-initio folding and docking the prepacked model: This is the main part of the protocol. In this step, the peptide backbone and its rigid-body orientation are optimized relative to the receptor protein using the Monte-Carlo with Minimization approach, including periodic on-the-fly side-chain optimization. The peptide backbone conformational space is extensivyly sampled using fragments derived from solved structures. The file abinitio_flags contains flags for running the ab-initio job. The run_abinitio script will run ab-initio modeling of the prepacked structure generated in the prepacking step located in the input directory.

After changing the Rosetta related paths run the run_abinitio script as:
    $./run_abinitio

The output will be a optimized structure (2A3I.ex.ppk_0001.pdb) located in the output directory; a scorefile named abintio.score.sc and a log file named abinitio.log file located in the output directory. This script has to be modified to run on a cluster during a production run.


Post Processing
---------------
Clustering is performed to select representative models. A clustering scripts is provided in the clustering directory. It will cluster top 500 decoys based on reweighted_sc with radius 2 A. Top scoring members from each cluster is reported in the cluster_list_reweighted_sc_sorted and similar file.
Runs as
 $bash cluster.sh 2.0 ../input/2A3I.ex.pdb ../../output/decoys.silent
The decoys.silent file a compress output formation convenient for production run. See https://www.rosettacommons.org/manuals/rosetta3.1_user_guide/app_silentfile.html for details.

Further information
-------------------
A detailed documentation on FlexPepDock is available at https://www.rosettacommons.org/docs/latest/application_documentation/docking/flex-pep-dock
Raveh B, London N, Zimmerman L, Schueler-Furman O (2011) Rosetta FlexPepDock ab-initio: Simultaneous Folding, Docking and Refinement of Peptides onto Their Receptors. PLoS ONE 6(4): e18934. doi: 10.1371/journal.pone.0018934

