#AbInitio fold-and-dock of peptides using FlexPepDock
----------------------------------------------------

KEYWORDS: PEPTIDES DOCKING

This demo illustrates how to run FlexPepDock ab-initio folding and docking of a peptide onto its receptor. The FlexPepDock ab-initio protocol is designed to generate high-resolution models of complexes between flexible peptides and globular proteins, given the approximate location of the peptide binding site. The ab-initio procol samples both rigid-body orientation and torsional space of the peptide extensively. No prior knowledge about the peptide backbone is necessary, as the protocol uses fragments to sample peptide backbone conformational space rigorously.

Protocol overview
-----------------
The input to the ab-initio protocol is a model of the peptide-protein complex in PDB format,starting from arbitrary (e.g., extended) peptide backbone conformation. It is required that the peptide is initially positioned in some proximity to the true binding pocket, but the exact starting orientation may vary.
Preliminary steps: (1) Generation of fragment libraries for the peptide sequence, including 3-mer, 5-mer and 9-mer fragments. (2) Pre-packing of the receptor and peptide to remove internal clashes that might confuse ranking.
Main protocol: Step 1: Monte-Carlo simulation for de-novo folding and docking of the peptide over the protein surface in low-resolution (centroid) mode, using a combination of fragment insertions, random backbone perturbations and rigid-body transformation moves. Step 2: The resulting low-resolution model is refined with FlexPepDock Refinement. As in the independent refinement protocol, the output models are then ranked based on their energy score, after their clustering for improved coverage of distinct conformations.

Refinement vs. ab-initio protocol
---------------------------------
The Refinement protocol is intended for cases where an approximate, coarse-grain model of the interaction is available that is close to the correct solution both in Cartesian and dihedral (phi, psi) space. The protocol iteratively optimizes the peptide backbone and its rigid-body orientation relative to the receptor protein including on-the-fly side-chain optimization (look at /demos/refinement_of_protein_peptide_complex_using_FlexPepDock/ to learn how to run refinement of protein-peptide complexes).
The ab-initio protocol extends the refinement protocol considerably, and is intended for cases where no information is available about the peptide backbone conformation. It simultaneously folds and docks the peptide over the receptor surface, starting from any arbitrary (e.g., extended) backbone conformation. It is assumed that the peptide is initially positioned close to the correct binding site, but the protocol is robust to the exact starting orientation. The resulting low-resolution models are refined using the FlexPepDock Refinement protocol.

Running the FlexPepDock ab-initio protocol
------------------------------------------
1. Generate an initial complex structure: An initial model can be built by placing the peptide in close promity to the binding site in an arbitary conformation. In this demo we have provided a starting structure with a peptide in extended conformation (2A3I.ex.pdb). Our goal is to optimize this structure using ab-initio FlexPepDock, towards a near-native model with a helical peptide conformation. Both the native structure (2A3I.pdb), as well as the starting structure (2A3I.ex.pdb) are provided in the input directory.

2. Prepack the input model: This step involves the packing of the side-chains in each monomer to remove internal clashes that are not related to inter-molecular interactions. The prepacking guarantees a uniform conformational background in non-interface regions prior to refinement.

Run this eample as:

```bash
$> $ROSETTA3/bin/FlexPepDocking.default.linuxgccrelease @prepack_flags
```

The output will be a prepacked structure, 2A3I.ex.ppk.pdb, located in the input directory; a scorefile named ppk.score.sc and a log file named prepack.log file located in the output directory. This prepacked structure will be used as the input for the ab-initio modeling step.

3. Create 3mer, 5mer & 9mer (peptide lingth >=9) fragment libraries: The scripts necessary for creating fragments are provided in the fragment_picking directory.
    a. Go to the fragment_picking directory.
    b. Save the peptide sequence in the xxxxx.fasta file.
    c. Run the make_fragments.pl script to generate the PSIPred secondary structure and PSI-Blast sequence profiles. You need to chnage the paths in the upper section of the make_fragments.pl file.
       Run as $perl make_fragments.pl -verbose -id xxxxx xxxxx.fasta
       This will create xxxxx.psipred_ss2, xxxxx.checkpoint along with other files.
    d. Run the executable fragment_picker.linuxgccrelease to create the frags. The flags are provided in the flags file and fragment scoring weights are provided in the psi_L1.cfg file.
    Run as $ROSETTA_BIN/fragment_picker.linuxgccrelease -database $ROSETTA_DB @flags >log
    e. Change the fragment numbering using shift.sh script.
    Run as $bash shift.sh frags.500.3mer X >frags.3mers.offset ; where X is the number of residues in the receptor. Do the same for 5mer and 9mer frags
    The offset fragment files will be used as input to the FlexPepDock ab-initio protocol. Put them in the input/frags directory.

4. Ab-initio folding and docking of the prepacked model: This is the main part of the protocol. In this step, the peptide backbone and its rigid-body orientation are optimized relative to the receptor protein using the Monte-Carlo with Minimization approach, including periodic on-the-fly side-chain optimization. The peptide backbone conformational space is extensively sampled using fragments derived from solved structures. The file abinitio_flags contains flags for running the ab-initio job. The run_abinitio script will run ab-initio modeling of the prepacked structure generated in the prepacking step located in the input directory.

After changing the Rosetta related paths run the run_abinitio script as:
    $./run_abinitio

The output will be an optimized structure (2A3I.ex.ppk_0001.pdb) located in the output directory; a scorefile named abintio.score.sc and a log file named abinitio.log, located in the output directory. This script has to be modified to run on a cluster during a production run (see below).


Specific changes needed for a production run
--------------------------------------------
For a production run it is recommended to generate large number of decoys (~10,000 to 50,000). In such a case you can run the job on a cluster. It is advided to use silent output format in such scenario to save space (See https://www.rosettacommons.org/manuals/rosetta3.1_user_guide/app_silentfile.html for details). Include the following lines to the abinitio_flags file:

-out:file:silent_struct_type binary
-out:file:silent decoys.silent

This will create the decoys.silent file containing data related to all the decoys in a compressed format. You can extract speicific decoy using the extract_pdbs.linuxgccrelease executable.
For example:
  $ROSETTA_BIN/extract_pdbs.linuxgccrelease -database $ROSETTA_DB -in:file:silent decoys.silent -in:file:tags 2A3I.ex.ppk_1234.pdb


Along with changes in the flags file you need to modify the run_abinitio file to run on a cluster. The file run_abinitio_slurm is the modified version of run_abinitio adapted to run on a slurm cluster. You should ask your cluster manager for relevent changes required.


Post Processing after a production run
--------------------------------------
In order to diversify our prediction, we cluster the results and select representative models. A clustering scripts is provided in the clustering directory. It will cluster the top 500 decoys based on a cutoff radius of 2.0 Angstrom, and select for each the top-scoring member (according to reweighted score,  reweighted_sc). A top scoring member from each cluster is reported in the file  cluster_list_reweighted_sc_sorted.

Runs as
$bash cluster.sh 2.0 ../input/2A3I.ex.pdb ../output/decoys.silent

Further information
-------------------
Detailed documentation on ab initio FlexPepDock is available under: https://www.rosettacommons.org/docs/latest/application_documentation/docking/flex-pep-dock.
Please cite: Raveh B, London N, Zimmerman L, Schueler-Furman O (2011) Rosetta FlexPepDock ab-initio: Simultaneous Folding, Docking and Refinement of Peptides onto Their Receptors. PLoS ONE 6(4): e18934. doi: 10.1371/journal.pone.0018934

