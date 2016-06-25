# StepWise Monte Carlo (examples for RNA)

KEYWORDS: NUCLEIC_ACIDS STRUCTURE_PREDICTION RNA

# Author
Rhiju Das, rhiju@stanford.edu

# Brief Description

Steps to build a model of a complex RNA fold

# Abstract

This code allows build-up of three-dimensional de novo models of RNAs of sizes up to ~300 nts, given secondary structure and experimental constraints. It can be carried out reasonably automatically, but human curation of submodels along the build-up path may improve accuracy. A fully automated pipeline is also in preparation.

# More docs
This documentation (and more) are available in the on-line docs at:

https://www.rosettacommons.org/docs/latest/rna-denovo-setup.html

# Example Rosetta Command Lines

## Before continuing, make sure you have the correct environment variables set up by 

```bash
   export ROSETTA_TOOLS=<path/to/Rosetta/tools>
   export ROSETTA_BINEXT=[executable extension, example: .default.linuxgccrelease]
$> source $ROSETTA_TOOLS/rna_tools/INSTALL
```

## Make helices

Example cd into step1_helix/

```bash
$> cd step1_helix/
```

and run

```
$> rna_helix.py  -o H2.pdb -seq cc gg -resnum 14-15 39-40 -rosetta_folder=$ROSETTA_TOOLS/../ -extension=$ROSETTA_BINEXT
$> replace_chain_inplace.py  H2.pdb 
```

## Use threading to build sub-pieces

Change to the step2_thread/rosetta_inputs/ directory,

```bash
$> cd ../step2_thread/rosetta_inputs/
```

In the problem above, there is a piece which is a well-recognized motif, 
the UUCG apical loop. Let's model it by threading from an exemplar
of the motif from the crystallographic database. In this directory you will find [1f7y.pdb](http://pdb.org/pdb/explore/explore.do?structureId=1f7y
), which has been downloaded from the RCSB PDB website.


Slice out the motif of interest:
```
$> pdbslice.py  1f7y.pdb  -subset B:31-38 uucg_
```

Thread it into our actual sequence:
```
$> $ROSETTA3/bin/rna_thread.$ROSETTA_BINEXT -s uucg_1f7y.pdb  -seq ccuucggg -o uucg_1f7y_thread.pdb
```

Let's get the numbering to match our actual test case:
```
$> renumber_pdb_in_place.py uucg_1f7y_thread.pdb 24-31
```

Done!

## Build models of a sub-piece denovo

In step3_farfar/, we will see how to setup the Rosetta job for motifs between H2 and H4, using our starting H2 and H4 helices as fixed boundary conditions. 

Change into the `step3_farfar/rosetta_inputs/` directory
```bash
$> cd ../../step3_farfar/rosetta_inputs/
```

There is currently a wrapper script that sets up the job for the rna_denovo executable, which actually runs fragment assembly of RNA with full atom refinement (FARFAR) is not yet equipped to map numbers from our full modeling problem into the subproblem. We have to create it a little sub-problem and map all the residue numberings into the local problem.

There's a file called README_SETUP which has the wrapper command to set up the job. For completeness, the command there is:

```
rna_denovo_setup.py -fasta RNAPZ11.fasta \
    -secstruct_file RNAPZ11_OPEN.secstruct \
    -working_res 14-25 30-40 \
    -s H2.pdb H4.pdb \
    -fixed_stems \
    -tag H2H3H4_run1b_openH3_SOLUTION1 \
    -native example1.pdb 
    -rosetta_folder $ROSETTA_TOOLS/../
    -extension $ROSETTA_BINEXT
```

You don't need to supply a native if you don't have it -- just useful
to compute RMSDs as a reference.

You can run the command by typing:

```
  source README_SETUP
```
or run

```bash
$> rna_denovo_setup.py -fasta RNAPZ11.fasta -secstruct_file RNAPZ11_OPEN.secstruct -working_res 14-25 30-40 -s H2.pdb H4.pdb -fixed_stems -tag H2H3H4_run1b_openH3_SOLUTION1 -native example1.pdb -rosetta_folder $ROSETTA_TOOLS/../ -extension $ROSETTA_BINEXT
```

Then try this:

```
 source README_FARFAR
```
To run a short version of this script, for testing purposes, run:
```bash
$> source README_FARFAR.short
```

Example output after a couple of structures is in example_output/.

[You should probably do a full cluster run -- some tools are available for
 condor, qsub, slurm queueing systems, documented here:

https://www.rosettacommons.org/docs/latest/RNA-tools.html

]

Extract 10 lowest energy models:

```
$> cd ../example_output
$> extract_lowscore_decoys.py H2H3H4_run1b_openH3_SOLUTION1.out 10 -rosetta_folder $ROSETTA_TOOLS/../
```

Inspect in pymol.
(For an automated workflow, you can also cluster these runs and just carry forward the top 5 clusters.)

## Graft together models for the full-length RNA

Change into the `step4_graft/` directory:
```bash
$> cd ../../step4_graft/rosetta_inputs/
```

These were threading and FARFAR solutions that we liked for each submotif -- now we can graft:

```
$> <path/to/Rosetta/main/source>/bin/rna_graft.default.linuxgccrelease -s H2H3H4_run1b_openH3_SOLUTION1.pdb  uucg_1f7y_thread.pdb  H1H2_run2_SOLUTION1.pdb -o full_graft.pdb
```

Done! 

