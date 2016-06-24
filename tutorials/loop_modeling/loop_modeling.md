Loop Modeling
=============

KEYWORDS: LOOPS HOMOLOGY_MODELING GENERAL   
Tutorial by Shourya S. Roy Burman (ssrb@jhu.edu)    
Created 23 June 2016

[[_TOC_]]

Summary
-------
Various loop modeling protocols can be used in Rosetta for various purposes. By the end of this tutorial, you should be able to understand:

* The variety of loop modeling methods in Rosetta
* How to join a break in the protein chain (without missing residues)
* How to model missing segments in proteins
* How to refine segments in proteins
* How to extend the termini
* How to remove peptide segments from a protein (and still get a closed conformation)
* How to combine loop modeling with other protocols

>**This tutorial will not cover the algorithmic details of the loop modeling methods. You will be directed to the documentation explaining the algorithms.**

Loop Modeling Methods
---------------------
>**Loop Modeling in not restricted to segments with a blank DSSP secondary structure assignment. It is more generally applicable to any fragment joining larger segments.**

There are several loop modeling methods present in Rosetta and more are being actively developed. The goal of all loop modeling methods is to sample conformational space of the peptide segment in such a manner that the endpoints of the peptide termini are connected. Here, we will present examples from the following protocols:

* _CCD (Cyclic coordinate descent)_    
This generates loops by fragment insertion from a pre-generated fragment library, and favorably scores conformations which close the loop.
* _KIC (Kinematic closure)_  
This generates loops by analytically calculating possible conformations subject to constraints of the anchored endpoints.
* _Remodel_  
This is not an algorithm in itself but an alternative, user-friendly executable to utilize CCD and KIC
* _Generalized KIC_
This uses the same algorithm as KIC, but can be used for arbitrary backbones, loops through sidechains and other biomolecules. (Links to detailed examples will be provided.)
* _Loop Hash_   
This rapidly searches for peptide conformations using a pre-generated hash map. (Since this is currently under development, we will not be discussing this.)

Navigating to the Demos
-----------------------
The demos are available at `$ROSETTA3/demos/tutorials/loop_modeling`. All demo commands listed in this tutorial should be executed when in this directory. All the demos here use the `linuxgccrelease` binary. You may be required to change it to whatever is appropriate given your operating system and compiler.

<a name="chainbreak_close"></a>Closing Breaks in Protein Chains
--------------------------------
Sometimes you have chain breaks in your protein, perhaps when threading from a homolog. In this case, there are no missing residues, just that the backbone itself is not closed. An example input PDB is provided at `$ROSETTA3/demos/tutorials/loop_modeling/input_files/3gbn_Ab.pdb` where the connection between residue numbers 127 and 128 is severed. To fix this, we will use the kinematic closure protocol (KIC) explained [here](https://www.rosettacommons.org/docs/latest/application_documentation/structure_prediction/loop_modeling/loopmodel-kinematic#purpose). First, we need to write a short _loop file_, detailing which residues are to be modeled as loops and where the cutpoint is.

    LOOP 125 130 0 0 1
    
This 2<sup>nd</sup> and the 3<sup>rd</sup> columns in this file correspond the the loop start and end residues. The 4<sup>th</sup> column indicates the cut point in the [[foldtree|fold_tree]] to allow motions in the loop without propagating them through the rest of the protein. It must be between the start and end residues (both included). 0 is the default option which allows Rosetta to pick a cut point. Note that **the residue numbering in the loop file is not based on PDB numbering but on Rosetta internal numbering.**. The 4th column represents the skip rate, which we have set to 0. Setting the last column to 1 makes Rosetta start building from an extended structure.

To close the loop, run:

    $> $ROSETTA3/bin/loopmodel.linuxgccrelease @flag_basic_KIC 

    
This should take about 2 minutes at the end of which you will produce one PDB and one score file in `output_files`. This PDB will now have a closed loop. In production runs, you should increase `-nstruct 1` to `-nstruct 500` and `-loops:max_kic_build_attempts 200` to `-loops:max_kic_build_attempts 250` in the option file.

A list of further options and documentation can be found [here](https://www.rosettacommons.org/docs/latest/application_documentation/structure_prediction/loop_modeling/loopmodel-kinematic).

<a name="missing_loop"></a>Modeling Missing Loops
----------------------
Modeling missing loops is a difficult problem. We will use cyclic coordinate descent (CCD), which is explained [here](https://www.rosettacommons.org/docs/latest/application_documentation/structure_prediction/loop_modeling/loopmodel-ccd#algorithm). While we can use the `loopmodel` executable to do this, it requires us to place the missing atoms manually in the input PDB file and provide a fragment library. We will instead use the [`remodel`](https://www.rosettacommons.org/docs/latest/application_documentation/design/rosettaremodel) executable to model the missing loop.

In the folder `input_files` you will find the file `3gbn_missing_loops.pdb` which had residues 13-15 of chain H missing. To model this loop, we first need to generate and modify a _blueprint_ file. To generate the file use:

    $> $ROSETTA_TOOLS/remodel/getBluePrintFromCoords.pl -pdbfile input_files/3gbn_missing_loops.pdb -chain H > input_files/3gbn_missing_loops.remodel
    
You should see a file `3gbn_missing_loops.remodel` in `input_files` which looks like:
    
```
...
11 V .
12 K .
13 S .
14 S .
...
```
The file will contain a file of all residues in the H chain with the corresponding Rosetta internal residue number, which is continuous unlike the PDB numbering. The `.` at the end of each line means do nothing to the backbone of this residue.

Now, we add the missing three residues (KPG) in the manner shown below, assigning them residue number 0, identity X and the preferred secondary structure (H: Helix, L: Loops, E: Extended). Next, we ask it to pick the amino acid (PIKAA) K to fill in the first spot, and so on. 

>**You must also specify the preferred secondary structure of the flanking residues and the identity of the residue to replace them with (i.e. themselves). This is crucial to provide backbone flexibility in these residues to close the loop.**

**Make sure that there is no empty line at the end of the bluprint file as it often causes `remodel` to crash with an uninformative error.**

```
...
11 V .
12 K L PIKAA K
0 X L PIKAA K
0 X L PIKAA P
0 X L PIKAA G
13 S L PIKAA S
14 S .
...
```

>**You can build multiple loops in the same simulation by altering the _blueprint file_ to indicate where these loops are.**

Now with this modified blueprint file, run:

    $> $ROSETTA3/bin/remodel.linuxgccrelease @flag_missing_loops


You should see something similar appearing in the log file:
    
```
core.fragment.picking_old.vall.vall_io: Reading Vall library from <path_to_Rosetta_directory>/main/database//sampling/filtered.vall.dat.2006-05-05 ... 
```

This indicates that it is reading in a database of pre-generated fragments from the database to bound the loops. The simulation should take ~1 minute to run and produce a score file and a PDB with a loop of the missing residues in the directory `output_files`. (It will also produce a file called `1.pdb` in the current working directory with the same structure as the output structure, but with different meta information.) Your output PDB will have chains renumbered _A, B_ etc.

This loop will likely not match the loop of the native `3gbn_Ab.pdb` in just one simulation. You need to run this multiple times by changing the `nstruc` option to `500` or more.


`remodel` can a multitude of applications, including design, which you can read about [here](https://www.rosettacommons.org/docs/latest/application_documentation/design/rosettaremodel#algorithm_basic-remodelling-tasks_extension).

>No loop modeling method gives consistent results on loops longer than 12 residues.

Refining Peptide Segments
-------------------------
Say you are not happy with the conformation of a peptide segment in your protein and you want to find a low-energy conformation for that protein. Such a case may arise if a loop has been added incorrectly. To refine the loop (residues 13-15) we inserted in the [Modeling Missing Loops](#missing_loop) section, we will use KIC with flags similar to those in [Closing Breaks in Protein Chains](#chainbreak_close). We will use the following loops file:


    LOOP 11 17 0 0 1

Now run:

    $> $ROSETTA3/main/source/bin/loopmodel.linuxgccrelease @flag_refine_loop
    
This takes ~2 minutes and produces an output PDB and a score file. The output PDB should be closer to the native `3gbn_Ab.pdb` as demonstrated in the PDB `output_files/expected_results/3gbn_refine_loop_0001.pdb`. You need to run this multiple times by changing the `nstruc` option to `500` or more.

Extending the Termini
---------------------
Once again we will use CCD using `remodel` to model the missing C-terminus strand in the H chain of 3GBN. The residues 115-120 of chain H are missing in `input_files/3gbn_missing_cterm.pdb`. In a fashion similar to the example above, we will generate a _blueprint_ file using:

    $> $ROSETTA_TOOLS/remodel/getBluePrintFromCoords.pl -pdbfile input_files/3gbn_missing_cterm.pdb -chain H > input_files/3gbn_missing_cterm.remodel

and the modify the bottom of the file to get:

```
...
113 K .
114 G L PIKAA G
0 X E PIKAA T
0 X E PIKAA T
0 X E PIKAA V
0 X E PIKAA T
0 X E PIKAA V
0 X L PIKAA S
```
Since we know from prior knowledge that the expected secondary structure is a strand for residues 115-119, we will specify them with E. We also need to make the flanking residue (114) backbone mobile as shown above.

Now with this modified blueprint file, run:

    $> $ROSETTA3/bin/remodel.linuxgccrelease @flag_missing_cterm

The simulation should take ~1 minute to run and produce a score file and a PDB with a C-terminus in the directory `output_files`. _This file may be missing the L-chain (a bug in the code for multiple chains), so you may have to manually go an enter it._

This segment will likely not match the C-terminus of the native `3gbn_Ab.pdb` in just one simulation. You need to run this multiple times by changing the `nstruc` option to `500` or more.

Removing a Loop from the Protein
--------------------------------
Removing an existing loop from a protein is a rather tricky thing. Depending on how far the end points of the segment to be deleted were in 3-D space, you may have to make a large portion of the flanking section just to make the gap close. This may distort the fold in some cases. In this example, we will delete residues 101-108 of chain H of the native `3gbn_Ab.pdb` and close the gap. These residues were specifically chosen as the termini of this segment are close in space, thus increasing the chance of gap closure. To use CCD using `remodel`, we will generate the blueprint file using:

    $> $ROSETTA_TOOLS/remodel/getBluePrintFromCoords.pl -pdbfile input_files/3gbn_Ab.pdb -chain H > input_files/3gbn_Ab_deletion.remodel
    
and simply delete the lines for residues 101-108 to get:
```
...
99 H L PIKAA H
100 M L PIKAA M
109 D L PIKAA D
110 V L PIKAA V
...

```
The residues flanking the deletions on both sides (residue numbers 99,100,109,110) must be made mobile so that the backbone can rearrange and close the gap. For your case, you may need to change the number of flanking residues whose backbone is mobile. Now run

    $> $ROSETTA3/bin/remodel.linuxgccrelease @flag_deletion

The simulation should take ~1 minute to run and produce a score file and a PDB with the loop deleted in the directory `output_files`. 

You need to run this multiple times by changing the `nstruc` flag to `500` or more to get the best gap closure.


Combining Loop Modeling with other Protocols
--------------------------------------------
Using the KIC algorithm through the `loopmodel` and `remodel` applications limits the use of the protocol to linear protein backbones. To make covalent bonds in artitrary backbones, sidechains, non-canonicals etc., we use the Generalized KIC protocol. A complete documentation of this protocol is available [here](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/composite_protocols/generalized_kic/GeneralizedKIC).

As this protocol is typically used for advanced use cases, it often requires a combination with other Rosetta protocols. We will not cover examples of using Generalized KIC in this tutorial. A detailed list of usage cases and examples are provided in the documentation linked above. 

>**Generalized KIC is not available as an executable; it needs to be used through RosettaScripts.**

A tutorial on RosettaScripts can here found [[here|rosetta_scripting]] and detailed documentation can be found [here](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts) which lists all the protocols and filters that can be used.

Tips
----
* There is an application called _Floppy Tail_ in Rosetta. It is **not** meant to be used for general purpose loop or termini modeling. It has only been tested on truly disordered termini.
* Make sure that the blueprint file does not have a blank line at the end.
* Make sure that the flanking residues in the blueprint file are set to be have a mobile backbone.
* Neither CCD nor KIC performs well on peptide segments longer than 12 resdiues.
