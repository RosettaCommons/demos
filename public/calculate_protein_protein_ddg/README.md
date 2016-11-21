Calculate Protein Protein ΔΔG
=============================

KEYWORDS: ANALYSIS INTERFACES

A common computational problem involves finding the binding energy of a
protein-protein complex. In this tutorial, we will calculate the change in 
binding energy caused by point mutations in the complex. 

We will use the ubiquitin ligase E3a - ubiquitin conjugating enzyme (UBC)
complex as an example. For instance, when we introduce a mutation in the E3a 
enzyme, the binding reactions for the wild-type and the mutant are as the 
following:

    (1) E + P --> EP    Delta_G   (wild-type or WT)

    (2) E' + P --> E'P  Delta_G'  (mutant or WT)

where E represents the E3a ligase and P represents the UBC protein. 
The binding energy change due to the mutation can be obtained by: 

    (3) Delta_Delta_G = Delta_G' - Delta_G

To obtain Delta_G in Rosetta, we calculate the energy of the complex
using the Rosetta scoring function. We then score the protein and ligase
that have been pulled apart. Subtracting these two scores gives 
the binding energy of the complex, Delta_G.

To obtain a change in binding energy, we create point mutations using resfiles
(explained below). Subtraction of the wild-type Delta_G from the mutant
Delta_G' yields the change in binding energy due to the point mutations.

In this example we will use Rosetta Scripts to peform all of the operations in
a stand-alone application. Using a protocol defined in an XML script we can
perform operations on the complex in a stepwise manner. The protocol will:

1. Relax the input structure to relieve possible clashes in the PDB.
2. Repack the structure.
3. Calculate the Delta_G of the wild type complex. 
4. Repack the structure with a resfile to make a point mutation.
5. Calculate the Delta_G' of the mutated complex.
6. Using Equation (3) to obtain the change in binding energy. The energy value
is in Rosetta energy units. 

Resfiles contain information for how the packer should behave, such as
telling the packer to make a point mutation. For full Resfile documentaion, 
see:

<https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/resfiles>

In root directory of this demo, you will find a sample resfile with the 
following content. It ask Rosetta to mutate residue 641 on chain A (E3a ligase) 
into a tryptophan.  The contents of this resfile are:

    NATAA
    USE_INPUT_SC
    EX1 EX23
    start 
    641 A PIKAA W

Running the demo
----------------
   
If you're doing a number of mutations, it will be more computationally efficient to pre-relax the input structures once, rather than re-optimizing the starting structure for each mutation.
(If you're just doing one mutation, there's a FastRelax option in the XML for the following step which you can enable, and skip the separate relax step.)

    $> $ROSETTA3/bin/relax.default.linuxgccrelease -s starting_files/1C4Z.pdb -ignore_unrecognized_res -out:path:pdb output_files -out:file:scorefile output_files/relax.sc -nstruct 1

The command line arguments are explained below. 

* $ROSETTA3/bin/relax.default.linuxgccrelease:
  the relax executable. You will need to adjust this for your installation path, operating system and compiler type. 

* `-s`:
  input PDB file.

* `-ignore_unrecognized_res`:
  this flag ignores lines in the PDB file that Rosetta doesn't recognize.

* `-out:path:pdb`:
  indicates the output directory for pdbs.

* `-out:file:scorefile`:
  the name for the relax score file.

* `-nstruct`:
  the number of models to output.
  You will likely want to create a number of different relaxed models
  to get an accurate representation of the energy landscape.
  How many depends on your particular system. (Possibly ~1000)
  You can then run the ddG protocol on the best scoring model, 
  or for better estimation of the energy landscape, run it on a selection (or all) of the top scoring models.
  During testing, you may set nstruct to 1. 


You can then run the ddG calculation on your selected structure(s).

    $> cp output_files/1C4Z_0001.pdb 1C4Z_relaxed.pdb

    $> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -parser:protocol mutation_script.xml -s 1C4Z_relaxed.pdb -ignore_unrecognized_res -out:path:pdb output_files -out:path:score output_files -nstruct 1

The command line arguments are explained below. 

* $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease:
  the rosetta_scripts executable. You will need to adjust this for your installation path, operating system and compiler type. 

* `-parser:protocol`:
  this flag indicates the XML file that contains the point mutation ddg protocol.

* `-s`:
  input PDB file.

* `-ignore_unrecognized_res`:
  this flag ignores lines in the PDB file that Rosetta doesn't recognize.

* `-out:path:pdb`:
  indicates the output directory for pdbs.

* `-out:path:score`:
  the output directory for the score file.

* `-nstruct`:
  the number of models to output.
  To get an accurate representation of the energy landscape, many models should be created,
  particularly if you include the FastRelax step in the script. How many depends on your particular system.
  Start small (5-10) and increase it until you get convergence. (Extra models giving the same results as the previous runs.) 
  During testing, you may set nstruct to 1. 

The output file will be a repacked and mutated PDB file with the dg_wt and
dg_mut lines at the end of the file. Subtration of these two numbers yields
the change in binding energy. 




