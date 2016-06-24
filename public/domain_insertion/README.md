Domain Insertion Demo
=====================

KEYWORDS: DESIGN GENERAL

In this tutorial, we will demonstrate how to perform domain insertion with 
Rosetta.  Domain insertion is when you have two well-folded domains, and you 
insert one (B) into a flexible loop of the other (A), such that domain A is 
split into halves in primary sequence, but the whole protein still folds into A 
and B.  We will also allow for redesign of the remodeled loop accepting the 
insertion.

The sample PDBs are 1EMA (GFP) and 2LCT (an SH2 domain).

A dirty little secret is that there is no proper way to do domain insertion in 
Rosetta.  Nobody's gotten around to writing a real mode for it.  We are instead 
performing what is technically known as an "epic hack" to use a totally 
different suite, AnchoredDesign, for the purpose.  Unfortunately, there are 
many strange nomenclature issues forced by this: the "anchor" is the inserted 
domain, and the "scaffold" is the domain receiving the insertion, and the 
"target" is non-existent for domain insertion but must pretend-exist for the 
purpose of the code.  `AnchoredDesign` and `AnchoredPDBCreator` are 
[[extensively documented|public/anchored_design/README]] with a 
protocol capture released in 3.3, please refer to that documentation for more 
details.

Where the insertion will occur
------------------------------

Open the structures in the Pymol. The goal is to insert SH2 domain into GFP. 
Presumably you will have a scientific problem where you are either interested 
in a particular insertion isomer/position, or you only care that it's a good 
structure but not exactly where the insertion occurs.

We identify the loop region 207-220 (chain A) in GFP arbitrarily as the target 
for insertion.  It would be good to sample many insertion positions and loop 
lengths to find the MOST stable insertion; that is beyond the scope of this 
tutorial but should be easy by extension.

Preparing input structures for the first step, AnchoredPDBCreator
-----------------------------------------------------------------

We will have to do some manual editing to get the input PDBs ready.  Some of 
these tweaks are just Rosetta idiosyncrasies, some are AnchoredDesign issues.

- Prepare the SH2 pdb
    - Delete all the PDB file head matter, up to the first ATOM record.
    - Delete all lines including and after the first ENDMDL card.  This happens 
      to be an NMR model, and we only want the first NMR sub-model in the PDB, 
      and it's best to just delete them now.
    - Delete the peptide chain B out of the pdb.  It is irrelevant to this 
      problem.
    - To make the loop closures more plausible, we chose to delete residues 
      661-668 out of the SH2 domain.  This brings its termini closer together, 
      so that the loop accepting the insertion will not have to stretch to 
      accommodate it.
    - Compare 2lct.pdb to 2lct_prepared.pdb to see these changes.

- Prepare the GFP pdb
    - We will use it as is.  It has a fluorophore (duh), we'll use the flag 
      `-ignore_unrecognized_res` to silently edit it out of the input.  It is 
      not near the surface and won't affect the modeling.

- Prepare the “target” pdb
    - The target has no meaning for domain insertion, but is a requirement of 
      the AnchoredDesign suite.  Here, take a glycine from the SH2 domain and 
      move it far from the SH2 domain.  The easy way to do this is to copy a 
      single GLY residue into a new file and translate it by manually adding 
      many (900) angstroms to its x/y/z coordinates.  You can compare 
      `rosetta_inputs/AnchoredPDBCreator/pseudotarget.pdb` to residue 691 of 
      the SH2 domain to see what we did.

Performing step 1, AnchoredPDBCreator
-------------------------------------

AnchoredPDBCreator will perform the mechanical part of the insertion operation 
(actually suturing the sequences together), but it does not attempt to model 
the new interface much.  We will run it briefly to get a rough inserted 
structure which we will refine later.

To run AnchoredPDBCreator, use its executeable (of the same name):

    export $ROSETTA3=<path/to/Rosetta/source>
    $> cd rosetta_inputs/AnchoredPDBCreator
    $> $ROSETTA3/bin/AnchoredPDBCreator.macosclangrelease @options

Replace `macosclangrelease` with your system settings.

This command will create two new files in that directory.  The output structure 
will be named `S_0001.pdb`; there will also be a scorefile `score.sc`.

In realistic usage, you would generate several hundred models and choose a 
subset to subject to more processing by analyzing the LAM score (reported in 
both the PDB file and score.sc).  Also note that you will increase the value of 
the APDBC_cycles argument to result in longer trajectories.  LAM score 
(LoopAnalyzerMover) attempts to capture how well-closed and formed the subject 
loop is; it is further explained in the AnchoredDesign documentation.  
AnchoredPDBCreator results need only be judged on that criterion; the 
AnchoredDesign protocol will refine it anyway.

To sort through many models, run scripts/sort_by_LAM.sh in the result folder.  
It will sort the scorefile to put the best (lowest) LAM scores at the top.  
Manually examine the best handful and pick your favorite.  It will be used as 
the input to the next step, AnchoredDesign.  Don't stress over your choice here 
– it doesn't have to be a great structure, it's an input not an output.

Creating inputs for AnchoredDesign
----------------------------------

The primary input to AnchoredDesign is the result from AnchoredPDBCreator, 
`S_0001.pdb`.  You will want to load it up in a viewer of your choice for the 
next step.  Note that all numbering from this point on is relative to 
`S_0001.pdb`, NOT the original PDBs.  Also note it is chain B, not chain A; 
chain A is the pseudotarget.

- Creating an anchor file

  The anchor file tells AnchoredDesign what the rigid inserted region is (the 
  insert domain).  Here, it is residues 213-305 in chain B of S_0001.pdb; that 
  is what used to be 2lct_prepared.pdb.  The anchor file is formatted B 213 
  305, see AnchoredDesign documentation for more details.

- Creating a loops file.

  The loops file tells AnchoredDesign what regions are flexible loops.  It will 
  actually treat what used to be one loop plus the insertion as one huge loop, 
  but leave the insertion rigid.  So, our loops file will specify a loop 
  running from the N-terminus of the insert loop to the C-terminus, going 
  through the whole insert domain.  The file is at 
  rosetta_inputs/AnchoredDesign/loopfile; the loop file format documentation is 
  in the manual.

- Creating a resfile

  The resfile is optional; it will allow you to mutate residues in the loop to 
  design a loop that best accepts the insertion.  (If you do not use a resfile, 
  use the flag -packing:repack_only instead to preclude design).  Resfile 
  format documentation is available in the manual.  In our resfile, we have 
  specified that the flexible positions in the loop (the loop, but not the 
  inserted domain) can be designed to any residue.

Running AnchoredDesign to refine the insertion
----------------------------------------------

To run AnchoredDesign, use its executable (of the same name):

    cd <demo directory>/rosetta_inputs/AnchoredDesign

or if you are still in the AnchoredPDBCreator directory,

    $> cd ../AnchoredDesign
    $> $ROSETTA3/bin/AnchoredDesign.macosclangrelease @options
    
Again, replace `macosclangrelease` with your system settings.

AnchoredDesign will remodel the loop containing the insertion and sample the 
pseudo-rigid-body degree of freedom between the SH2 and GFP, while leaving the 
cores of each domain rigid.  It will also (optionally) design the loop region 
to create a loop that best accepts the insertion.

In this tutorial, the settings nstruct, refine_cycles, and perturb_cycles are 
set fairly low for speed.  In production, you will want to turn these flags up 
higher for better results; please see the options file for more details.

Interpreting results
--------------------

(The AnchoredDesign results will contain a meaningless chain A from the 
pseudotarget – delete or ignore it at your leisure.  Also, AnchoredDesign's 
interface metrics refer to the chain A – chain B interface, which won't exist; 
you should ignore those too.)

AnchoredDesign will create PDB files (here of the form S_0001_*.pdb) and a 
scorefile, score.sc.  Interpreting the results requires all your scientific 
intuition.  For a first pass, you can sort the models by total score with the 
command sort_by_score.sh in the scripts directory.  As before, low scores are 
better.  You can examine the other score terms (reported in the score file), 
and even per-residue scores (reported at the end of each PDB), to help you 
decide which model you think is most physically plausible.  You will want to 
examine your models individually in a viewer to pick the best.  Look for 
well-formed interfaces between the two domains and well-closed loops with good 
geometry.
