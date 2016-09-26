Enzyme Design Demo
==================

KEYWORDS: DESIGN ENZYMES

Tutorial for a complete *de novo* enzyme design run, using the TIM reaction as an 
example, as published in 

* Richter F, Leaver-Fay A, Khare SD, Bjelic S, Baker D (2011) De Novo Enzyme 
  Design Using Rosetta3. PLoS ONE 6(5): e19230. 
  doi:10.1371/journal.pone.0019230

Tutorial written at RosettaCon2011 by Florian Richter (floric at uw dot edu), 
with help from Patrick Conway, Amanda Loshbaugh, Neil King, and Gert Kiss. 

The contents of the demo directory should be:

    starting_files/
     -- directory in which the raw input files are given - these
        are provided for you and serve as the base for your
        tutorial

    rosetta_inputs/
     -- directory in which the modified starting files should
        be placed which will be used as inputs to Rosetta.
        You may need to make modifications like stripping
        extra chains from the input PDB; store the modified
        PDB here and leave the unaltered one in starting_files

    scripts/
     -- python scripts, shell scripts, used in the workflow
     -- awk, grep, sed, etc. command lines

    README.dox
     -- A prose or list description of how to perform the protocol

    FOR_AUTHORS.txt
     -- A description for the demo creators of what their demo
        should achieve.
     -- Most of what starts in this file should end up in the
        README file as well.

Relevant documentation
----------------------

1. The above cited PLoS ONE paper,
2. Documentation for the enzyme design app  
   https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d6/dbc/enzyme_design.html
3. Documentation for the match app  
   https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d7/dfc/match.html)
4. Documentation about the enzdes cstfile format used for both matching and 
   enforcing catalytic geometries during design  
   https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d5/dd4/match_cstfile_format.html)
5. Familiarize yourself with how to generate a .params file and rotamer library 
   for your ligand of interest, as described in the ligand docking app 
   documentation.  
   https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d4/d47/ligand_dock.html  
   (in particular the section "Preparing the small-molecule ligand for docking")

The purpose of this tutorial is to reenact all steps described in the
PLoS ONE paper.

Step 1: Defining a theozyme in Rosetta format
---------------------------------------------
1.  Inputs required:

    A .params file for all residues used in the theozyme. The Rosetta3 database 
    has .params files for all amino acids, but they will have to be generated 
    for the ligand/reaction substrate. Refer to ligand docking documentation 
    under the link above.

2.  Outputs generated:

    A .cst (constraints) file. This file is written by the user and defines the 
    geometry of the theozyme in Rosetta format (see enzdes cstfile 
    documentation linked above) for use in subsequent steps of the design 
    process. 

3.  Defining a theozyme:

    A theozyme is not defined using Rosetta, but usually in one of the
    following ways: 1) through quantum mechanical calculations (see
    references of Ken Houk's work in the PLoS ONE paper), 2) from chemical 
    intuition, or 3) by stealing or using as inspiration a naturally occurring 
    enzyme’s active site. In this tutorial, we will design a novel triose 
    phosphate isomerase (TIM) based upon a naturally occurring TIM active site. 
    Our theozyme will consist of three catalytic residues and the DHAP ligand 
    from the S. cerevisiae TIM (PDB code 1ney).

4.  Making and checking a .cst file:

    Now that we have dreamed up our theozyme, it needs to be expressed in a 
    format that Rosetta can read. In general, the most unambiguous or precisely 
    defined interaction should come first in the .cst file. We use the Rosetta 
    executable CstfileToTheozymePDB.linuxgccrelease to generate .pdb format models 
    from our .cst file so that we can visually check that it defines our 
    theozyme correctly. The command:

        $> $ROSETTA3/main/source/bin/CstfileToTheozymePDB.linuxgccrelease -extra_res_fa rosetta_inputs/1n1.params -match:geometric_constraint_file rosetta_inputs/mocktim_first_2interactions_example.cst

    produces a file called 
    `PDB_Model_mocktim_first_2interactions_example.cst.pdb` in the working 
    directory which can then be visualized with PyMOL. Inspecting this output 
    .pdb file will ensure that the theozyme geometry that is given to Rosetta 
    is as the user intends. See the additional comments on running 
    CstfileToTheozymePDB in the .cst file documentation linked above.

Step 2: Matching
----------------
(i.e. finding suitable sites for the active site theozyme in a library of scaffold proteins)

1.  Inputs required

    * The .params file from Step 1.

    * The .cst file from Step 1.

    * A library of scaffold protein structures in .pdb format. The scaffold 
      library should be as big as possible. Refer to the matcher documentation 
      linked above for instruction on how to prepare a scaffold for matching 
      (which is mainly deciding where the binding site is and what positions 
      will be considered for theozyme residue attachment). In this tutorial, we 
      will match our theozyme into one scaffold, PDB code 1tml, which can be 
      found here: rosetta_inputs/scaffolds/1tml_11.pdb.

    * A scaffold position file for each scaffold that defines which residues in 
      the scaffold structure will be considered for theozyme residue 
      attachment. The format of scaffold position files, and instructions on 
      how to generate them, can be found in the matcher documentation. In this 
      tutorial, our scaffold position file can be found here: 
      rosetta_inputs/scaffolds/1tml_11.pos.

    * Optionally, a ligand grid file that defines where in three-dimensional 
      space the ligand should be placed during matching. In this tutorial we 
      will not use a ligand grid file, but in case one wants to make sure that 
      the ligand is confined to a certain region of space, it is recommended 
      that one be used. More information on ligand grid files can be found in 
      the matcher documentation linked above.

2.  Outputs generated

    * A .pdb file for each “match”. Each match file contains the 
      three-dimensional coordinates of both the scaffold protein and the 
      theozyme, including the ligand, and will be used as an input in step

    * The number of matches found in the scaffold depends on the complexity of 
      the theozyme, and be anywhere between 0 and hundreds. 

3.  Performing matching

    The command line

        $> $ROSETTA3/main/source/bin/match.linuxgccrelease @rosetta_inputs/general_match.flags @rosetta_inputs/1tml_sys.flags

    finds a bunch of matches (~11 in this tutorial) and writes them to the 
    working directory. In a real-life enzdes project, one should look at a few 
    of the matches in PyMOL to make sure that they look roughly as envisioned. 
    Alternatively, the matches could be ranked by similarity to the ideal 
    theozyme (watch for Scott Johnson’s / Ken Houk's EDGE publication).

Step 3: Design
--------------

1.  Inputs required
    * The .params file from step 1.
    * The .cst file from step 1.
    * The matches from step 2.

2.  Outputs generated
    * Designed 'enzymes' in .pdb format.
    * A score file that contains scoring information for each designed enzyme, 
      which will be used in step 4 to evaluate and rank the designs. In this 
      tutorial, our score file generated can be found here: scorefile.txt, 
      while an example score file from a larger design run is provided here: 
      rosetta_inputs/mocktim_all_design_scores.out.

3.  Performing design

    Usually, every match from step 2 is designed several times (between 10-100, 
    depending on computational resources) with the Rosetta enzyme_design 
    executable. For a detailed explanation of what this executable does (and 
    potential options), refer to documentation/paper linked above. Briefly, all 
    residues that are within a given (8A) sphere of the ligand (excepting 
    theozyme residues) are identified and considered changeable in design. 
    These residues are mutated to alanine, and the resulting structure 
    (scaffold with theozyme residues and ligand in a poly-ala cavity) is 
    minimized with respect to the theozyme geometry as specified in the .cst 
    file. Then, 2-3 cycles of constrained sequence design and minimization are 
    carried out, and the final sequence is subjected to an unconstrained repack 
    and minimization. In this tutorial, we will run this stage for one of the 
    matches (rosetta_inputs/UM_1_D41H116K189_1tml_11_mocktim_1.pdb). The 
    command line:

        $> $ROSETTA3/main/source/bin/enzyme_design.linuxgccrelease @rosetta_inputs/general_design.flags -s rosetta_inputs/UM_1_D41H116K189_1tml_11_mocktim_1.pdb -out:file:o scorefile.txt

    generates a designed protein in .pdb format and a score file that has one 
    line of values for several score terms and other metrics. In a real life 
    example, there would be a .pdb file for every design as well as a score 
    file that contains one line of values for each design. An example of such a 
    score file with information about all designs is 
    rosetta_inputs/mocktim_all_design_scores.out, which was taken from the runs 
    done for the calculations in the PLoS paper.


Step 4: Evaluating and ranking designs
--------------------------------------

The purpose of this stage is to reduce the number of candidate designs to an 
amount tractable for visual inspection, and get rid of designs that have 
obvious defects. In short, one is looking for designs that have 1) good 
catalytic geometry, 2) good ligand binding score, 3) a preformed/preorganized 
active site and 4) a well behaved protein (expressible, soluble, stably folded, 
etc).

For 1), the constraint energy is taken as a metric, for 2) there is
the straight rosetta ligand binding score, for 3), the designed site
was repacked without the ligand in the design calculation, and for 4)
the designed protein is compared to the scaffold it came from at the
end of the design calculation. For each of these 4 criteria, there are
terms in the scorefile (refer to documentation linked above). The script 
DesignSelect.pl, which is part of the Rosetta3 distribution, can read the 
output scorefile, as well as a file that specifies required values for certain 
columns, and will then output only those designs in the scorefile that satisfy 
all required values:

    $> $ROSETTA3/src/apps/public/enzdes/DesignSelect.pl -d rosetta_inputs/mocktim_all_design_scores.out -c rosetta_inputs/mocktim_design_selectreqs.txt > selected_designs.txt

These commands will output 44 designs from the 3720 produced for the PLoS ONE 
paper. These 44 would then be visually examined for whether any of them look 
promising enough to be expressed and characterized.
