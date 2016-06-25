Design with flexible loops
==========================

KEYWORDS: DESIGN LOOPS 

Method 1: RosettaRemodel
------------------------

Documentation: <https://www.rosettacommons.org/docs/latest/application_documentation/design/rosettaremodel>


1.  Preparing the PDB

    Take the PDB 3k2m.pdb. Select chain B and C. Get rid of HETATM. Do fast 
    relax. (where `$ROSETTA3`=path-to-Rosetta/main/source)

        $> $ROSETTA3/bin/relax.default.linuxgccrelease -s input_files/3k2m_bc.pdb -ignore_unrecognized_res -use_input_sc -constrain_relax_to_start_coords -relax:fast -out:file:renumber_pdb

    Let's say the PDB is `3k2m_bc_relax.pdb`. Delete the Chain ID column.

2.  Generating BluePrint

        $> cp scripts/getBluePrintFromCoords.pl .
        $> ./getBluePrintFromCoords.pl -pdbfile input_files/3k2m_bc_relax.pdb > test.blueprint

    The blueprint file has the information to direct the protocol on which 
    residue to design and remodel.

3.  Editing BluePrint

    The interface on the monobody is Residue 182-187. We decided to rebuild with 
    1 residue insertion starting at 183. The blueprint has to be modified in 
    the following way

        ...
        178 G .
        179 E .
        180 D .
        181 S .
        182 A L
        183 G L
        0 x L
        184 Y E
        185 M E
        186 F E
        187 M E
        188 Y .
        189 S .
        190 P .
        191 I .
        ...

    In the above example, "0 x L" will mean eXtension and the secondary 
    structure assined for the inserted region in Loop.

    * Column 1 is the residue postion
    * Column 2 is the residue identity
    * Column 3 is the backbone behavior 

4.  Running the remodel application (you should replace input_files/* with your files)

        > $ROSETTA3/bin/remodel.default.linuxgccrelease -s input_files/3k2m_bc_relax.pdb -remodel:blueprint input_files/test.blueprint -extrachi_cutoff 1 -ex1 -ex2 -use_input_sc -num_trajectory 3 -save_top 1 -use_clusters false -find_neighbors


Method 2: Loopmodel and Fixbb
-----------------------------
We can use `loopmodel` with fragment files and `fixbb` with a resfile to design 
the monobody part of interface.

1.  Preparing the starting PDB

    Same as Step 1 in METHOD 1. We do not need to get rid of the Chai ID column 
    in this case. The relaxed PDB `3k2m_bc_relax.pdb` will be the input PDB.

2.  Creating fragment libraries

    Take the fasta file of 3k2m_bc_relax.pdb. Create  fragment libraries of 
    sizes 9 and 3 locally or through [Robetta Server](robetta.bakerlab.org).

3.  Other input files

    Loop file: `3k2m.loop_file`

    The format is as follows:

        #LOOP  start end cutpoint skip-rate extend
        LOOP 85 89 0 0.0 0
        LOOP 179 185 0 0.0 0

    where

        column1  "LOOP":     The loop file identify tag
        column2  "integer":  Loop start residue number
        column3  "integer":  Loop end residue number
        column4  "integer":  Cut point residue number, >=startRes, <=endRes. default - let LoopRebuild choose cutpoint
        column5  "float":    Skip rate. default - never skip
        column6  "boolean":  Extend loop. Default false.

4.  Running the loopmodel application. For this run, you need to have your fragments ready in the input directory).

        > $ROSETTA3/bin/loopmodel.default.linuxgccrelease  @input_files/flags

    where the flags file consist of following options (edit path to database!):

        
        -in:file:fullatom
        -loops:input_pdb input_files/3k2m_bc_relax.pdb
        -loops:loop_file 3k2m.loop_file
        -loops:frag_sizes 9 3 1
        -loops:frag_files aat000_09_05.200_v1_3 aat000_03_05.200_v1_3 none
        -loops:remodel quick_ccd
        -loops:ccd_closure
        -loops:random_loop
        -out:prefix 3k2m_
        -mut core.io.database
        -nstruct 5

    The actual experiment should have 1000-10000 nstruct.

4.  Running the fixbb application with a resfile:

    Create a list of the output PDBs from loop modeling.

        > $ROSETTA3/bin/fixbb.default.linuxgccrelease -l list -resfile resfile -extrachi_cutoff 1 -ex1 -ex2 -nstruct 5

    Resfile format:

        NATAA
        start
        179 B ALLAAxc
        180 B ALLAAxc
        181 B ALLAAxc
        182 B ALLAAxc
        183 B ALLAAxc
        184 B ALLAAxc
        185 B ALLAAxc

5.  Select PDBs based on total_score

6.  Optional run for binding energy.

