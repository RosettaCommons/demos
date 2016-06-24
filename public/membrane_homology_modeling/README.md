Membrane Homology Modeling
==========================

KEYWORDS: MEMBRANES STRUCTURE_PREDICTION

This demo shows how to generate a homology model of a membrane protein.  It 
only works for helical membrane proteins.  In the example, a dopamine receptor 
structure is used as the template for a homologous dopamine receptor.  This is 
a relatively easy homology modeling problem.

Preparation of input files
--------------------------

There is a non-trivial amount of preparatory work that is required before 
running a homology modeling task in Rosetta.

1.  Find a structural template for homology modeling.

    The first step is to find a structural alignment to a template protein.  In 
    this demo, the fasta sequence file for the target (`DXDR_.fasta`) was 
    pasted into the input panel of the HHpred web server.  The server gives you 
    a list of structural homologs ranked by the quality of alignment.  For 
    each, you get a secondary structure prediction and a template structure. 
    The next task is to manipulate this information to create a combined fasta 
    file that shows the alignment between the target and the template file.  In 
    this demo this combined fasta file is named `DXDR_D3DR.HHp`.

2.  Process the template structure to accept the target sequenced.

    Next you need to acquire the pdb file for the template, and you need to 
    renumber the residues in the file to match the target protein numbering. 
    Scripts for performing this are provided.  First, the perl script 
    fastAln2zones.pl is invoked on the combined fasta file:

        fastAln2zones.pl DXDR_D3DR.HHp <zones_file>

    This creates a "zones" file that shows the correspondence between residues 
    in the target and residues in the template.

    Using several of these pieces of information, the script 
    `createTemplate.pl` creates the renumbered template structure:

        createTemplate.pl -zonesfile <zones_file> -fastafile DXDR_.fasta -parentpdb 3PBL_A_renum.pdb -outpdb DXDR_D3DR.pdb

    3PBL_A_renum.pdb was previously created by taking 3PBL.pdb and cutting out 
    the T4 lysozyme insertion from a loop in the membrane domain, and 
    renumbering the remaining residues with the awk script `renum.awk`. 
    Similar manipulations may be required for each specific application.

    The end product of step 2 is the file DXDR_D3DR.pdb

3.  Create a fragment library for the target protein.

    You also need to generate a fragments file as for any other loop/homology 
    modeling application.  See the appropriate demo for fragment generation.

4.  Identify the regions of the protein that must be rebuilt by Rosetta.

    Now we need to tell Rosetta the regions that must be rebuilt by Rosetta. 
    These are the regions that are poorly defined.  This is defined in 
    DXDR.loopfile.  The format is:  first line, a comment, other lines are 
    `LOOP start-res end-res cutpoint skip-rate extend` for each flexible region 
    to rebuild.  The cutpoint column is used to tell Rosetta where to cut the 
    loop.  A zero tells Rosetta to use its default.  Skip rate allows you to 
    spend more time on one loop versus the others.  Here we don't ask for this. 
    Extend "X" indicates not to use any information from the starting template 
    as part of the flexible regions.

5.  Create membrane specific input files.

    There are two membrane protein-specific input files that are required.  The 
    first is a weight file with specialized scoring function weights.  It also 
    triggers the initialization of a number of membrane-specific scoring 
    function capabilities.  The second membrane specific file tells Rosetta 
    where the membrane is.  This file contains more information than is 
    necessary.  It is called DXDR.span.  The format is:

    * Line 1 is a comment.
    * Line 2 is the number of transmembrane helices and total number of 
      residues.
    * Line 3 tells about the topology, only important for ab initio folding, 
      but part of the format.
    * Line 4+ tell the bounds of the transmembrane helices in residue pairs. 
      Only the first pair is used for homology modeling.  The number of lines 
      must match the number of helices given above.

To sum up, we now have:

* `DXDR_D3DR.pdb` – the template pdb file generated above
* `DXDR.span` – the membrane helix information
* `Frags/` a directory with fragment data for the target protein
* `DXDR.loopfile` – the file that tells which regions are to be rebuilt

Running the demo
----------------

The following command line options are used as arguments to the minirosetta 
executable (minirosetta.(os)(options)(mode), for example:

    minirosetta.linuxgccrelease
    minirosetta.macosclangrelease

Provide the path to your own rosetta database (not necessary if you would like to use the database located in the same Rosetta/main/ directory as the executable):

    -database /Users/patrickbarth/RosettaCon2011/tutorial/trunk_r43621/minirosetta_database

Membrane proteins use the same protocol as regular proteins:

    -run:protocol looprelax

Information for fragments files:

    -loops:frag_sizes 9 3 1
    -loops:frag_files ./rosetta_inputs/frags/aaDXDR_09_05.200_v1_3 ./rosetta_inputs/frags/aaDXDR_03_05.200_v1_3

Some input information.  For benchmark purposes, replace the following PDB by 
the experimentally-determined native structure:

    -in:file:native ./rosetta_inputs/DXDR_D3DR.pdb
    -in:file:fullatom
    -s DXDR.pdb

The following PDB corresponds to the starting template:

    -loops:input_pdb /Users/patrickbarth/RosettaCon2011/tutorial/dopamine/input/DXDR_D3DR.pdb

Loop file generated from the zone file:

    -loops:loop_file /Users/patrickbarth/RosettaCon2011/tutorial/dopamine/input/DXDR.loopfile

Membrane specific input information:

    -score:weights ./rosetta_inputs/membrane_highres_t1.wts
    -in:file:spanfile ./rosetta_inputs/DXDR.span

Output format:

    -out:file:fullatom
    -out:file:silent_struct_type binary
    -out:file:silent DXDR_lprlxmb.out

Number of structures generated:

    -nstruct 1

Loop remodeling step: for additional details, see the regular looprelax 
protocol options:

    -loops:remodel  quick_ccd
    #-loops:refine  refine_ccd
    -loops:random_order
    -loops:idealize_before_loop_close

If you want FULL STRUCTURE relax then do this, otherwise don't:

    -loops::relax    fastrelax

Fail on bad H-bonds:

    -fail_on_bad_hbond false

The expected output is a silent file with the model coordinates.  In this demo, 
the name of the file is:

    DXDR_lprlxmb.out.

For example, your full commandline could look like this:

	$> $ROSETTA3/bin/minirosetta.default.linuxgccrelease -run:protocol looprelax -loops:frag_sizes 9 3 -loops:frag_files ./rosetta_inputs/frags/aaDXDR_09_05.200_v1_3 ./rosetta_inputs/frags/aaDXDR_03_05.200_v1_3 -in:file:native ./rosetta_inputs/DXDR_D3DR.pdb -in:file:fullatom -s ./rosetta_inputs/DXDR_D3DR.pdb -loops:input_pdb ./rosetta_inputs/DXDR_D3DR.pdb -loops:loop_file ./rosetta_inputs/DXDR.loopfile -score:weights ./rosetta_inputs/membrane_highres_t1.wts -in:file:spanfile ./rosetta_inputs/DXDR.span -out:file:fullatom -out:file:silent_struct_type binary -out:file:silent DXDR_lprlxmb.out -nstruct 1 -loops:remodel quick_ccd -loops:random_order -loops:idealize_before_loop_close -fail_on_bad_hbond false
