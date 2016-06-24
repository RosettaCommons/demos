Basic Homology Modeling Demo
============================

KEYWORDS: STRUCTURE_PREDICTION GENERAL

Conceptual steps
----------------

1. Create fragments for target sequence
2. Find template structure
3. Make a sequence alignment between template sequence and target sequence
4. Download and clean the template pdb
5. Make Rosetta flags file
6. Execute Rosetta
7. Analyze results

Detailed steps
--------------

1.  Create fragments for target sequence

    Obtain a FASTA file of your target sequence. You can get this from 
    [[NCBI|http://www.ncbi.nlm.nih.gov/guide/proteins/]] or you can make the 
    file manually using the following format:

        >template
        AAAAKDLLEKDF

    Open a web browser and go to 
    http://robetta.bakerlab.org/fragmentsubmit.jsp. Paste the target sequence 
    into the sequence input box and click submit. Wait until the server 
    finishes and download all the files produced to your working directory

2.  Find template structure

    Finding a template can be done in multiple ways. One suggestion is to 
    submit the target sequence to HHpred or similar. HHpred can be found here: 
    http://toolkit.tuebingen.mpg.de/hhpred.

3.  Make a sequence alignment between the template sequence and the target 
    sequence

    In order to generate the alignment, you need the primary sequences of both 
    the target and the template. For the template, go to www.pdb.org and search 
    for the PDB code of your template (obtained in Step 2). In this case, we 
    are using 2AST.  In fact, we are using part of chain B from 2AST. On the 
    top right, click the down arrow by "Download Files" and click on "FASTA 
    sequence". Alternatively, you can get the FASTA sequence from NCBI (see 
    Step 1) or make it manually.

    Align the sequences using ClustalW  or the alignment tool of your choice 
    (e.g., EMBOSS, FASTA, Nexus, etc.). Create a FASTA file 
    (template_target.fasta) containing both primary sequences. Go to 
    http://www.ebi.ac.uk/Tools/msa/clustalw3/ and upload the newly created 
    fasta file. The default settings are a good place to start to generate a 
    generally good sequence alignment.  Use the slow alignment for better 
    results Enter your email so you can results. Click submit at the bottom of 
    the page. Depending on the size of the sequences being aligned, it will 
    take anywhere from 5 minutes to some hours or days to run

    After you generate the alignment file, put it in a format like this:

        score 123.456
        t000_			1 VIAFRCPRSFMDQPLAEHFSPFRVQHMDLSNS------VIEVSTL
        2astB_66-105_renumbered	1 ILSLRRSLSYVIQGMANIESLNLSGCYNLTDNGLGHAFVQEIGSL

    where the score is a floating point number, column 1 starting on line 2 is 
    the name of the target (or template, on line 3), column 2 is the beginning 
    sequence position for the threading, and the rest of the line is the 
    sequence.

4.  Download and clean the template PDB

    Once you have decided upon a template, search for the PDB ID at www.pdb.org 
    and go to "Download Files" and select to download the PDB.  This file then 
    needs to be "cleaned" to run properly in Rosetta. To avoid errors when 
    Rosetta reads in the PDB file, the protein must be formatted correctly or 
    “cleaned”. A correctly formatted PDB file includes removed non-ATOM 
    records, renumbered residues from 1, renumbered atoms from 1, and corrected 
    chain ID inconsistencies. The script clean_pdb.py located in the scripts 
    directory will be used to format the template PDB file.

    The clean_pdb.py script requires that you have python2.2 (go to 
    www.python.org for instructions on download and installation) or higher 
    installed, as well as biopython (biopython.org)..  

    Execute the script by typing:

        ./scripts/clean_pdb.py template.pdb A

    The clean_pdb.py script will output two files:

        template_A.pdb
        template_A.fasta

5.  Make the Rosetta flags file

    Rosetta supports a threading protocol (minirosetta application in the bin) 
    which needs a set of input files and specific commandline flags. The 
    commandline flags can be stated in a file in which case the commandline 
    reduces to 

        minirosetta @flags

    The format of the flagsfile is key-value pairs and an example of a working 
    flagsfile to perform thrading looks like:

        -run:protocol threading # run the rosetta threading, loopbuilding, and refinement protocol
        -in:file:alignment ./starting_files/template_target_short.aln # path to alignment file
        -cm:aln_format general # format of alignment file
        -frag3 ./starting_files/fragments/aat000_03_05.200_v1_3.gz # path in 9mer fragment file
        -frag9 ./starting_files/fragments/aat000_09_05.200_v1_3.gz # path to 3mer fragment file
        -in:file:fasta ./starting_files/fragments/t000_.fasta # path to target fasta file
        -in:file:fullatom # we have fullatom format for the input template
        -loops:frag_sizes 9 3 1 # which size fragments are you using?
        -loops:frag_files ./starting_files/fragments/aat000_09_05.200_v1_3.gz ./starting_files/fragments/aat000_03_05.200_v1_3.gz none  # paths to 9mer file, 3mer file, and say "none" for the 1mer file
        -in:file:psipred_ss2 ./starting_files/t000_.psipred_ss2 # path to psipred secondary structure prediction file
        -in:file:fullatom
        -out:nstruct 1 # number of structures you want to build.  Should build at least 1000, but 10,000 would be better
        -in:file:template_pdb ./starting_files/2astB_66-105_renumbered.pdb # path to template pdb
        -database ./rosetta-3.3/rosetta_database/ # path to rosetta database
        -loops:extended # Force extended on loops, independent of loop input file
        -loops:build_initial # Precede loop-modeling with an initial round of just removing the missing densities and building simple loops
        -loops:remodel quick_ccd # closing loops by quick_ccd
        -loops:refine refine_ccd # small movements to remodel loop
        -silent_decoytime # Add time since the last structure was written to score line
        -random_grow_loops_by 4 # randomly grow loops by up to this number of residues on either side.
        -select_best_loop_from 1 # Keep building loops until N and choose best (by score)
        -out:file:fullatom # output in fullatom mode
        -out:output
        -out:file:silent threaded_model.out # silent file stores internal coordinates of the PDB
        -out:file:silent_struct_type binary # output the silent file in binary format
        -out:file:scorefile threaded_model.fasc # output a table of Rosetta scores
        -run:constant_seed # Use a constant seed (1111111 unless specified)
        -run:jran 1111111 # this is good for testing since you should always get the same result
        -overwrite # overwrite any already-existing results having the same name

6.  Execute Rosetta

        $> <path/to/Rosetta>/main/source/bin/minirosetta.macosgccrelease @flags

7.  Analyze results

    There are several options for things to do next, see the following demos:

    * [[analyzing_structure_quality|public/analyzing_structure_quality/readme]]
    * [[clustering|public/clustering/readme]]
    * [[homology_modeling_with_end_extension|public/homology_modeling_with_end_extension/readme]]
    * [[relax_a_large_structure|public/relax_a_large_structure/readme]]

