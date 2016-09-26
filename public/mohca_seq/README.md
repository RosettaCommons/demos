MOHCA-Seq: Modeling RNA with fragmentation data
===============================================

KEYWORDS: NUCLEIC_ACIDS STRUCTURE_PREDICTION

Author: Clarence Cheng (cyucheng at stanford dot edu)

---

The topic of this demo is Fragment Assembly of RNA with Full-Atom Refinement 
(FARFAR), guided by pseudoenergy constraints from Multiplexed •OH Cleavage 
Analysis by paired-end Sequencing (MOHCA-seq) and secondary structure 
information from one-dimensional chemical mapping and mutate-and-map (M2).

This demo was taken from the appendix of the Methods in Enzymology 
chapter: "Modeling complex RNA tertiary folds with Rosetta".  The abstract for 
that chapter is included below:

Reliable modeling of RNA tertiary structures is key to both understanding these 
structures’ roles in complex biological machines and to eventually facilitating 
their design for molecular computing and robotics. In recent years, a concerted 
effort to improve computational prediction of RNA structure through the 
RNA-Puzzles blind prediction trials has accelerated advances in the field. 
Among other approaches, the versatile and expanding Rosetta molecular modeling 
software now permits modeling of RNAs in the 100 to 300 nucleotide size range 
at consistent sub-helical (~1 nanometer) resolution. Our lab’s current 
state-of-the-art methods for RNAs in this size range involve Fragment Assembly 
of RNA with Full-Atom Refinement (FARFAR), which optimizes RNA conformations in 
the context of a physically realistic energy function, as well as hybrid 
techniques that leverage experimental data to inform computational modeling. In 
this chapter, we give a practical guide to our current workflow for modeling 
RNA three-dimensional structures using FARFAR, including strategies for using 
data from multidimensional chemical mapping experiments to focus sampling and 
select accurate conformations.

Running the demo
----------------

Documentation for setting up Rosetta and RNA tools:
* https://www.rosettacommons.org/docs/latest/build_documentation/Build-Documentation
* https://www.rosettacommons.org/docs/latest/application_documentation/rna/RNA-tools

Before attempting the demo, make sure to set your environment variables:
```bash
   export ROSETTA_TOOLS=<path/to/Rosetta/tools>
$> source $ROSETTA_TOOLS/rna_tools/INSTALL
```

Step-by-step instructions for running the demo:

1.  Example FASTA file:

        >3P49_RNA.pdb
        ggauaugaggagagauuucauuuuaaugaaacaccgaagaaguaaaucuuucagguaaaaaggacucauauuggacgaaccucuggagagcuuaucuaagagauaacaccgaaggagcaaagcuaauuuuagccuaaacucucagguaaaaggacggag

    The RNA sequence must be lowercase.
    
    Copy this file, located in 1_helix_preassembly/1_generate_commands/rosetta_inputs/ to this directory:
    ```bash
    $> cp 1_helix_preassembly/1_generate_commands/rosetta_inputs/fasta .
    ```

2.  Example secondary structure file:

        .((((((((......((((((....)))))).(((...((((.....))))..)))........))))))))........(((((......((((((...)))))).(((...((((....((((....)))).....))))..))).......)))))
        
        Copy this file, located in 1_helix_preassembly/1_generate_commands/rosetta_inputs/ to this directory:
        ```bash
        $> cp 1_helix_preassembly/1_generate_commands/rosetta_inputs/secstruct .
        ```

3.  Generate command lines for helix pre-assembly:

        helix_preassemble_setup.py -secstruct [secondary structure file] -fasta [FASTA file]
        
        ```bash
        $> $ROSETTA_TOOLS/rna_tools/bin/helix_preassemble_setup.py -secstruct 1_helix_preassembly/1_generate_commands/rosetta_inputs/secstruct -fasta 1_helix_preassembly/1_generate_commands/rosetta_inputs/fasta
        ```

4.  Example command line for helix pre-assembly:

        rna_denovo -nstruct 100 -params_file helix0.params -fasta helix0.fasta  -out:file:silent helix0.out -include_neighbor_base_stacks  -minimize_rna true -rna::corrected_geo -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -score:weights stepwise/rna/rna_helix -cycles 1000 -output_res_num 2-9 65-72

5.  Run command lines for helix pre-assembly (local):

        source CMDLINES
        
        Short version:
        ```bash
        source CMDLINES.short ## had to make by hand.
        ```

6.  Prepare native/reference structure for Rosetta, if available:

        make_rna_rosetta_ready.py 3P49.pdb
        
        ```bash
        $> make_rna_rosetta_ready.py 2_prepare_reference_model/1_slice_pdb/rosetta_inputs/3P40.pdb
        ```
        
    Outputs reformatted native model as `3p49_RNA.pdb`, to be input to 
    README_SETUP. In the glycine riboswitch example presented here, the 3P49 
    crystal structure includes a protein-binding loop that is not part of the 
    construct used for de novo modeling. Command lines [5] through [14] show 
    how to replace the extraneous residues with a tetraloop matching the 
    experimentally probed construct using a short FARFAR modeling run.

7.  Cut out a segment of a model:

        pdbslice.py [3p49_RNA.pdb] -subset [1-21 36-169] [slice_]
        
        ```bash
        $> pdbslice.py 3p49_RNA.pdb -subset 1-21 36-169 slice_
        ```

    The first input is the model from which you want to excise regions of 
    interest. The second input is the range of nucleotides that you want to 
    keep in your model. The third input is the prefix that will be added to the 
    beginning of the input model’s filename. Here, the protein-binding loop is 
    excised by specifying the range of residues given in the command line.

8.  Renumber a PDB:

        renumber_pdb_in_place.py [slice_3p49_RNA.pdb] [1-21 26-159]
        
        ```bash
        $> renumber_pdb_in_place.py slice_3p49_RNA.pdb 1-21 26-159
        ```

    The first input is the PDB file to be renumbered and the second input is 
    the desired final ranges of sequence positions. Gaps may be intentionally 
    left in the final sequence range to allow for remodeling in the middle of 
    the RNA. Here, a UUUA tetraloop will be built in place of the excised 
    protein-binding loop.

9.  Example README_SETUP for de novo remodeling with a sliced input PDB:

        rna_denovo_setup.py \
            -fasta fasta \
            -secstruct_file secstruct \
            -tag native \
            -working_res 1-159 \
            -s slice_3p49_RNA.pdb \
            -cycles 20000 \
            -ignore_zero_occupancy false \

    Options:

        -fasta [fasta]               input FASTA file
        -secstruct_file [secstruct]  input secondary structure file
        -tag                         name for output files
        -working_res                 specify range of residues to model
        -s slice_3p49_RNA.pdb        see below

    The `-s` flag allows users to input a list of PDB files to use in the 
    modeling; the residues that are part of the input PDB files will not be 
    moved relative to each other, though if multiple PDB files are input, the 
    orientations of the residues in the separate files may change. In this 
    example, the full-atom refinement algorithm will be applied in the same run 
    as fragment assembly.
    
    ```bash
    $> source 2_prepare_reference_model/2_run/rosetta_inputs/README_SETUP
    ```

10. Generate command line for FARFAR modeling:

        source README_SETUP
        
        ```bash
        source 

11. Example README_FARFAR:

        rna_denovo -nstruct 500 -params_file native.params -fasta native.fasta  -out:file:silent native.out -include_neighbor_base_stacks  -minimize_rna true -s slice_3p49_RNA.pdb -input_res  1-21 26-159 -cycles 20000 -ignore_zero_occupancy false -output_res_num  1-159

12. Test command line for FARFAR modeling:

        source README_FARFAR

    This command runs a single local job on your computer. Wait until sampling 
    begins successfully (command line output similar to "Picked Fragment 
    Library for sequence u and sec. struct H ... found 2308 potential 
    fragments"), then cancel the run and submit the job to the cluster.

13. Submit jobs to the cluster:

    rosetta_submit.py README_FARFAR out [16] [1]

    The first number states how many processors to use for the run, while the 
    second number states the maximum time each job will be allowed to run 
    (walltime, in hours). Note that certain supercomputers only allow requests 
    specific multiples of processors (e.g. the Stampede cluster requires a 
    multiple of 16).

14. Concatenate all models from the out folder:

        easy_cat.py out

    Also outputs the number of models in the final silent file to the screen.

15. Extract lowest-energy models to .pdb files for viewing in PyMOL:

    extract_lowscore_decoys.py native.out [1]

    Input the number of lowest-scoring models to extract from the silent file. 
    Here, extract the single lowest-scoring model to use as the native model 
    input for comparison to the de novo models.

16. Rename lowest-score model for use as reference model

    mv native.out.1.pdb 3p49_native_RNA.pdb

17. Example pseudo-energy constraint file:

        [ atompairs ]
        O2' 2 C4' 38 FADE   0 30 15 -4.00  4.00
        O2' 2 C4' 38 FADE -99 60 30 -36.00 36.00
        O2' 1 C4' 44 FADE   0 30 15 -4.00  4.00
        O2' 1 C4' 44 FADE -99 60 30 -36.00 36.00
        O2' 5 C4' 60 FADE   0 30 15 -4.00  4.00
        O2' 5 C4' 60 FADE -99 60 30 -36.00 36.00
        O2' 2 C4' 64 FADE   0 30 15 -4.00  4.00
        O2' 2 C4' 64 FADE -99 60 30 -36.00 36.00
        O2' 25 C4' 54 FADE   0 30 15 -4.00  4.00
        O2' 25 C4' 54 FADE -99 60 30 -36.00 36.00
        O2' 45 C4' 64 FADE   0 30 15 -4.00  4.00
        O2' 45 C4' 64 FADE -99 60 30 -36.00 36.00
        O2' 45 C4' 75 FADE   0 30 15 -4.00  4.00
        O2' 45 C4' 75 FADE -99 60 30 -36.00 36.00
        O2' 32 C4' 88 FADE   0 30 15 -4.00  4.00
        O2' 32 C4' 88 FADE -99 60 30 -36.00 36.00
        O2' 42 C4' 84 FADE   0 30 15 -4.00  4.00
        O2' 42 C4' 84 FADE -99 60 30 -36.00 36.00
        O2' 48 C4' 84 FADE   0 30 15 -4.00  4.00
        O2' 48 C4' 84 FADE -99 60 30 -36.00 36.00
        O2' 55 C4' 88 FADE   0 30 15 -4.00  4.00
        O2' 55 C4' 88 FADE -99 60 30 -36.00 36.00
        O2' 55 C4' 108 FADE   0 30 15 -4.00  4.00
        O2' 55 C4' 108 FADE -99 60 30 -36.00 36.00
        O2' 58 C4' 118 FADE   0 30 15 -4.00  4.00
        O2' 58 C4' 118 FADE -99 60 30 -36.00 36.00
        O2' 67 C4' 119 FADE   0 30 15 -4.00  4.00
        O2' 67 C4' 119 FADE -99 60 30 -36.00 36.00
        O2' 67 C4' 121 FADE   0 30 15 -4.00  4.00
        O2' 67 C4' 121 FADE -99 60 30 -36.00 36.00
        O2' 78 C4' 113 FADE   0 30 15 -4.00  4.00
        O2' 78 C4' 113 FADE -99 60 30 -36.00 36.00
        O2' 78 C4' 135 FADE   0 30 15 -4.00  4.00
        O2' 78 C4' 135 FADE -99 60 30 -36.00 36.00
        O2' 42 C4' 157 FADE   0 30 15 -4.00  4.00
        O2' 42 C4' 157 FADE -99 60 30 -36.00 36.00
        O2' 74 C4' 156 FADE   0 30 15 -4.00  4.00
        O2' 74 C4' 156 FADE -99 60 30 -36.00 36.00
        O2' 100 C4' 148 FADE   0 30 15 -4.00  4.00
        O2' 100 C4' 148 FADE -99 60 30 -36.00 36.00
        O2' 100 C4' 145 FADE   0 30 15 -4.00  4.00
        O2' 100 C4' 145 FADE -99 60 30 -36.00 36.00
        O2' 113 C4' 153 FADE   0 30 15 -4.00  4.00
        O2' 113 C4' 153 FADE -99 60 30 -36.00 36.00
        O2' 135 C4' 154 FADE   0 30 15 -4.00  4.00
        O2' 135 C4' 154 FADE -99 60 30 -36.00 36.00
        O2' 5 C4' 119 FADE   0 30 15 -4.00  4.00
        O2' 5 C4' 119 FADE -99 60 30 -36.00 36.00
        O2' 25 C4' 88 FADE   0 30 15 -0.80  0.80
        O2' 25 C4' 88 FADE -99 60 30 -7.20  7.20
        O2' 37 C4' 62 FADE   0 30 15 -0.80  0.80
        O2' 37 C4' 62 FADE -99 60 30 -7.20  7.20
        O2' 79 C4' 103 FADE   0 30 15 -0.80  0.80
        O2' 79 C4' 103 FADE -99 60 30 -7.20  7.20
        O2' 15 C4' 88 FADE   0 30 15 -0.80  0.80
        O2' 15 C4' 88 FADE -99 60 30 -7.20  7.20
        O2' 32 C4' 108 FADE   0 30 15 -0.80  0.80
        O2' 32 C4' 108 FADE -99 60 30 -7.20  7.20
        O2' 9 C4' 138 FADE   0 30 15 -0.80  0.80
        O2' 9 C4' 138 FADE -99 60 30 -7.20  7.20
        O2' 25 C4' 118 FADE   0 30 15 -0.80  0.80
        O2' 25 C4' 118 FADE -99 60 30 -7.20  7.20

18. Example README_SETUP:

        rna_denovo_setup.py \
            -fasta fasta \
            -secstruct_file secstruct \
            -fixed_stems \
            -no_minimize \
            -tag glycine_riboswitch \
            -working_res 1-159 \
            -native 3p49_native_RNA.pdb \
            -cst_file constraints \
            -staged_constraints \
            -cycles 20000 \
            -ignore_zero_occupancy false \
            -silent helix0.out helix1.out helix2.out helix3.out helix4.out helix5.out helix6.out helix7.out \
            -input_silent_res 2-9 65-72 16-21 26-31 33-35 54-56 39-42 48-51 81-85 155-159 92-97 101-106 108-110 145-147 114-117 139-142 \

    Options:

        -fasta [fasta]	input FASTA file
        -secstruct_file [secstruct]	input secondary structure file
        -fixed_stems	specify whether helices should be fixed
        -no_minimize	specify not to perform full-atom refinement; minimization will be performed in the next stage of modeling
        -tag	name for output files
        -working_res	specify range of residues to model
        -native [native.pdb]	input reference or native model; used for benchmarking cases and will return rms calculations for all models (see command line [5])
        -cst_file [constraints]	input file with pseudoenergy constraints
        -staged_constraints	apply constraints
        -ignore_zero_occupancy false	
        -silent [helix0.out helix1.out …]	input silent files with pre-assembled helices
        -input_silent_res [2-9 65-72 16-21 26-31 …]		specify position ranges of helices in silent files

19. Generate command line for FARFAR modeling:

        source README_SETUP

20. Example README_FARFAR:

        rna_denovo -nstruct 500 -params_file glycine_riboswitch.params -fasta glycine_riboswitch.fasta  -out:file:silent glycine_riboswitch.out -include_neighbor_base_stacks  -minimize_rna false -native glycine_riboswitch_3p49_native_RNA.pdb  -in:file:silent helix0.out helix1.out helix2.out helix3.out helix4.out helix5.out helix6.out helix7.out -input_res  2-9 65-72 16-21 26-31 33-35 54-56 39-42 48-51 81-85 155-159 92-97 101-106 108-110 145-147 114-117 139-142 -cst_file glycine_riboswitch_constraints -staged_constraints -cycles 20000 -ignore_zero_occupancy false -output_res_num  1-159

21. Test command line for FARFAR modeling:

        source README_FARFAR

22. Submit jobs to the cluster:

        rosetta_submit.py README_FARFAR out [96] [16]

23. Concatenate all models from the out folder:

        easy_cat.py out

24. Extract lowest-energy models to .pdb files for viewing in PyMOL:

        extract_lowscore_decoys.py glycine_riboswitch.out [15]

25. Example MINIMIZE:

        parallel_min_setup.py -silent glycine_riboswitch.out -tag glycine_riboswitch_min -proc [96] -nstruct [2000] -out_folder min_out -out_script min_cmdline "-native glycine_riboswitch_3p49_native_RNA.pdb -cst_fa_file glycine_riboswitch_constraints -params_file glycine_riboswitch.params -ignore_zero_occupancy false -skip_coord_constraints"

    The first number states how many processors to use for the run, while the 
    second number is 1/6 the total number of previously generated FARNA models. 
    If you are running on a supercomputer that only allows specific multiples 
    of processors, use an appropriate number for the first input.

26. Generate command lines for full-atom refinement:

        source MINIMIZE

27. Example command line from min_cmdline to run as test:

        rna_minimize -native glycine_riboswitch_3p49_native_RNA.pdb -cst_fa_file glycine_riboswitch_constraints -params_file glycine_riboswitch.params -ignore_zero_occupancy false -skip_coord_constraints -in:file:silent min_out/0/0.silent -out:file:silent min_out/0/glycine_riboswitch_min.out

28. Submit jobs to the cluster:

        rosetta_submit.py min_cmdline min_out [1] [16]

    The first number states how many processors to use for each line in 
    min_cmdline. Here, enter 1 for the first input so that the total number of 
    processors used will be equal to the number of processors entered with the 
    `-proc` flag in command line [12], above. The second number states the 
    maximum time each job will be allowed to run (walltime).

29. Concatenate all models from the min_out folder:

        easy_cat.py min_out

30. Sort models by Rosetta energy and select a subset for clustering:

        silent_file_sort_and_select.py [glycine_riboswitch_min.out] -select [1-60] -o [glycine_riboswitch_min_sort.out]

    The range of models under the -select tag includes 0.5% of the total number 
    of FARNA models generated previously. Outputs a new silent file containing 
    selected number of lowest-energy models.

31. Cluster models:

        cluster -in:file:silent glycine_riboswitch_min_sort.out -in:file:fullatom -out:file:silent_struct_type binary -export_only_low false -out:file:silent cluster.out -cluster:radius [radius]

    Select a radius so that 1/6 of the models in the input sorted silent file 
    are in the largest cluster (cluster0) of models.

32. Copy clustered .out file to a new file to isolate cluster0:

        cp cluster.out cluster0.out

33. Extract lowest-energy models to .pdb files for viewing in PyMOL:

        extract_lowscore_decoys.py cluster0.out [15] –no_replace_names

    Input the number of models in cluster0. The -no_replace_names tag preserves 
    the filenames of the cluster members to reflect their order in the cluster, 
    rather than renaming them in order of Rosetta energy score.

34. Cut out a segment of a model:

        pdbslice.py [3p49_native_RNA.pdb] -subset [2-72 81-159] [slice_kinkturn_]

    Here, the 3P49 crystal structure includes an additional G at position 0, 
    which must be excised to allow the leader sequence to be added to the 5´ 
    end, and the internal linker that forms the kink-turn motif with the leader 
    sequence is also excised to allow remodeling.

35. Renumber a PDB:

        renumber_pdb_in_place.py [slice_kinkturn_3P49_native_RNA.pdb] [10-80 89-167]

    Here, the PDB is renumbered to allow the leader sequence to be added at the 
    5´ end.

36. Example revised FASTA file:

        >3P49_RNA_kinkturn.pdb
        ucggaugaagauaugaggagagauuucauuuuaaugaaacaccgaagaaguaaaucuuucagguaaaaaggacucauauuggacgaaccucuggagagcuuaucuaagagauaacaccgaaggagcaaagcuaauuuuagccuaaacucucagguaaaaggacggag

37. Example revised secondary structure file:

        (((......((((((((......((((((....)))))).(((...((((.....))))..)))........))))))))...)))..(((((......((((((...)))))).(((...((((....((((....)))).....))))..))).......)))))

38. Example README_SETUP for de novo remodeling with a sliced input PDB:

        rna_denovo_setup.py -fasta fasta2 -secstruct_file secstruct2 \
         -fixed_stems \
         -tag glycine_rbsw_kinkturn \
         -working_res 1-167 \
         -s slice_kinkturn_3P49_native_RNA.pdb \
         -cycles 20000 \
         -ignore_zero_occupancy false \

39. Thread an RNA sequence into a template structure:

        rna_thread –in:file:fasta [fasta] -in:file:s [template PDB] –o [output PDB]

    The first input is a FASTA file containing two RNA sequences: 1. the 
    sequence of interest, onto which the structure of the template sequence 
    will be threaded, and 2. the template sequence. The template sequence 
    should be truncated to the regions into which the sequence of interest will 
    be threaded; use hyphens (‘-’) to align the template sequence with the 
    target sequence in the FASTA file. The second input, the template structure 
    in PDB format, should be similarly truncated, using pdbslice.py if 
    necessary. If the template PDB is not correctly formatted for Rosetta 
    modeling, use make_rna_rosetta_ready.py to reformat it. The last input is 
    the name of the output PDB.

Further documentation for RNA threading in Rosetta can be found on the 
RosettaCommons website:  
https://www.rosettacommons.org/docs/latest/rna-thread.html


