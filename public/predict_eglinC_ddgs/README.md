# Predict DDGs for the eglinC protein

KEYWORDS: DOCKING ANALYSIS

The entire workflow for this demo should be described in a file
named README.dox.  It should describe an entire work flow, with
command lines, tested if possible.

```
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
  extract_chains.pl
	-extracts specified chain of pdb and prints to standard out. execute without arguments for usage.
  sequentialPdbResSeq.pl
	-renumbers specified pdb and prints to standard out. execute without arguments for usage.
  (we didn't use the scripts topN_average.scr or compute_top3avg_energies.scr)

README.dox
  -- A prose or list description of how to perform the protocol

FOR_AUTHORS.txt
  -- A description for the demo creators of what their demo
     should achieve.
  -- Most of what starts in this file should end up in the
     README file as well.
```

## Running
The given starting structure was called starting_files/1CSE.pdb.

We started by:
1. removed the chain we didn't want to do ddG calculations on. the command for doing so is:

    ```
	$> ./scripts/extract_chains.pl starting_files/1CSE.pdb I  > 1CSEi.pdb
    ```
this outputs chain I to the file starting_files/1CSEi.pdb

2. renumbered the crystal structure starting from 1.
    ```
   $> ./scripts/sequentialPdbResSeq.pl  -pdbfile 1CSEi.pdb -res1 1 > 1CSEi.ren.pdb
    ```
    this script takes chain I of 1CSE and renumbers starting from 1, then outputs to 1CSEi.ren.pdb
    **WARNING:** if your pdb has a chainbreak (missing part of the poly-peptide chain), then your numbering will be inconsistent. For example, if you have a chain-break between 12 and 23, it will be renumbered as : 12 13 and so on..

3. minimized the starting input file with harmonic constraints on all C-alpha atoms within 9 Angstrom rmsd.
   this must be run from the directory which contains the pdb-file, otherwise you might get an error. 
he minimization protocol only takes in lists of files, so you need to do the following:
    ```
    $> ls 1CSEi.ren.pdb > lst
    $> $ROSETTA3/bin/minimize_with_cst.default.macosgccrelease -in:file:l lst -in:file:fullatom -ddg::out_pdb_prefix minimize_with_cst        
    ```

(where `$ROSETTA3`=path-to-Rosetta/main/source)

4. if you want to double-check that the minimization worked, you can score the structures as follows:
    ```
    $> ls 1CSEi.ren.pdb > test.lst 
    $> ls minimize_with_cst.1CSEi_0001.pdb >> test.lst
    $> $ROSETTA3/bin/score.default.macosgccrelease -in:file:l test.lst -in:file:fullatom -out:file:scorefile score.chk.fsc 
    ```
    and the score for minimize_with_cst.1CSEi.ren_0001.pdb should be lower than 1CSEi.ren.pdb.
    In this case, 1CSEi.ren.pdb has a score of 35.294 and minimize_with_cst.1CSEi.ren_0001.pdb has a score of -64.426. And for a given input structure you should always converge on a score (you should get the same score for each minimized-input structure).

5. prepare the mutation file:
    mutations should be in the form:
    (Wild-type-residue)(residue-position)(mutant-residue)
    with no spaces in between.
    For example, we have the file: mutations.multiples.txt which has the contents:
    ```
    V13A
    V14G
    V18G
    
    A21F
    E23A
    
    F25A
    ```
    output will be as follows:
    ```
    total 6
    3
    V 6 A
    V 7 G
    V 11 G
    2
    A 14 F
    E 16 A
    1
    F 18 A
    ```

    the new-lines mean that this is a complete batch of mutations. ( In this case we would make a triple mutant: V 13 -> A, V 14 -> G, and V 18 -> G. We would also make a double mutant: A 21 -> F, and E 23 -> A. And finally, we make the single mutant: F 25 -> A).
    **REMEMBER:** keep track of the offset between your initial pdb and renumbered pdb.
    To format the mutations.multiples.txt for input into ddgs run the script mutation_format_ddgs.pl as follows:
    ```
    perl mutation_format_ddgs.pl mutations.multiples.txt offset output_path.mut
    ```
    Where offest is the offset used to prepare the pdb for input into ddgs and output_path is the title you want for the .mut file. The output path is an optional variable. If no output path is provided, the script will print to standard output.

    The full explanation for the mutation format is here:
    http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/app_ddg_prediction.html
 
    Basically, the total keyword specifies how many *total* mutations you want to make. For example, if you make a triple mutant, and double mutant, and a single mutant, the value for total will be 6. Then in order to specify how many mutations you make in each 'batch' you specify with a number followed by a newline. 


6. run the ddg prediction application as follows: 
    ```
    $> $ROSETTA3/bin/ddg_monomer.default.linuxclangrelease -in:file:s minimize_with_cst.1CSEi_0001.pdb -ddg::weight_file soft_rep -ddg::iterations 1 -ddg::dump_pdbs true -ddg::mut_file mutations.multiples.txt -ddg::local_opt_only false -ddg::min_cst false -ddg::mean true -ddg::min -ignore_unrecognized_res 
    ```
    This repacks the wild-type and the mutant structures 5 times, and in order to compute the ddG (which is Emutant - Ewt ) it averages the scores of the mutant and wild-type ensembles of structures.  It uses the soft_rep_design scoring function and it is a fixed backbone protocol. Normally I would run this protocol for 20 iterations , but for the demo purposes I'm only testing with 5 iterations.
    The full explanation of all options is here:
    http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/app_ddg_prediction.html
	The reason we choose these set of options is that it gives pretty reliable results in terms of correlation and stability predictions (in terms of destabilizing, neutral, or stabilizing) and it is relatively quick to run.     

    It dumps out the pdbs in the format:
     repacked_wt_round_X.pdb (where X is from 1-5)
     mut_(mutation_label).pdb where the mutation is listed for the triple mutant (for example) as mut_V6AV7GV11G_round_X.pdb (where X is from 1-5)

    Most importantly, it dumps out ddg-predictions in the file: ddg_predictions.out
     If you take a look at the file you will see a header:
    ```
    ddG: description total fa_atr fa_rep fa_sol fa_intra_rep pro_close fa_pair hbond_sr_bb hbond_lr_bb hbond_bb_sc hbond_sc dslf_ss_dst dslf_cs_ang dslf_ss_dih dslf_ca_dih fa_dun p_aa_pp ref 
    ddG: V6AV7GV11G     8.842    16.575    -2.875    -1.825    -1.496     0.000     0.029     0.000     0.000    -0.089     0.446     0.000     0.000     0.000     0.000    -0.557     0.564    -1.930 
    ```

    The description corresponds to the mutation information, the total is the predicted overall ddg. Positive means the mutation was destabilizing, and negative means the mutation was stabilizing. The numbers after the total correspond to the components that contribute to the overall score. Adding them, should sum to the total. 

     
