Calculate relative binding affinities for a RNA-protein complex
=====================================================================================

KEYWORDS: NUCLEIC_ACIDS RNA INTERFACES

Written in April 2018 by Kalli Kappel (kappel at stanford dot edu).  

This demo shows how to use the Rosetta-Vienna ddG method to calculate relative binding affinites for a RNA-protein complex.  


## Setting up the demo: 

1. Make sure that you have python (v2.7) installed.
2. Install Rosetta RNA tools. See instructions and documentation [here](https://www.rosettacommons.org/docs/latest/application_documentation/rna/RNA-tools).
3. Install the ViennaRNA package. See instructions [here](https://www.tbi.univie.ac.at/RNA/).

## Brief explanation of required input files: 

1. A structure of a RNA-protein complex, here `start_structure.pdb`. This is the MS2 coat protein-RNA hairpint complex. 
2. A list of sequences for which to calculate binding affinities relative to the sequence found in the starting structure, here `mutant_list.txt`. This is a text file specifying the sequences for which we want to calculate relative binding affinities. One sequence should be specified per line. These can either be the full sequence of the complex (RNA and protein), or just the RNA sequence. If the protein sequence is not specified, then no mutations to the protein will be made. Here, `mutant_list.txt` contains:

```
ugaggcucaccca
ugaggagcaccca
```

So, we will be calculating relative binding affinities for two sequences specifying mutations in the RNA.  

## Step 1: Relax the starting structure.

First, we need to "relax" the starting structure into the Rosetta score function that we'll be using for the ddG calculations.  

&nbsp;&nbsp;&nbsp;&nbsp;**1.1**&nbsp;&nbsp;&nbsp;&nbsp; Set up the relaxation run. We need to provide the starting struture (this will be our "wildtype" -- mutant binding affinities will be calculated relative to the sequence in this PDB file). Here, we're only relaxing the structure once (--nstructs 1), but for actual runs it is recommended that the relaxation is performed 100 times (--nstructs 100). Type:  

```
python PATH_TO_ROSETTA/main/source/src/apps/public/rnp_ddg/relax_starting_structure.py --start_struct start_structure.pdb --nstructs 1 --rosetta_prefix PATH_TO_ROSETTA/main/source/bin/
```
**Note** that you need to replace PATH_TO_ROSETTA with your path to the Rosetta code.

This will create a directory named `relax_start_structure` and command files `ALL_RELAX_COMMANDS` and `RELAX_COMMAND_1`. `relax_start_structure` contains sub-directories for each relaxation run. There should be `nstructs` numbered sub-directories (here `nstructs` = 1, so there is only 1 numbered sub-directory).

&nbsp;&nbsp;&nbsp;&nbsp;**1.2**&nbsp;&nbsp;&nbsp;&nbsp; To run the relaxation calculations that we just set up in the previous step, type: 

```
source ALL_RELAX_COMMANDS
```

This will take a few minutes to run. The final relaxed structure is in `relax_start_structure/1/min_again_start_structure_wildtype_bound.pdb`. If we had specified `nstructs` > 1, then there would be a file named `min_again_start_structure_wildtype_bound.pdb` in each of the numbered directories in `relax_start_structure`.

&nbsp;&nbsp;&nbsp;&nbsp;**1.3**&nbsp;&nbsp;&nbsp;&nbsp; To get a ranked list of the lowest scoring relaxed structures, type:

```
python PATH_TO_ROSETTA/main/source/src/apps/public/rnp_ddg/get_lowest_scoring_relaxed_models.py --relax_dir relax_start_structure/
```

**Note:** you need to replace PATH_TO_ROSETTA with your path to the Rosetta code.  
This will print the following message to the screen:  

```
The lowest scoring 1 models:
Rank 0, Score -999.996: relax_start_structure//1/min_again_start_structure_wildtype_bound.pdb
```

(Here, this wasn't completely necessary since the relaxation was only performed once.)

## Step 2: Perform ddG calculations.

&nbsp;&nbsp;&nbsp;&nbsp;**2.1**&nbsp;&nbsp;&nbsp;&nbsp; Set up the ddG calculations. This will create a directory and all the necessary files for the ddG calculation runs. Type:

```
python PATH_TO_ROSETTA/main/source/src/apps/public/rnp_ddg/general_RNP_setup_script.py --low_res --tag demo_run --start_struct relax_start_structure//1/min_again_start_structure_wildtype_bound.pdb --seq_file mutant_list.txt --rosetta_prefix PATH_TO_ROSETTA/main/source/bin/ --Nreps 1 --protein_pack_reps 1
```

**Note:** you need to replace PATH_TO_ROSETTA with your path to the Rosetta code.

Let's walk through each of the options:  

`--low-res`: Use the low-res method of calculating RNA-protein relative binding affinities. This is recommended. Alternatively, you can specify `--med-res`, which will use `rna_denovo` to build mutant structures when they introduce a non-canonical RNA base pair (i.e. if the WT structure contains a G-C base pair and you specify a mutation to G-G).  

`--tag`: Here we can specify a string that will be used to create the name of the output directory.  

`--start_struct`: This is the relaxed starting structure that will be used to calculate energies of mutants. This will be considered the wildtype sequence. This should be one of the lowest-scoring structures that was printed to the screen by `get_lowest_scoring_relaxed_models.py`.  

`--seq_file`: This is a text file specifying the sequences for which we want to calculate relative binding affinities. One sequence should be specified per line. These can either be the full sequence of the complex (RNA and protein), or just the RNA sequence. If the protein sequence is not specified, then no mutations to the protein will be made.  

`--rosetta_prefix`: The path to the Rosetta executables. 

`--Nreps`: The number of times the mutation and subsequent relaxation of surrounding residues should be performed. **For actual runs, it is recommended that this is set to 10.**

`--protein_pack_reps`: The number of times the "unbound" protein structure should be repacked, to calculate the energy of the unbound protein. **For actual runs, it is recommended that this is set to 10.**  

   
This setup command will create a directory named `ddG_demo_run_low-res`. In this directory, there is a numbered subdirectory corresponding to each of the mutant sequences that were specified in the `mutant_list.txt` file. `0/` corresponds to the wildtype, `1/` corresponds to the first sequence listed in `mutant_list.txt`, `2/` corresponds to the second sequence listed in `mutant_list.txt`, etc. `ddG_demo_run_low/general_setup_settings.txt` lists the options that will be used for the run. `ddG_demo_run_low/ALL_COMMANDS` contains the actual command lines to run all of the ddG calculations. `ddG_demo_run_low/COMMAND_0`, `ddG_demo_run_low/COMMAND_1`, and `ddG_demo_run_low/COMMAND_2` contain the commands to run the ddG calculations for the wildtype, first, and second mutations, respectively. All of these commands are contained within `ddG_demo_run_low/ALL_COMMANDS`; these individual files are just useful if you're running a lot of mutants on a cluster -- each command can be run simultaneously on a different core.  


&nbsp;&nbsp;&nbsp;&nbsp;**2.2**&nbsp;&nbsp;&nbsp;&nbsp; Run the ddG calculations. Type:

```
source ddG_demo_run_low-res/ALL_COMMANDS
```

This will take several minutes to run.

&nbsp;&nbsp;&nbsp;&nbsp;**2.3**&nbsp;&nbsp;&nbsp;&nbsp; Get the ddG results. Type:

```
python PATH_TO_ROSETTA/main/source/src/apps/public/rnp_ddg/get_final_ddG_scores.py --run_dir ddG_demo_run_low-res/ --seq_file mutant_list.txt
```
**Note:** you need to replace PATH_TO_ROSETTA with your path to the Rosetta code.

`--run_dir` should specify the directory that the ddG calculations were run in (the one that was created by `general_RNP_setup_script.py`). `--seq_file` should list the same `seq_file` that was provided to `general_RNP_setup_script.py`.  
 
This will print the following to the screen:  

```
####################################
RESULTS:

WT ddG: 0.00 kcal/mol
ugaggcucaccca ddG: 2.39 kcal/mol
ugaggagcaccca ddG: 0.62 kcal/mol

####################################
Results are also written to the ddG_score.txt files in 
each mutant directory in ddG_demo_run_low-res/
The format is:
ddG dG complex_score protein_score rna_score
####################################
```

These calculations are not deterministic, so the actual numbers might differ slightly. **Normally, these calculations should be performed on the top 20 relaxed structures, then the final ddG values should be averaged over the 20 results.**  
As noted in the message above, the ddG results are also listed in the `ddG_score.txt` files in each of the run directories. The first column in a `ddG_score.txt` specifies the calculated ddG value in kcal/mol.

For reference, example output for this demo is provided in the `example_output` directory.   

## Additional information

See the Rosetta-Vienna ddG documentation [here](https://www.rosettacommons.org/docs/latest/application_documentation/rna/rnp-ddg).
