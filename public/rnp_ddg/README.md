Calculate relative binding affinities for a RNA-protein complex
=====================================================================================

KEYWORDS: NUCLEIC_ACIDS RNA 

Written in April 2018 by Kalli Kappel (kappel at stanford dot edu).  

This demo shows how to use the Rosetta-Vienna ddG method to calculate relative binding affinites for an RNA-protein complex.  


## Setting up the demo:  

## Brief explanation of input files:  

## Step 1: Relax the starting structure.

First, we need to "relax" the starting structure into the Rosetta score function that we'll be using for the ddG calculations.  

1. Set up the relaxation run. We need to provide the starting struture (this will be our "wildtype" -- mutant binding affinities will be calculated relative to the sequence in this PDB file). Here, we're only relaxing the structure once (--nstructs 1), but for actual runs it is recommended that the relaxation is performed 100 times (--nstructs 100). Type:  

```
python PATH_TO_ROSETTA/main/source/src/apps/public/rnp_ddg/relax_starting_structure.py --start_struct start_structure.pdb --nstructs 1 --rosetta_prefix PATH_TO_ROSETTA/main/source/bin/
```
*Note* that you need to replace PATH_TO_ROSETTA with your path to the Rosetta code.

This will create a directory named `relax_start_structure` and command files `ALL_RELAX_COMMANDS` and `RELAX_COMMAND_1`. `relax_start_structure` contains sub-directories for each relaxation run. There should be `nstructs` numbered sub-directories (here `nstructs` = 1, so there is only 1 numbered sub-directory).

2. To run the relaxation calculations that we just set up in the previous step, type: 

```
source ALL_RELAX_COMMANDS
```

The final relaxed structure is in `relax_start_structure/1/min_again_start_structure_wildtype_bound.pdb`. If we had specified `nstructs` > 1, then there would be a file named `min_again_start_structure_wildtype_bound.pdb` in each of the numbered directories in `relax_start_structure`.

3. To get a ranked list of the lowest scoring relaxed structures, type:

```
python PATH_TO_ROSETTA/main/source/src/apps/public/rnp_ddg/get_lowest_scoring_relaxed_models.py --relax_dir relax_start_structure/
```

This will print the following message to the screen:  

```
The lowest scoring 1 models:
Rank 0, Score -999.996: relax_start_structure//1/min_again_start_structure_wildtype_bound.pdb
```

(Here, this wasn't completely necessary since the relaxation was only performed once.)

## Step 2: Perform ddG calculations.

For reference, example output is provided in the `example_output` directory.   

## Additional information

See the Rosetta-Vienna ddG documentation [here]().
