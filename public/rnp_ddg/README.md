Calculate relative binding affinities for an RNA-protein complex
=====================================================================================

KEYWORDS: NUCLEIC_ACIDS RNA 

Written in April 2018 by Kalli Kappel (kappel at stanford dot edu).  

This demo shows how to use the Rosetta-Vienna ddG method to calculate relative binding affinites for an RNA-protein complex.  


## Setting up the demo:  

## Brief explanation of input files:  

## Step 1: Relax the starting structure.

```
python PATH_TO_ROSETTA/main/source/src/apps/public/rnp_ddg/relax_starting_structure.py --start_struct 1ZDH_chainA_B_R.pdb --nstructs 1 --rosetta_prefix PATH_TO_ROSETTA/main/source/bin/
```

```
source ALL_RELAX_COMMANDS
```

## Step 2: Perform ddG calculations.

For reference, example output is provided in the `example_output` directory.   

## Additional information

See the Rosetta-Vienna ddG documentation [here]().
