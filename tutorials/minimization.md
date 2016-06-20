# Minimization

## Introduction
Proteins are not static structures, but rather they undergo fluctuations in their conformations and exist as ensembles of states. 

Each snapshot of a protein's repertoire of conformations can be associated with an energy, where some conformations will have high energies and some will have low energies. In molecular modeling, it is usually desirable to find the global minimum (representing the lowest-energy conformation) of this energy function. This, however, is a very difficult task given the vast energy landscape that needs to be searched, so we'll settle for the next best thing: **a local minimum**.

Minimization is a sampling technique for the purpose of finding the nearest local minimum in the energy function given your starting structure's conformation and energy.

In Rosetta, we use the `minimize` executable to perform minimization on protein structures, which in general carries out a _gradient-based minimization_ to find the nearest local minimum in the energy function. 

All minimization algorithms choose a vector as the descent direction, determine a step along that vector, then choose a new direction and repeat. However, there are many different minimization algorithms that the user can choose that work under-the-hood of the `minimize` executable. In this tutorial, we will use `lbfgs_armijo_nonmonotone` algorithm, which is a multi-step algorithm that needs to only be called once to reach the local mimimum of a function (rather than 

## Goals
What are the goals of this tutorial?
In this tutorial, you will learn to use one of the many minimization algorithms, both with and without constraints.
* Learn to run the `minimize` and `minimize_with_cst` executables via the command line.
* Use MoveMaps to control the movement of residues during minimization.
* 


## How-To: Minimize
### Analyzing the input structure
Take a look at the crystal structure file `3hon.pdb` in your favorite PDB Viewer (eg. PyMOL). This structure was taken directly from the RCSB PDB and contains many physical _errors_. 

Let's get an initial score for this conformation. In your command line, type:

```bash
$> score.default.linuxgccrelease -database <database> -s 3hon.pdb
```

This command will output a `default.sc` file, which contains the Rosetta score for this conformation. 

[picture]

The first line of this file is the header line that tells us which columns correspond to which score terms.

On the second line of this file, we see that this conformation's score is total score is 240.074, and that it has a particularly high repulsive score, _fa_rep_, which means there are clashes between atoms in the structure, and a high dunbrack score, _fa_dun_, which means many of the rotamers in this structure are of low probability.

Let's try and fix these issues using the minimizer.

### Setting up the flags file

First, we will need to specify how to run the minimizer. Open the `minimizer_flags` file and analyze its contents:

 - -s 3hon.pdb
 - -run:min_type lbfgs_armijo_nonmonotone
 - -run:min_tolerance 0.001

The first flag specifies our input file, in this case, the crystal structure of 3hon.

The second flag specifies the type of minimization algorithm to use, in this case, lbfgs_armijo_nonmonotone.

The third flag specifies the convergence tolerance for the minimization algorithm. Rosetta has at least two kinds of "tolerance" for function minimization, "regular" (for lack of a better name) tolerance and absolute tolerance. "Regular" tolerance is _fractional_ tolerance for the _value_ of the function being minimized; i.e. a tolerance of 0.01 means the minimum function value found will be within 1% of the true minimum value. Absolute tolerance is specified without regard to the current function value; i.e. an absolute tolerance of 0.01 means that the minimum function value found will be equal to the actual minimum plus or minus 0.01, period. Minimizers use "regular" fractional tolerance by default. (If you would like to learn more about absolute tolerance, see here (somethingsomething)[]). In general, setting the tolerance to 0.01 is very loose, and it is recommended to specify a tolerance setting of something less than 0.01. Therefore, for this tutorial we have set the tolerance to 0.001.

### Running the minimization command

In your terminal window, run the minimization executable as follows:

```bash
$> minimization.default.linuxgccrelease @minimization_flags
```
If this executable runs with no errors and the terminal output ends with something like this:

[picture]

then you have successfully minimized the input structure. Check to make sure a `score.sc` file and a `3hon_0001.pdb` output structure have also been generated. Now let's analyze the output scores and structure.

### Analyzing the output

The `score.sc` file contains the NEW scores for the minimized structure. Open this file and compare the initial scores to the minimized scores:

| Score Term | Xtal Scores | Minimized Scores 	|
| ------------- | -----------: |  -----------: |
| total_score	| 240.074		| -40.93		|
| fa_atr		| -213.265		| -213.025		|
| fa_rep		| 131.638		| 22.807		|
| fa_sol		| 111.304		| 105.466		|
| fa_intra_rep	| 0.82			| 0.649			|
| fa_elec		| -9.08			| -11.103		|
| pro_close		| 15.694		| 0.434			|
| hbond_sr_bb	| -2.736		| -7.094		|
| hbond_lr_bb	| -7.547		| -13.808		|
| hbond_sc		| -0.757		| -1.842		|
| rama			| 4.53			| -2.326		|
| omega			| 0.404			| 3.677			|
| fa_dun		| 215.208		| 91.005		|
| p_aa_pp		| -3.773		| -12.833		|

Most of the score terms have gone down in value (which is good!). But how has the minimization affected the structure? Open the `3hon_0001.pdb` to see what has changed:

[picture]

The minimized structure has moved out of alignment with the native structure, so first let's align the two structures. If you are using PyMOL, type `align 3hon_0001, 3hon` and hit `Enter`.

[picture]

Now that the structure are aligned, notice that the last nine residues of the loop region have moved quite significantly. Furthermore, most of the minimized rotamers are no longer in their native conformations.

Sometimes it may be undesirable to allow large movements in the starting conformation. To this end, we can use minimization with constraints to minimize our input structure, while movements in certain atoms will be penalized by the score function.


## How-To: Minimization with Constraints

For this walk-through, we will use the same crystal structure as before, `3hon.pdb`, but this time we will apply harmonic coordinate constraints on the backbone heavy atoms of the nine C-terminal tail residues.

### Setting up the flags file and constraints file

We will need to add extra options to our flagsfile in order to tell Rosetta to read the constraints that specifies how to hold the atoms in place and to use a scorefunction that has non-zero weights for the constraint terms. The new flags file looks like this:

- -s 3hon.pdb
- -run:min_type lbfgs_armijo_nonmonotone
- -run:min_toleranace 0.001
- -constraints:cst_file cstfile
- -score:weights talaris2014_cst
- out:suffix _minwithcsts

The flag `constraints:cst_file` specifies the name of the constraints file to use.

The flag `score:weights` specifies which set of weights to use for the energy function. In this case, we have specified the talaris2014_cst weights file, which is identical to talaris2014 except for the addition of the constraint terms _chainbreak_, _coordinate_constraint_, _atom_pair_constraint_, _angle_constraint_, _dihedral_constraint_, _res_type_constraint_, each with a weight of 1.0. 

The flag `out:suffix` will prevent us from overwriting our previous output by appending the suffix "_minwithcsts" to the end of our output file names.


### Running the minimization command

Run the minimization executable in the same way as before, but now with the new flags file:

```bash
$> minimize.cc @minwithcsts_flags
```

This time, look for a few lines coming from the `core.scoring.constraints` tracers in the log output. They should look similar to this:

[picture]

These lines indicate that our constraints were read and applied to the pose. If the executable ended without errors and generated a `3hon_minwithcsts_0001.pdb` file as well as a `score_minwithcsts.sc` file, then you have succesfully run minimization with coordinate constraints.

### Analyzing the output

As before, the `score_minwithcsts.sc` file contains the score of the minimized structure. Open this file and compare its energy to the crystal structure's score as

| Score Term	  | Xtal Scores | Minimized Scores | Minimized With CST Scores |
| -----------  | ----------: |  --------------: | ------------------------: |
| total_score	 | 240.074 | -40.93		         | -28.538                   |
| fa_atr		     | -213.265 | -213.025		       | -207.931 |
| fa_rep		     | 131.638 | 22.807		         | 21.68 |
| fa_sol		     | 111.304 | 105.466		        | 104.359 |
| fa_intra_rep | 0.82			     | 0.649			         | 0.714 |
| fa_elec		    | -9.08			    | -11.103		        | -12.245 |
| pro_close		  | 15.694		    | 0.434			         | 0.781 |
| hbond_sr_bb	 | -2.736		    | -7.094		         | -5.19 |
| hbond_lr_bb	 | -7.547		    | -13.808		        | -13.425 |
| hbond_sc		   | -0.757		    | -1.842		         | -1.239 |
| rama			      | 4.53			     | -2.326		         | -3.515 |
| omega			     | 0.404			    | 3.677			         | 3.989 |
| fa_dun		     | 215.208		   | 91.005		         | 96.731 |
| p_aa_pp		    | -3.773		    | -12.833		        | -10.393 |
| coordinate_constraint | | | 0.727 |

As before, most of the new scores from the minimized structure with constraints are lower than those in the crystal structure. Notice also the addition of the coordinate constraint term to the list of energy terms for the newly minimized structure. Comparing the minimized structure to the minimized-with-csts structure, we see an increase in total energy cause predominantly by differences in the fa_atr and fa_dun terms. Let's look at the structure now to see if those differences make sense.

Open the `3hon_minwithcsts_0001.pdb` file and compare it to the crystal structure.

[picture]

This time we see little to no movement in the nine C-terminal residues.

There is another way to prevent motion all or certain residues of the protein structure, and this is by using a MoveMap.

## How-To: Minimization with a MoveMap

### The MoveMap (https://www.rosettacommons.org/docs/wiki/rosetta_basics/structural_concepts/Rosetta-overview#scoring_movemap)
Certain protocols accept a user-specific move map that tells the algorithm what torsion angles and rigid-body degrees of freedom (DOFs) are allowed to move. For example, one may not want to move highly-conserved sidechains in modeling applications, or one may want to preserve certain interactions in design applications.

#### Description of the file format
Each line in the movemap file identifies a jump, residue, or residue range, followed by the allowed degrees of freedom. These entities may be specified as follows:
  RESIDUE 


### Setting up the flagsfile and movemap file

### Running the minimizer with a movemap

### Analysing output



######Notes
Molecular modeling of proteins most often requires the resolution of sterically clashing atoms as a form of quality control to ensure that our input structures are physically accurate.

Molecular models with physically overlapping (or unreasonably close) atoms are generally bad.
