# Minimization

## Introduction
Proteins are not static structures, but rather they undergo fluctuations in their conformations and exist as an ensemble of states. 

Each snapshot of a protein in its repertoire of conformations can be associated with an energy, where some conformations will have high energies and some will have low energies. In molecular modeling, it is usually desirable to find the global minimum (representing the lowest-energy conformation) of this energy function. This, however, is a very difficult task given the vast energy landscape that needs to be searched, so we'll settle for the next best thing: **a local minimum**.

Minimization is a sampling technique for the purpose of finding the nearest local minimum in the energy function given a starting structure's conformation and energy.

In Rosetta, we use the `minimize` executable to perform minimization on protein structures, which in general carries out a _gradient-based minimization_ to find the nearest local minimum in the energy function. There are many different [minimization algorithms](https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/minimization-overview#flavors-of-minimization-in-rosetta) that the `minimize` application can use, but essentially all minimization algorithms choose a vector as the descent direction, determine a step along that vector, then choose a new direction and repeat. In this tutorial, we will use the `lbfgs_armijo_nonmonotone` algorithm, which is a multi-step algorithm that needs to only be called once to reach the local mimimum of a function (rather than invoking repeated iterations to reach convergence).

## Goals
In this tutorial, you will learn to use one of the many minimization algorithms in three different ways, first by allowing all residues in the structure to move during the minimation and then by using two methods that allow only a subset of the degrees of freedom (DOFs) to be sampled. Specifically, we will:
* Learn to run the `score` and `minimize` executables via the command line.
* Create and utilize a **constraints file** to prevent the movement of CA atoms from ocurring in part of the structure.
* Create and utilize a **move map** to control the degrees of freedom allowed to move during the minimization.
* Compare and analyze the output files and structures from each type of minimization.

## How-To: Minimize

### Analyzing the input structure

Take a look at the crystal structure file `3hon.pdb` in your favorite PDB Viewer (e.g. PyMOL). This structure was taken directly from the [RCSB PDB](http://www.rcsb.org/pdb/explore.do?structureId=3HON) and contains many physical "errors". 

Let's get an initial score for this conformation. In your command line, type:

```bash
$> <path/to/Rosetta/bin/>score.default.linuxgccrelease -s 3hon.pdb
```

This command will output a `default.sc` file, which contains the Rosetta score for this conformation. Inside this file, we see:

![default scorefile](https://github.com/RosettaCommons/demos/blob/XRW2016_kmb/tutorials/default_score.png)

The first line of this file is the header line that tells us which columns correspond to which score terms.

On the second line of this file, we see that the crystal conformation has a total Rosetta score of 240.074, and that it has a particularly high repulsive score, _fa_rep_, which means there are clashes between atoms in the structure, and a high dunbrack score, _fa_dun_, which means many of the rotamers in this structure are of low probability.

Let's try and fix these issues using the minimizer.

### Setting up the flags file

First, we will need to specify how to run the minimizer. Open the `minimizer_flags` file and analyze its contents:

```
 -s 3hon.pdb
 -run:min_type lbfgs_armijo_nonmonotone
 -run:min_tolerance 0.001
```

The first flag, `s`, specifies our input file, in this case, the crystal structure of 3hon.

The second flag, `run:min_type`, specifies the type of minimization algorithm to use, in this case, lbfgs_armijo_nonmonotone.

The third flag, `run:min_tolerance`, specifies the convergence tolerance for the minimization algorithm. Rosetta has at least two kinds of "tolerance" for function minimization, "regular" (for lack of a better name) tolerance and absolute tolerance. "Regular" tolerance is _fractional_ tolerance for the _value_ of the function being minimized; i.e. a tolerance of 0.01 means the minimum function value found will be within 1% of the true minimum value. Absolute tolerance is specified without regard to the current function value; i.e. an absolute tolerance of 0.01 means that the minimum function value found will be equal to the actual minimum plus or minus 0.01, period. Minimizers use "regular" fractional tolerance by default. (Click [here](https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/minimization-overview#the-meaning-of-tolerance) for more information on absolute tolerance). In general, setting the fractional tolerance to 0.01 is very loose, and so it is recommended to specify a tolerance setting of something less than 0.01. Therefore, for this tutorial we have set the tolerance to 0.001.

### Running the minimization command

In your terminal window, run the minimization executable by typing,

```bash
$> <path/to/Rosetta/bin/>minimize.default.linuxgccrelease @minimizer_flags
```
If this executable runs with no errors and the terminal output ends with something like this,

![min log](https://github.com/RosettaCommons/demos/blob/XRW2016_kmb/tutorials/min.png)

then you have successfully minimized the input structure. Check to make sure a `score.sc` file and a `3hon_0001.pdb` output structure have also been generated. Now let's analyze the output scores and structure.

### Analyzing the output

The `score.sc` file contains the new scores for the minimized structure. Open this file and compare the initial crystal (xtal) structure scores to the minimized scores:

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

Most of the score terms have gone down in value (which is good!). But how has the minimization affected the structure? Open the `3hon_0001.pdb` to see what has changed.

![3hon minimized](https://github.com/RosettaCommons/demos/blob/XRW2016_kmb/tutorials/3hon_min_on_xtal.png)

The minimized structure (in green) has moved out of alignment with the native structure (in cyan), so first let's align the two structures. If you are using PyMOL, type `align 3hon_0001, 3hon` and hit `Enter`.

![3hon minimized_aligned](https://github.com/RosettaCommons/demos/blob/XRW2016_kmb/tutorials/3hon_minaligned_on_xtal.png)

Now that the structures are aligned, notice that the last nine residues of the loop region have moved quite significantly from their original conformation. Furthermore, most of the minimized rotamers are no longer in their native conformations.

![3hon minimized_aligned sticks](https://github.com/RosettaCommons/demos/blob/XRW2016_kmb/tutorials/3hon_minaligned_on_xtal_sticks.png)

_Something something about what it all means blah blah._

Sometimes it may be undesirable to allow such large movements in the starting conformation. To this end, we can use minimization with constraints to minimize our input structure in which movements of certain atoms will be penalized by the score function.

## How-To: Minimization with Constraints

For this section of the tutorial, we will use the same crystal structure as before, `3hon.pdb`, but this time we will apply harmonic coordinate constraints (csts) on the backbone heavy atoms of the nine C-terminal tail residues.

### Setting up the flags file and constraints file

We will need to add extra options to our flagsfile in order to tell Rosetta to read the constraints that specify how to hold the atoms in place. Furthermore, we need to be sure we are using a scorefunction that has non-zero weights for the constraint terms. The new flags file, `minwithcsts_flags`, looks like this:

```
-s 3hon.pdb
-run:min_type lbfgs_armijo_nonmonotone
-run:min_toleranace 0.001
-constraints:cst_file cstfile
-score:weights talaris2014_cst
-out:suffix _minwithcsts
```

The flag `constraints:cst_file` specifies the name of the constraints file to use.

The flag `score:weights` specifies which set of weights to use for the energy function. In this case, we have specified the talaris2014_cst weights file, which is identical to talaris2014 except for the addition of the constraint terms _chainbreak_, _coordinate_constraint_, _atom_pair_constraint_, _angle_constraint_, _dihedral_constraint_, and _res_type_constraint_, each with a weight of 1.0. 

The flag `out:suffix` will prevent us from overwriting our previous output by appending the suffix "_minwithcsts" to the end of our new output file names.

### Running the minimization command

Run the minimization executable in the same way as before, but now with the new flags file:

```bash
$> <path/to/Rosetta/bin/>minimize.cc @minwithcsts_flags
```

This time, look for a few lines coming from the `core.scoring.constraints` tracers in the log output. They should look similar to this:

![mincsts log](https://github.com/RosettaCommons/demos/blob/XRW2016_kmb/tutorials/minwithcsts.png)

These lines indicate that our constraints were read and applied to the pose. If the executable ended without errors and generated a `3hon_minwithcsts_0001.pdb` file as well as a `score_minwithcsts.sc` file, then you have succesfully run minimization with coordinate constraints.

### Analyzing the output

As before, the `score_minwithcsts.sc` file contains the score of the minimized structure subject to constraints. Open this file and compare its energy to the crystal structure's and the structure generate from constraints-free minimization:

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

Again, most of the new scores from the minimized-with-constraints structure are lower than those in the crystal structure. Notice also the addition of the coordinate constraint term to the list of energy terms for the newly minimized structure. 

Comparing the minimized structure to the minimized-with-csts structure, we see an increase in total energy caused predominantly by differences in the fa_atr and fa_dun terms. (Why?) 

Now lets take a look at the minimized-with-csts structure to see how it compares to our previous structures.

Opening the `3hon_minwithcsts_0001.pdb` file and comparing it to the crystal structure,

![3hon mincsts](https://github.com/RosettaCommons/demos/blob/XRW2016_kmb/tutorials/3hon_minwithcsts_onxtal.png)

we immediately see little to no movement in the nine C-terminal residues.

(Insert nice segue here: There is another way to prevent motion in all or a subset of residues of the protein structure and this is by using a MoveMap.)

## How-To: Minimization with a MoveMap

### The [MoveMap](https://www.rosettacommons.org/docs/wiki/rosetta_basics/structural_concepts/Rosetta-overview#scoring_movemap)
Certain protocols accept a user-defined move map file that tells the algorithm what torsion angles and rigid-body degrees of freedom (DOFs) are allowed to move. For example, one may not want to move highly-conserved sidechains in modeling applications, or one may want to preserve certain interactions in design applications.

In the context of the minimizer, a move map allows the user to specify if the backbone (BB) torsions angles (phi, psi) or the sidechain torsions angles (CHI) are allowed to be moved during the minimization of the energy function. In addition, if the input structure has more than one chain (separated by one or more JUMPS), the move map can also specify if rigid-body movements between the different chains are allowed.

##### Caveat: Even if a residue's backbone and sidechain torsion movements are turned off in a move map, its relative position with respect to other residues may still change depending on the motion of residues upstream in the FoldTree.

#### Description of the move map file format
Each line in the move map file identifies a jump, residue, or residue range, followed by the allowed degrees of freedom. These entities may be specified as follows:
```
RESIDUE <#> <BB/CHI/BBCHI/NO>         # a single residue <#> followed by a single option (pick one of the four)
RESIDUE <#1> <#2> <BB/CHI/BBCHI/NO>   # a range of residues from <#1> to <#2> followed by a single option (again, pick one)
JUMP <#> <YES/NO>                     # a jump <#> followed by a single option
```

For example,
```
RESIDUE 28 BB        # allows BB movements at residue 28
RESIDUE 32 48 BBCHI  # allows BB and sidechain CHI movements from residue 32 to residue 48
JUMP 1 YES           # allows rigid-body movements between the structures separated by jump 1
```
Is is also possible to choose all residues or all jumps with the a `*` symbol:
```
RESIDUE *  CHI    # allows sidechain movements at all residues
JUMP * YES        # allows rigid body movements between all structures separated by jumps
```

If a residue appears more than once, the last appearance in the file determines the movement (i.e. move map lines are NOT additive). For example, the movemap speficied here
```
RESIDUE * CHI
RESIDUE * BB
```
will only allow backbone (BB) movements for all residues and will disallow sidechain (CHI) movements for all residues, which is probably not what the user meant.

##### Note: If a residue or jump is not specified in the move map, it will revert to the default behavior, which is protocol-specific.

### Setting up the flagsfile and movemap file

For the next minimization walkthrough, we will need to add one option to our flags file to tell Rosetta to read the move map file. The contents of the new flags file, `minwithmm_flags`, should look like this:

```
-s 3hon.pdb
-run:min_type lbfgs_armijo_nonmonotone
-run:min_toleranace 0.001
-movemap movemapfile
-out:suffix _minwithmm
```

The flag `movemap` specifies the movemap file to apply to the pose.

Our move map file `movemapfile` (expanded below) tells the minimizer to first set all backbone (BB) and sidechain (CHI) torsion angles to movable. Then, it reverts residues 47 through 55 (in pose numbering) to be fixed.
```
RESIDUE * BBCHI
RESIDUE 47 55 NO
```

### Running the minimizer with a movemap

Run the minimization executable in the same way as before, but now with the new flags file:

```bash
$> minimize.cc @minwithmm_flags
```

The log output should be similar to the first minimization protocol except for a difference in score. You should also have a `3hon_minwithmm_0001.pdb` file and a `score_minwithmm.sc` file after running this command.

### Analyzing the output

Let's compare the score of the minimize-with-movemap structure to our previous crystal, minimized, and minimized-with-csts structures:

| Score Term	  | Xtal Scores | Minimized Scores | Minimized With CST Scores | Minimized with MoveMap |
| -----------  | ----------: |  --------------: | ------------------------: | ---------------------: |
| total_score	 | 240.074     | -40.93		         | -28.538                   | -4.217
| fa_atr		     | -213.265    | -213.025		       | -207.931                  | -212.704
| fa_rep		     | 131.638     | 22.807		         | 21.68                     | 23.41
| fa_sol		     | 111.304     | 105.466		        | 104.359                   | 105.753
| fa_intra_rep | 0.82			     | 0.649			         | 0.714                     | 0.658
| fa_elec		    | -9.08			    | -11.103		        | -12.245                   | -10.932
| pro_close		  | 15.694		    | 0.434			         | 0.781                     | 12.29
| hbond_sr_bb	 | -2.736		    | -7.094		         | -5.19                     | -7.069
| hbond_lr_bb	 | -7.547		    | -13.808		        | -13.425                   | -13.833
| hbond_sc		   | -0.757		    | -1.842		         | -1.239                    | -1.82
| rama			      | 4.53			     | -2.326		         | -3.515                    | -1.645
| omega			     | 0.404			    | 3.677			         | 3.989                     | 3.554
| fa_dun		     | 215.208		   | 91.005		         | 96.731                    | 113.143
| p_aa_pp		    | -3.773		    | -12.833		        | -10.393                   | -11.999
| coordinate_constraint | | | 0.727 | |

The total score of the minimized-with-movemap structure is lower still than the crystal structure but has the highest energy of the minimized structures, caused primarily by differences in the _pro_close_ term, which is an energy associated with the proline ring closure, and the _fa_dun_ term.

Let's open the minimized-with-movemap structure and compare it visually to the crystal structure.

![3hon minmm](https://github.com/RosettaCommons/demos/blob/XRW2016_kmb/tutorials/3hon_minwithmm_xtal.png)


After aligning the minimized-with-movemap structure to the crystal structure as a whole, it appears that the nine C-terminal residues have moved. However, when we align only the nine C-terminal residues against each other, it becomes clear that the minimizer has not changed the backbone or sidechain angles:

![3hon minmm](https://github.com/RosettaCommons/demos/blob/XRW2016_kmb/tutorials/3hon_minwithmm_9resicterm.png)

#### This is an important point to reiterate: The movemap can prevent internal geometries from changing, but not necessarily the global position of the residues.

## Summary

In this tutorial, we learned how to run the `minimize` executable in three ways: 
  1. By using the default behavior, which allowed movement in all backbone and chi angles,
  2. By applying coordinate constraints to the input structure, which disallowed global movements a subset of atoms,
  3. And by applying a MoveMap, which disallowed local movements in the backbone and sidechain torsion angles.
