# Generalized Kinematic Closure Tutorial 1:
# Basic loop closure with GeneralizedKIC
======================================

KEYWORDS: LOOPS SCRIPTING_INTERFACES

Tutorial by Vikram K. Mulligan (vmullig@uw.edu).  Created on 28 March 2017 for the Baker lab Rosetta Tutorial Series.

[[_TOC_]]

## Goals

At the end of this tutorial, you will understand:

- What the kinematic closure algorithm is, and what problem it solves
- How to use the PeptideStubMover to add loop residues to a pose lacking a loop
- How to use the DeclareBond mover to create a bond
- How to use the GeneralizedKIC mover to sample loop conformations
- How to set up GeneralizedKIC perturbers, filters, and selectors
- How to instruct GeneralizedKIC to close a peptide bond

## Introduction to Kinematic Closure

Kinematic closure algorithms were originally developed for the robotics field to solve the problem of determining the necessary joint angles that would place a robot's hand or foot in a desired place.  We have adapted them for use within the Rosetta software suite to sample conformations of chains of atoms with well-defined start and end points.

A molecular kinematic closure problem may be described as follows: given a covalently-contiguous chain of atoms within a molecule, with covalent linkages at the ends of the chain fixing the start and end of the chain, what possible conformations maintain the integrity of bond length, bond angle, and dihedral angle restrictions within the chain?  To solve such a problem, we divide the chain of atoms into two segments, and define "pivot points" at the start of the chain (the first pivot), the end of the chain (the last pivot), and the breakpoint between the two segments (the middle pivot).

Having done this, the degrees of freedom within the two segments may be held fixed, randomized, perturbed, or otherwise altered as one sees fit.  These degrees of freedom include bond lengths, bond angles, and dihedral angles.  Whatever one does to these degrees of freedom, one ends up with two segments that still have well-defined rigid body transforms from the first pivot to the middle pivot (in the first segment), and from the middle pivot to the last pivot (in the second segment).  It is then possible to solve a system of equations for the six torsion angles adjacent to the three pivots in order to keep the system closed.  The matrix math that gives rise to the solution(s) is extremely fast as compared to alternative loop closure methods (which typically rely on iterative gradient-descent minimization); however, for a given system, this step may yield anywhere from 0 to 16 solutions.  It then becomes necessary to choose a solution for downstream molecular design or conformational refinement.

The GeneralizedKIC mover in Rosetta gives a user full control over pre-closure sampling, post-closure filtering, and selection of a closure solution.  It is fully accessible to the RosettaScripts scripting language, and interfaces nicely with other Rosetta movers and filters, allowing arbitrary protocols to be carried out on closure solutions before choosing a final solution.  It also allows closure of chains that do not consist solely of polypeptide backbones: that is, it is fully compatible with closure of atomic chains that run through disulfide bonds, arbitrary side-chain cross-links or cross-linkers, and non-canonical backbones.

## A Note on Using GeneralizedKIC

A loop that is open may be thought of as a continuous loop containing a bond that is badly stretched, and which likely has very strange bond angles and torsion angles at the cutpoint.  GeneralizedKIC can close an open loop by using perturbations that set the bond length, bond angles, and, possibly, the torsion angle of the cutpoint to reasonable values prior to solving for pivot torsion values.

## Exercise 1: Building and Closing a Polypeptide Loop Using RosettaScripts

### Inputs

For this exercise, we will be using an NMR structure of an artificial mini-protein designed by Dr. Chris Bahl (PDB ID 2ND2).  This mini-protein is a 44-residue 3-helix bundle.  For the purposes of this tutorial, the structure has been stripped of its amino acid sequence (_i.e._ it has been mutated to poly-glycine), and the loop connecting the second and third helices has been deleted.  This is meant to simulate many common design cases, in which one might arrange secondary structure elements first and build loops later (_e.g._ in the case of parametric design approaches), as well as certain structure prediction cases, in which one might wish to model loops that are missing in crystal structures.  We will rebuild this loop and sample its possible conformations.

**The input structure, an edited version of PDB structure 2ND2 (`2ND2_state1_glyonly_loop_removed.pdb`):**
![The input structure](images/Example1_input_structure.png)

Additionally, we will use the following Rosetta flags file.  Briefly, this instructs Rosetta to run the input script 10 times to produce 10 sampled loop conformations, to use the `beta_nov15` score function, and to include all chemical bonds in the output PDB files (which can be convenient when debugging bad geometry, since bonds are drawn even if bonded atoms are too far apart).

**File `rosetta.flags`:**
```
-nstruct 10
-beta_nov15
-in:file:s inputs/2ND2_state1_glyonly_loop_removed.pdb
-in:file:fullatom
-write_all_connect_info
-parser:protocol xml/exercise1.xml
-jd2:failed_job_exception false
-mute protocols.generalized_kinematic_closure.filter.GeneralizedKICfilter core.chemical.AtomICoor core.conformation.Residue
```

### Step 1: Building loop geometry

The GeneralizedKIC mover is only capable of sampling conformations of existing geometry.  It can neither add amino acid residues to a pose, nor create new bonds between residues.  For this reason, we must use the [PeptideStubMover](https://www.rosettacommons.org/docs/latest/PeptideStubMover) to build the new loop, and the [DeclareBond mover](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/DeclareBond) to add a chemical bond across the loop cutpoint.

Run the `rosetta_scripts` application with no commandline options to create a template script, or copy and paste the example below:

```xml
<ROSETTASCRIPTS>
	<SCOREFXNS>
	</SCOREFXNS>
	<RESIDUE_SELECTORS>
	</RESIDUE_SELECTORS>
	<TASKOPERATIONS>
	</TASKOPERATIONS>
	<FILTERS>
	</FILTERS>
	<MOVERS>
	</MOVERS>
	<APPLY_TO_POSE>
	</APPLY_TO_POSE>
	<PROTOCOLS>
	</PROTOCOLS>
	<OUTPUT />
</ROSETTASCRIPTS>
```

Now let's add a [PeptideStubMover](https://www.rosettacommons.org/docs/latest/PeptideStubMover) in the `<MOVERS>` section.  We will append three residues to the end of the second helix (residue 28), and prepend two residues to the start of the third helix (which was residue 29, but which becomes residue 32 after appending three residues).  Note that the `Insert` command is used instead of the `Append` command because the added residues are in the middle of the sequence.  We'll use a poly-alanine sequence for sampling, but will cheat a little bit just for the purposes of this tutorial by keeping a glycine at the second loop position, since this is present in the original structure.

```xml
<PeptideStubMover name="add_loop_residues" >
	<Insert anchor_rsd="28" resname="ALA" />
	<Insert anchor_rsd="29" resname="GLY" />
	<Insert anchor_rsd="30" resname="ALA" />
	<Prepend anchor_rsd="32" resname="ALA" />
	<Prepend anchor_rsd="32" resname="ALA" />
</PeptideStubMover>

```

Be sure to add the above mover to the `<PROTOCOLS>` section as well (`<Add mover="add_loop_residues" />`).

We now have a five-residue loop, albeit one with several problems.  First, there is no bond between residues 30 and 31.  Second, these residues are too far apart in space to form a bond.  And third, all dihedral angles in the new loop are set to 0 degrees, including omega angles (_i.e._ all peptide bonds are _cis_ instead of _trans_).  To solve the first problem, we will add a [DeclareBond mover](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/DeclareBond) to the `<MOVERS>` section of the script, telling Rosetta that there ought to be a chemical bond between the C atom of residue 30 and the N atom of residue 31.  Note that this does not move any geometry, but it does mean that the bond is present in the output PDB file if we run the script at this point.

```xml
<DeclareBond name="new_bond" atom1="C" atom2="N" res1="31" res2="32" />
```

As before, this must also be added to the `<PROTOCOLS>` section (`<Add mover="new_bond" />`);

### Step 2: Preparing for sampling

We will create a GeneralizedKIC mover for sampling in a moment.  Before we do so, though, we want to mutate two flanking residues, which will also be sampled as part of the loop, to alanine.  We do this because the sampling that we will use is biased by the Ramachandran preferences of the amino acid type.  Alanine is a good residue to use to represent the generic L-alpha amino acid, but glycine (which is currently present) is not because it has as much preference for the positive-phi region of Ramachandran space as for the negative.  (Indeed, our glycine Ramachandran tables are biased to _favour_ the positive-phi region, since glycine is disproportionately found in the positive-phi region in the PDB structures used to train our scorefunction.)  Let us add two [MutateResidue movers](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/MutateResidueMover) to change the residue identities at this position.  In the `<MOVERS>` section add:

```xml
<MutateResidue name="mut1" target="28" new_res="ALA" />
<MutateResidue name="mut2" target="34" new_res="ALA" />
```

As before, add these two movers to the `<PROTOCOLS>` section after the movers already added.  At this point, the `<PROTOCOLS>` section should look like this:

```xml
<PROTOCOLS>
	<Add mover="add_loop_residues" />
	<Add mover="new_bond" />
	<Add mover="mut1" />
	<Add mover="mut2" />
</PROTOCOLS>
```

### Step 3: Initial GeneralizedKIC setup

We're now ready to add the GeneralizedKIC mover that will close the gap in the loop and sample loop conformations.  In the movers section, add the following.

```xml
<GeneralizedKIC name="genkic" />
```

Add this mover to the end of the `<PROTOCOLS`> section, too, as before.

GeneralizedKIC gives the user many, many options to control loop sampling, filtering, and solution selection.  Typically, the first option that one wants to set is the number of attempts that the mover will make to find a closed solution.  Each attempt consists of perturbing loop degrees of freedom (in a manner that we will define), solving for closed solutions, and applying GeneralizedKIC filters.  Because the KIC algorithm is so fast, one can easily attempt hundreds of solutions per second.  For our purposes, let's set the number of attempts to 5000 by modifying our GeneralizedKIC block as follows:

```xml
<GeneralizedKIC name="genkic" closure_attempts="5000" />
```

Each attempt could yield anywhere from 0 to 16 solutions.  By default, every solution found will be stored until we've made the specified number of attempts (in our case 5000).  This could be far too many solutions, though.  It makes more sense to stop looking for solutions after we've found a small number.  That number could be as low as 1 (_i.e._ GeneralizedKIC stops as soon as a solution is found), but for our purposes, let's set that number at 5:


```xml
<GeneralizedKIC name="genkic" closure_attempts="5000" stop_when_n_solutions_found="5" />
```

Even if we had set this numer at 1, a single attempt might have yielded up to 16 solutions.  We always need to tell GeneralizedKIC how to pick a single solution from among the solutions.  Here, we'll choose our solution by energy -- but there is an important caveat.  Since we are sampling backbone conformations, with no consideration of side-chains, we should use a scoring function that consists primarily of backbone-only terms to pick the best solution.  (Later we will see how we can apply an arbitrary mover -- _e.g._ a full repack and minimization -- to every solution, in which case it might make sense to use the full Rosetta scoring function to pick the best solution).  Let's set up a backbone-only scoring function using weights from the `beta_nov15` scoring function.  In the `<SCOREFXNS>` section of your script, add the following:

```xml
<ScoreFunction name="bnv" weights="beta_nov15.wts" />
<ScoreFunction name="bb_only" weights="empty.wts" >
	<Reweight scoretype="fa_rep" weight="0.1" />
	<Reweight scoretype="fa_atr" weight="0.2" />
	<Reweight scoretype="hbond_sr_bb" weight="2.0" />
	<Reweight scoretype="hbond_lr_bb" weight="2.0" />
	<Reweight scoretype="rama_prepro" weight="0.45" />
	<Reweight scoretype="omega" weight="0.4" />
	<Reweight scoretype="p_aa_pp" weight="0.6" />
</ScoreFunction>

```

The "bnv" scoring function is now the full `beta_nov15` scoring function, which may be useful if we later add sidechain refinement to our protocol (which we will do in the third GeneralizedKIC tutorial).  For now, we will use the "bb_only" scoring function, adding two more options to our already-declared GeneralizedKIC mover to tell it how to pick a solution:

```xml
<GeneralizedKIC name="genkic" closure_attempts="5000" stop_when_n_solutions_found="5"
	selector="lowest_energy_selector" selector_scorefunction="bb_only"
/>
```

## Conclusion

**TODO**

## Further Reading

Bhardwaj G, Mulligan VK, Bahl CD, Gilmore JM, Harvey PJ, Cheneval O, Buchko GW, Pulavarti SV, Kaas Q, Eletsky A, Huang PS, Johnsen WA, Greisen PJ, Rocklin GJ, Song Y, Linsky TW, Watkins A, Rettie SA, Xu X, Carter LP, Bonneau R, Olson JM, Coutsias E, Correnti CE, Szyperski T, Craik DJ, Baker D.  (2016).  Accurate de novo design of hyperstable constrained peptides.  _Nature_ 538(7625):329-335.

Mandell DJ, Coutsias EA, Kortemme T. (2009).  Sub-angstrom accuracy in protein loop reconstruction by robotics-inspired conformational sampling.  _Nat. Methods_ 6(8):551-2.

Coutsias EA, Seok C, Jacobson MP, Dill KA.  (2004).  A kinematic view of loop closure.  _J. Comput. Chem._ 25(4):510-28.

[GeneralizedKIC documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/composite_protocols/generalized_kic/GeneralizedKIC)

