# Generalized Kinematic Closure Tutorial 2:
# Perturbing loops from a starting conformation with GeneralizedKIC
======================================

KEYWORDS: LOOPS SCRIPTING_INTERFACES

Tutorial by Vikram K. Mulligan (vmullig@uw.edu).  Created on 28 March 2017 for the Baker lab Rosetta Tutorial Series.

[[_TOC_]]

## Goals

At the end of this tutorial, you will understand:

- How to use the GeneralizedKIC mover to sample small perturbations of an existing loop.
- How to use the GeneralizedKIC mover with the GenericMonteCarlo mover to perform a Monte Carlo search of loop conformations.

It is highly recommended that you complete the [first tutorial](generalized_kinematic_closure_1.md) before proceeding.

## Perturbing loops with GeneralizedKIC

In the [first tutorial](generalized_kinematic_closure_1.md), we saw how one can construct an entirely new loop using the [PeptideStubMover](https://www.rosettacommons.org/docs/latest/PeptideStubMover), the [DeclareBond mover](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/DeclareBond), and GeneralizedKIC.  The GeneralizedKIC mover was used to sample entirely new conformations of the loop.  However, there are many situations in which one may wish to perturb a loop slightly from a starting conformation without randomizing its conformation completely.  This is particularly true when using GeneralizedKIC to explore variant conformations of a starting structure, or when one has used fragment insertion to build a loop crudely and wishes to refine the starting model.  In this tutorial, we will learn how to perturb loops with GeneralizedKIC.  We will do this in the context of a Monte Carlo search of loop conformational space, in which the moves are small perturbations of the loop conformation (using GeneralizedKIC), and the acceptance criterion is the effect on a backbone-only score function.

## Exercise 2:  A RosettaScripts-Scripted Monte Carlo Search of Loop Conformational Space

### Inputs

For this exercise, we will be using one of the imperfect loop conformations from the first tutorial.  This simulates the design case, in which one may have built an initial, imperfect loop using a fragment-based or other method, and one wishes to improve one's initial model.

**The starting model for this tutorial:**
![The starting model for this tutorial](images/Exercise2_startingstruct.png)

The following `rosetta.flags` file will be used for this tutorial:

```
-nstruct 1
-beta_nov15
-in:file:s inputs/2ND2_exercise1_solution4.pdb
-in:file:fullatom
-write_all_connect_info
-parser:protocol xml/exercise2.xml
-jd2:failed_job_exception false
-mute all
-unmute protocols.simple_moves.GenericMonteCarloMover
```

We will also start with the script from the previous tutorial, and modify it to our needs.  Here it is for reference:

```xml
<ROSETTASCRIPTS>
	<SCOREFXNS>
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
	</SCOREFXNS>
	<RESIDUE_SELECTORS>
	</RESIDUE_SELECTORS>
	<TASKOPERATIONS>
	</TASKOPERATIONS>
	<FILTERS>
	</FILTERS>
	<MOVERS>	
		<PeptideStubMover name="add_loop_residues" >
			<Insert anchor_rsd="28" resname="ALA" />
			<Insert anchor_rsd="29" resname="GLY" />
			<Insert anchor_rsd="30" resname="ALA" />
			<Prepend anchor_rsd="32" resname="ALA" />
			<Prepend anchor_rsd="32" resname="ALA" />
		</PeptideStubMover>
		
		<DeclareBond name="new_bond" atom1="C" atom2="N" res1="31" res2="32" />

		<MutateResidue name="mut1" target="28" new_res="ALA" />
		<MutateResidue name="mut2" target="34" new_res="ALA" />

		<GeneralizedKIC name="genkic" selector="lowest_energy_selector" selector_scorefunction="bb_only"
			closure_attempts="5000" stop_when_n_solutions_found="5" >
			<AddResidue res_index="28" />
			<AddResidue res_index="29" />
			<AddResidue res_index="30" />
			<AddResidue res_index="31" />
			<AddResidue res_index="32" />
			<AddResidue res_index="33" />
			<AddResidue res_index="34" />
			<SetPivots res1="28" res2="31" res3="34" atom1="CA" atom2="CA" atom3="CA" />
			<AddPerturber effect="set_dihedral" >
				<AddAtoms res1="28" atom1="C" res2="29" atom2="N" />
				<AddAtoms res1="29" atom1="C" res2="30" atom2="N" />
				<AddAtoms res1="30" atom1="C" res2="31" atom2="N" />
				<AddAtoms res1="31" atom1="C" res2="32" atom2="N" />
				<AddAtoms res1="32" atom1="C" res2="33" atom2="N" />
				<AddAtoms res1="33" atom1="C" res2="34" atom2="N" />
				<AddValue value="180.0" />
			</AddPerturber>
			<CloseBond res1="31" res2="32" atom1="C" atom2="N" bondlength="1.328685" angle1="121.699997" angle2="116.199993" torsion="180.0" />
			<AddPerturber effect="randomize_backbone_by_rama_prepro" >
				<AddResidue index="28" />
				<AddResidue index="29" />
				<AddResidue index="30" />
				<AddResidue index="31" />
				<AddResidue index="32" />
				<AddResidue index="33" />
				<AddResidue index="34" />
			</AddPerturber>
			<AddFilter type="backbone_bin" residue="28" bin_params_file="ABBA" bin="A" />
			<AddFilter type="backbone_bin" residue="34" bin_params_file="ABBA" bin="A" />
			<AddFilter type="loop_bump_check" />
			<AddFilter type="rama_prepro_check" residue="28" rama_cutoff_energy="0.5" />
			<AddFilter type="rama_prepro_check" residue="31" rama_cutoff_energy="0.5" />
			<AddFilter type="rama_prepro_check" residue="34" rama_cutoff_energy="0.5" />
		</GeneralizedKIC>
		
	</MOVERS>
	<APPLY_TO_POSE>
	</APPLY_TO_POSE>
	<PROTOCOLS>
		<Add mover="add_loop_residues" />
		<Add mover="new_bond" />
		<Add mover="mut1" />
		<Add mover="mut2" />
		<Add mover="genkic" />
	</PROTOCOLS>
	<OUTPUT />
</ROSETTASCRIPTS>

```

### Step 1: Deleting unnecessary parts of the starting script

Because we are starting with an input structure that already has a loop, we do not need to build one from scratch.  We can therefore delete the PeptideStubMover, DeclareBond mover, and the MutateResidue movers from the `<MOVERS>` and `<PROTOCOLS>` sections of the starting script.

### Step 2: Modifying the GeneralizedKIC mover

We now want to modify the GeneralizedKIC mover so that, rather than completely randomizing the conformation of the loop, it merely perturbs it slightly from its starting conformation.  Let's start at the top and work our way down.  The first thing to change is the GeneralizedKIC selector.  Rather than selecting the lowest-energy conformation that we find, if there are multiple KIC solutions, we want to be sure that we're not picking one that drastically alters the loop conformation (which some of the many solutions to the system of equations that KIC solves might do).  For this reason, we will use a `lowest_delta_torsion_selector`, which chooses the solution that has the smallest RMSD, in torsion space, to the input structure.  We will also set `stop_when_n_solutions_found` to `"1"`, since we have no interest in finding large numbers of solutions and discarding most of them.  Because it is easier to find KIC solutions when one starts from a viable solution, we will set `closure_attempts` to `"100"`.  The GeneralizedKIC tag should now look like this:

```xml
<GeneralizedKIC name="genkic" selector="lowest_delta_torsion_selector" selector_scorefunction="bb_only"
	closure_attempts="100" stop_when_n_solutions_found="1" >
	...
</GeneralizedKIC>
```

We can remove all of the `<AddPerturber>` and `<CloseBond>` tags, since we no longer want to perturb in this way.  Instead, we will add a `perturb_dihedral` GeneralizedKIC perturber.  This perturber takes the input dihedral value, adds a small, random value to it, and returns the slighly perturbed dihedral.  We want to perturb both phi and psi of all loop residues:

```xml
<GeneralizedKIC ...>
	...
	<AddPerturber effect="perturb_dihedral" >
		<AddAtoms res1="28" atom1="N" res2="28" atom2="CA" />
		<AddAtoms res1="29" atom1="N" res2="29" atom2="CA" />
		<AddAtoms res1="30" atom1="N" res2="30" atom2="CA" />
		<AddAtoms res1="31" atom1="N" res2="31" atom2="CA" />
		<AddAtoms res1="32" atom1="N" res2="32" atom2="CA" />
		<AddAtoms res1="33" atom1="N" res2="33" atom2="CA" />
		<AddAtoms res1="34" atom1="N" res2="34" atom2="CA" />
		<AddAtoms res1="28" atom1="CA" res2="28" atom2="C" />
		<AddAtoms res1="29" atom1="CA" res2="29" atom2="C" />
		<AddAtoms res1="30" atom1="CA" res2="30" atom2="C" />
		<AddAtoms res1="31" atom1="CA" res2="31" atom2="C" />
		<AddAtoms res1="32" atom1="CA" res2="32" atom2="C" />
		<AddAtoms res1="33" atom1="CA" res2="33" atom2="C" />
		<AddAtoms res1="34" atom1="CA" res2="34" atom2="C" />
		<AddValue value="10.0" />
	</AddPerturber>
</GeneralizedKIC>
```

The `<AddValue>` tag, above, in this case sets the bredth of the Gaussian for the randomly-chosen value added to each dihedral value; this may roughly be thought of as the maximum size of the perturbation (though it is more accurately its standard deviation).

Next, we should review the filters that we have in place.  We can delete all of these except the `loop_bump_check` filter.

The GeneralizedKIC mover setup should now look like this:

```xml
<GeneralizedKIC name="genkic" selector="lowest_delta_torsion_selector" selector_scorefunction="bb_only"
	closure_attempts="100" stop_when_n_solutions_found="1" >
	<AddResidue res_index="28" />
	<AddResidue res_index="29" />
	<AddResidue res_index="30" />
	<AddResidue res_index="31" />
	<AddResidue res_index="32" />
	<AddResidue res_index="33" />
	<AddResidue res_index="34" />
	<SetPivots res1="28" res2="31" res3="34" atom1="CA" atom2="CA" atom3="CA" />
	<AddPerturber effect="perturb_dihedral" >
		<AddAtoms res1="28" atom1="N" res2="28" atom2="CA" />
		<AddAtoms res1="29" atom1="N" res2="29" atom2="CA" />
		<AddAtoms res1="30" atom1="N" res2="30" atom2="CA" />
		<AddAtoms res1="31" atom1="N" res2="31" atom2="CA" />
		<AddAtoms res1="32" atom1="N" res2="32" atom2="CA" />
		<AddAtoms res1="33" atom1="N" res2="33" atom2="CA" />
		<AddAtoms res1="34" atom1="N" res2="34" atom2="CA" />
		<AddAtoms res1="28" atom1="CA" res2="28" atom2="C" />
		<AddAtoms res1="29" atom1="CA" res2="29" atom2="C" />
		<AddAtoms res1="30" atom1="CA" res2="30" atom2="C" />
		<AddAtoms res1="31" atom1="CA" res2="31" atom2="C" />
		<AddAtoms res1="32" atom1="CA" res2="32" atom2="C" />
		<AddAtoms res1="33" atom1="CA" res2="33" atom2="C" />
		<AddAtoms res1="34" atom1="CA" res2="34" atom2="C" />
		<AddValue value="10.0" />
	</AddPerturber>
	<AddFilter type="loop_bump_check" />
</GeneralizedKIC>

```



## Conclusion

**TODO**

## Further Reading

Bhardwaj G, Mulligan VK, Bahl CD, Gilmore JM, Harvey PJ, Cheneval O, Buchko GW, Pulavarti SV, Kaas Q, Eletsky A, Huang PS, Johnsen WA, Greisen PJ, Rocklin GJ, Song Y, Linsky TW, Watkins A, Rettie SA, Xu X, Carter LP, Bonneau R, Olson JM, Coutsias E, Correnti CE, Szyperski T, Craik DJ, Baker D.  (2016).  Accurate de novo design of hyperstable constrained peptides.  _Nature_ 538(7625):329-335.

Mandell DJ, Coutsias EA, Kortemme T. (2009).  Sub-angstrom accuracy in protein loop reconstruction by robotics-inspired conformational sampling.  _Nat. Methods_ 6(8):551-2.

Coutsias EA, Seok C, Jacobson MP, Dill KA.  (2004).  A kinematic view of loop closure.  _J. Comput. Chem._ 25(4):510-28.

[GeneralizedKIC documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/composite_protocols/generalized_kic/GeneralizedKIC)

