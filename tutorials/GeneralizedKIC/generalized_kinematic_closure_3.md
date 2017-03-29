# Generalized Kinematic Closure Tutorial 3:
# Using pre-selection movers within GeneralizedKIC
======================================

KEYWORDS: LOOPS SCRIPTING_INTERFACES

Tutorial by Vikram K. Mulligan (vmullig@uw.edu).  Created on 28 March 2017 for the Baker lab Rosetta Tutorial Series.

[[_TOC_]]

## Goals

At the end of this tutorial, you will understand:

- What a GeneralizedKIC pre-selection mover is
- How to apply an arbitrary protocol to every solution found by GeneralizedKIC before choosing a solution
- How to use GeneralizedKIC with flexible-backbone minimization
- How to use GeneralizedKIC with design

It is highly recommended that you complete the [first](generalized_kinematic_closure_1.md) and [second](generalized_kinematic_closure_2.md) tutorials before proceeding.

## Using pre-selection movers within GeneralizedKIC

In the first two tutorials, we build, randomized, and perturbed loops, but did little with the many solutions that GeneralizedKIC might return.  A major limitation of the way in which we have been using GeneralizedKIC so far is the fact that our selection of a "best" closure solution has been limited to criteria that only consider backbone geometry.  In this tutorial, we will learn how to carry out arbitrary protocols on every closure solution _prior_ to picking a "best" solution, so that the selection may be based on less arbitrary criteria.

## Exercise 3: Designing each GeneralizedKIC solution before selecting the "best" solution

### Inputs

For this exercise, we will be using the same starting point as in the first exercise.  We will also be starting with the same script, which we will modify to allow us to design and relax every solution before picking one to return.

**The input structure, an edited version of PDB structure 2ND2 (`2ND2_state1_glyonly_loop_removed.pdb`):**
![The input structure](images/Example1_input_structure.png)

The following `rosetta.flags` file will be used for this tutorial:

```
-nstruct 10
-beta_nov15
-in:file:s inputs/2ND2_state1_glyonly_loop_removed.pdb
-in:file:fullatom
-write_all_connect_info
-parser:protocol xml/exercise3.xml
-jd2:failed_job_exception false
-mute protocols.generalized_kinematic_closure.filter.GeneralizedKICfilter core.chemical.AtomICoor core.conformation.Residue core.kinematics.AtomTree
```

Here is the script from the first tutorial for reference.  We will modify it to our needs:

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

Since we will be carrying out design in this tutorial, we will need a resfile defining the amino acids with which we can design.  We will not discuss the details of design here; instead, please refer to the [[Protein Design Tutorial|protein_design_tutorial]].  Here is the resfile that we will use, defining a limited palette of amino acids:

```
PIKAA AGPILYVEKR
start
```

Finally, we need to define a suitable FoldTree to allow loop relaxation without distorting the structure as a whole.  Again, we will not discuss the details of FoldTrees here; instead, refer to the [[FoldTree tutorial|fold_tree]].

### Step 1:

**TODO**

## Running the example script

The above script is provided in the `demos/tutorials/GeneralizedKIC/exercise3/xml/` directory.  To run this, navigate to the `demos/tutorials/GeneralizedKIC` directory and type the following:

```bash
$> cd exercise3
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @inputs/rosetta.flags
$> cd ..
```

In the above, `$ROSETTA3` is the path to your Rosetta directory.  You may need to replace `linuxgccrelease` for your operating system and compilation (_e.g._ `macosclangrelease` on a Mac).

## Expected output

When tested with Rosetta 3.8 SHA 3cad483ccac973741499159e12989a7143bf79de (nightly build from Tuesday, March 28th, 2017), the script produced various loop conformations and sequences, some of which packed well or had favourable interactions between the loop and the rest of the protein.  Because of the stringent filtering used, some replicates failed to return a result.  Note that this tutorial has not been optimized to yield _good_ designs.

**An example designed structure**
![An example designed structure](images/Tutorial3_example_output.png)

## Conclusion

In this tutorial, we have covered loop conformational perturbation and Monte Carlo searches using GeneralizedKIC.  The reader should experiment with settings to learn how different perturbation magnitudes and Monte Carlo temperatures affect the trajectory.

## Further Reading

Bhardwaj G, Mulligan VK, Bahl CD, Gilmore JM, Harvey PJ, Cheneval O, Buchko GW, Pulavarti SV, Kaas Q, Eletsky A, Huang PS, Johnsen WA, Greisen PJ, Rocklin GJ, Song Y, Linsky TW, Watkins A, Rettie SA, Xu X, Carter LP, Bonneau R, Olson JM, Coutsias E, Correnti CE, Szyperski T, Craik DJ, Baker D.  (2016).  Accurate de novo design of hyperstable constrained peptides.  _Nature_ 538(7625):329-335.

Mandell DJ, Coutsias EA, Kortemme T. (2009).  Sub-angstrom accuracy in protein loop reconstruction by robotics-inspired conformational sampling.  _Nat. Methods_ 6(8):551-2.

Coutsias EA, Seok C, Jacobson MP, Dill KA.  (2004).  A kinematic view of loop closure.  _J. Comput. Chem._ 25(4):510-28.

[GeneralizedKIC documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/composite_protocols/generalized_kic/GeneralizedKIC)

[GeneralizedKIC Tutorial 1](generalized_kinematic_closure_1.md)

[GeneralizedKIC Tutorial 2](generalized_kinematic_closure_2.md)

[GeneralizedKIC Tutorial 4](generalized_kinematic_closure_4.md)
