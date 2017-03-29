# Generalized Kinematic Closure Tutorial 4:
# Using GeneralizedKIC to close through disulfide bonds
======================================

KEYWORDS: LOOPS SCRIPTING_INTERFACES

Tutorial by Vikram K. Mulligan (vmullig@uw.edu).  Created on 28 March 2017 for the Baker lab Rosetta Tutorial Series.

[[_TOC_]]

## Goals

At the end of this tutorial, you will understand:

- How broadly GeneralizedKIC defines the concept of "loop"
- What a "tail residue" is, and why it is important
- How to close a loop that passes through connections other than mainchain peptide bonds (in this case, through a sidechain disulfide linkage)
- How to keep segments rigid during GeneralizedKIC sampling

It is highly recommended that you complete the [first](generalized_kinematic_closure_1.md), [second](generalized_kinematic_closure_2.md), and [third](generalized_kinematic_closure_3.md) tutorials before proceeding.

## Using GeneralizedKIC to close through non-backbone connections

In all examples seen so far, we have used GeneralizedKIC to close a segment of alpha-amino acid backbone.  GeneralizedKIC is written to be completely general, however: it can close through arbitary chains of atoms, including segments of canonical or non-canonical backbone, side-chain connections, covalent cross-linkers, metal ions, _etc._.  The only requirements are that the chain have at least three atoms flanked by freely rotatable bonds (and separated by at least one atom), and that all residues in the chain are covalently linked.  In this exercise, we will close a loop that passes through a segment of backbone, an alpha helix, and a side-chain disulfide connection, closing gaps both in the backbone loop and in the side-chain sulphur-sulphur bond.

## Exercise 4: Using GeneralizedKIC to close through disulfide bonds

### Inputs

For this exercise, we will be using the same starting point as in the first exercise.  We will also be starting with the same script.  However, we will modify it so that it not only builds and closes a loop between helix 2 and helix 3, but also forms a disulfide between helix 3 and the helix 1-2 loop.  We will define our loop to be closed as one starting at residue 28, continuing through the alpha helix to residue 42 (which we will mutate to cysteine), and then passing through the cysteine sidechain to the CB atom of another cysteine added at position 19.  As such, the entire third alpha helix will be able to shift its position during sampling; the fixed points are now residue 27 and residue 19.

**The input structure, an edited version of PDB structure 2ND2 (`2ND2_state1_glyonly_loop_removed.pdb`):**
![The input structure](images/Example1_input_structure.png)

The following `rosetta.flags` file will be used for this tutorial:

```
-nstruct 10
-beta_nov15
-in:file:s inputs/2ND2_state1_glyonly_loop_removed.pdb
-in:file:fullatom
-write_all_connect_info
-parser:protocol xml/exercise4.xml
-jd2:failed_job_exception false
-mute protocols.generalized_kinematic_closure.filter.GeneralizedKICfilter core.chemical.AtomICoor core.conformation.Residue
```

Here is the script from the first tutorial for reference.  We will modify it to fit our needs:

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

### Step 1:

**TODO**

## Running the example script

The above script is provided in the `demos/tutorials/GeneralizedKIC/exercise4/xml/` directory.  To run this, navigate to the `demos/tutorials/GeneralizedKIC` directory and type the following:

```bash
$> cd exercise4
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @inputs/rosetta.flags
$> cd ..
```

In the above, `$ROSETTA3` is the path to your Rosetta directory.  You may need to replace `linuxgccrelease` for your operating system and compilation (_e.g._ `macosclangrelease` on a Mac).

## Expected output

When tested with Rosetta 3.8 SHA 3cad483ccac973741499159e12989a7143bf79de (nightly build from Tuesday, March 28th, 2017), the script rapidly produced many possible placements for helix 3 relative to helices 1 and 2, with ideal geometry both for the loop between helices 2 and 3 and for the disulfide bond between residues 42 and 19.

**Some of the sampled conformations (green -- helices 1 and 2; cyan -- helix 3)**

![Some of the sampled conformations](images/Tutorial4_loop_closure.gif)

## Conclusion

**TODO**

## Further Reading

Bhardwaj G, Mulligan VK, Bahl CD, Gilmore JM, Harvey PJ, Cheneval O, Buchko GW, Pulavarti SV, Kaas Q, Eletsky A, Huang PS, Johnsen WA, Greisen PJ, Rocklin GJ, Song Y, Linsky TW, Watkins A, Rettie SA, Xu X, Carter LP, Bonneau R, Olson JM, Coutsias E, Correnti CE, Szyperski T, Craik DJ, Baker D.  (2016).  Accurate de novo design of hyperstable constrained peptides.  _Nature_ 538(7625):329-335.

Mandell DJ, Coutsias EA, Kortemme T. (2009).  Sub-angstrom accuracy in protein loop reconstruction by robotics-inspired conformational sampling.  _Nat. Methods_ 6(8):551-2.

Coutsias EA, Seok C, Jacobson MP, Dill KA.  (2004).  A kinematic view of loop closure.  _J. Comput. Chem._ 25(4):510-28.

[GeneralizedKIC documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/composite_protocols/generalized_kic/GeneralizedKIC)

[GeneralizedKIC Tutorial 1](generalized_kinematic_closure_1.md)

[GeneralizedKIC Tutorial 2](generalized_kinematic_closure_2.md)

[GeneralizedKIC Tutorial 3](generalized_kinematic_closure_3.md)
