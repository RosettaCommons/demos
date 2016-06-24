#Relax

Written by Frank Teets
Last Modified Jun 21 2016

KEYWORDS: STRUCTURE_PREDICTION GENERAL

[[_TOC_]]

##The Relax protocol

Relax is the main protocol for relaxing a structure in Rosetta; that is, it samples conformations of a given structure close to it in 3D space to find the lowest-scoring variant, running both the [[packer|Optimizing_Sidechains_The_Packer]] and the [[minimizer|minimization]] alternately. This is usually done to enable an apples-to-apples comparison between disparate structures, including crystal structures and the output of Rosetta's sampling protocols, by first minimizing them in local space according to the same score function. It is therefore advisable to run relax on any structures you intend to compare to each other.

##Navigating to the Demos
The demos are available at `<path_to_Rosetta_directory>/demos/tutorials/Relax_Tutorial`. All demo commands listed in this tutorial should be executed when in this directory. All the demos here use the `linuxclangrelease` binary. You may be required to change it to whatever is appropriate given your operating system and compiler.

##Demo

To demonstrate how to run the `relax` executable on an input structure, run the following command in your terminal: 

	$>../../../main/source/bin/score_jd2.default.linuxclangrelease -s 1ubq.pdb -out:suffix _crystal @crystal_score_flags

and note the score in the `score_crystal.sc` file. Now, run

	$>../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -out:suffix _relaxed @general_relax_flags

where the contents of the `general_relax_flags` file are expanded below:

```
-nstruct 2
-relax:default_repeats 5
-out:path:pdb ./tutorial_output
-out:path:score ./expected_output
```

These flags instruct Rosetta to run two separate relax trajectories on the input structure (`-nstruct 2`), to perform the relaxation algorithm with five cycles of sidechain repack and minimization (`-relax:range:cycles 5`), and to output the relaxed structures in a directory that is in this folder named `tutorial_output/` (`-out:path:pdb`) and the scorefile in a directory that is in this folder named `expected_output/` (`-out:path:score`).

The previous command should have generate two structures (as per the `-nstruct 2` option) in the `tutorial_output/` folder, namely `1ubq_relaxed_0001.pdb` and `1ubq_relaxed_0002.pdb`. It will also create a `score_relaxed.sc` file in the `expected_output/` directory. Open this score file and note the dramatic difference in score compared to the relatively minor difference in structure; this may be 300-500 REU depending on the success of the relaxation. In the first case, 1ubq.pdb was pulled directly from the [Protein Data Bank](http://www.rcsb.org/pdb/home/home.do) and optimized according to some crystallographic energy function. In the second, Rosetta has moved the protein until it is optimal according to its internal score function.  

To further explore this, change "N" iteratively in the following command to the integers between 1 and 10

	>../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -out:suffix _N_relax_cycles -relax:default_repeats N -nstruct 10 -out:path:pdb ./tutorial_output -out:path:score ./expected_output

which will run the relax protocol with 1 to 10 cycles. Then, plot the relationship between cycle number and average score. You should see diminishing returns, particularly for very large N, as well as increasing divergence from the starting structure. Production runs generally include between 5 and 15 cycles of Relax; 5 is most often sufficient. Five repeats is the default, so the `-relax:default_repeats` option can normally be omitted.

## Modifying the scope of Relax

### Restricting the conformations that Relax can sample

By default, Relax is permitted to select new side chain rotamers, move the protein backbone, and move protein subunits relative to each other. While this allows the protocol to find a more optimal solution, it can be useful to restrict Relax from modifying a structure in ways that run counter to biological data. Relax may be provided with a [MoveMap] by use of the option

```
-in:file:movemap
```

In lieu of a specified MoveMap, the options

```
-relax:chi_move false
-relax:bb_move false
-relax:jump_move false
```

will disable side chain, backbone, and interdomain motion, respectively. It can be useful, for example to prevent motion between a designed protein and its native binding partner.

To explore how a movemap influences the `relax` application, run

	$>../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb @general_relax_flags -out:suffix _lever_arm -in:file:movemap lever_arm_movemap

and observe that everything C-terminal to the region allowed to move has also moved. This is endemic to [[movemaps|minimization]] using [[internal coordinates|Core_Concepts]], and is called the [[lever-arm effect|Core_Concepts]]; while some protocols in Rosetta are written with this in mind, Relax allows lever-arm effects if not specifically prohibited from doing so within its MoveMap.

To demonstrate how the `-relax:chi_move`, `-relax:bb_move`, and `-relax:jump_move` flags affect the final structure, run

	$>../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -relax:bb_move false -out:suffix _no_sidechain_motion @general_relax_flags

and align one or both of the two output structures in `tutorial_output/` directory to the original 1ubq.pdb. You should see a close alignment between the backbones before and after the run -- and a correspondingly higher final energy.

It can also be useful to disfavor dramatic movements in Relax without completely disallowing them. This may be done by adding constraints via

	-relax:constrain_relax_to_start_coords
	-relax:constrain_relax_to_native_coords -in:file:native <filename.pdb>

The former option disfavors output that is structurally dissimilar to the input; the latter similarly disfavors divergence from a provided input file. As the name implies, this is particularly useful for ensuring fidelity to some kind of native structure. These are implemented as harmonic constraints, so a linear divergence is reflected in a quadratic increase in score. If small changes are not to be disfavored, these constraints can be modified with 

	-relax:coord_cst_width    <width> 

which replaces the normal harmonic constraints with flat harmonic constraints out to a distance of *width*, such that any changes that leave the output within *width* of the input see no change in score at all.

To demonstrate this, run

	$>../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -out:suffix _constrained_relax -relax:constrain_relax_to_start_coords @general_relax_flags 

and compare to the results of the original, unconstrained relax run.

###Changing the behavior of constraints

By default, Relax gradually decreases the weight of any constraints given as the run progresses in order to explore more space in which the constraints hold while also optimizing the final structure. This behavior may be deactivated with

	-relax:ramp_constraints false

if the constraints must be absolutely maintained. 

This option is particularly useful when using `-relax:constrain_relax_to_start_coords` and similar options, as it keeps the output coordinates closer to the specified coordinates than a ramped constraint run does.

###Restricting the sequence it can sample

By default, Relax will not change the input sequence. It can be allowed to do so in a controlled way via [[resfiles|Optimizing_Sidechains_The_Packer]] and the options

	-relax:respect_resfile -packing:resfile *resfile*

which will set its internal packer to respect the provided resfile. This only controls packing behavior, not the minimizer; it can also be used to increase rotamer sampling around critical residues, but will not by itself ensure that particular rotamers are preserved.

To demonstrate this, run 

	$> ../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -out:suffix _relaxed_with_resfile -relax:respect_resfile -packing:resfile 1ubq.resfile @general_relax_flags

and compare to 1ubq.pdb. Note that the two structures are different; these differences arise during the minimization step.

Relax may also be allowed to change the sequence globally through the use of the -disable_design false option. 

## More information

For more information about Relax, including the ability to precisely control relax via script, see the [relax application documentation](https://www.rosettacommons.org/docs/latest/application_documentation/structure_prediction/relax).

