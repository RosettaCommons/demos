#Relax
##The Relax protocol
Relax is the main protocol for relaxing a structure in Rosetta; that is, it samples conformations of a given structure close to it in 3d space to find the lowest-scoring variant, running both the [packer] and the [minimizer]. This is usually done to enable an apples-to-apples comparison between disparate structures, including crystal structures and the output of Rosetta's sampling protocols, by first minimizing them in local space according to the same score function. It is therefore advisable to run relax on any structures you intend to compare to each other.

To demonstrate this, run 

	$> ../../../main/source/bin/score_jd2.default.linuxclangrelease -s 1ubq.pdb @crystal_score_flags

and note the score. Now, run

	$> ../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -out:file:name relaxed.pdb @general_relax_flags

and note the dramatic difference in score compared to the relatively minor difference in structure; this may be 300-500 REU depending on the success of the relaxation. In the first case, 1ubq.pdb was pulled directly from the Protein Databank and optimized according to some crystallographic energy function. In the second, Rosetta has moved the protein until it is optimal according to its internal score function.  
To further explore this, run

	> ../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -relax:cycles N -nstruct 10

for N between 1 and 10 cycles, plotting the relationship between cycle number and average score. You should see diminishing returns, particularly for very large N, as well as increasing divergence from the starting structure. Production runs generally include between 5 and 15 cycles of Relax; 5 is most often sufficient.
## Modifying the scope of Relax
###Restricting the conformations it can sample
By default, Relax is permitted to select new side chain rotamers, move the protein backbone, and move protein subunits relative to each other; while this allows the protocol to find a more optimal solution, it can be useful to restrict Relax from modifying a structure in ways that run counter to biological data. Relax may be provided with a [MoveMap] by use of the option
	-in:file:movemap
In lieu of a specified MoveMap, the options
	-relax:chi_move false
	-relax:bb_move false
	-relax:jump_move false
will disable side chain, backbone, and interdomain motion, respectively. It can be useful, for example. to prevent motion between a designed protein and its native binding partner.
To demonstrate this, run

	$> ../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -relax:bb_move false @general_relax_flags

and align it to the original 1ubq.pdb. You should see a close alignment between the backbones before and after the run -- and a correspondingly higher final energy.
To further explore this, run

	$> ../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb @general_relax_flags -in:file:movemap lever_arm_movemap

and observe that everything c-terminal to the region allowed to move has also moved. This is endemic to [movemaps] using [relative coordinates], and is called the [lever-arm effect]; while some protocols in Rosetta are written with this in mind, Relax allows lever-arm effects if not specifically prohibited from doing so within its MoveMap.

It can also be useful to disfavor dramatic movements in Relax without completely disallowing them. This may be done by adding constraints via

	-relax:constrain_relax_to_start_coords
	-relax:constrain_relax_to_native_coords -in:file:native

The former option disfavors output that is structurally dissimilar to the input; the latter similarly disfavors divergence from a provided input file. These are implemented as harmonic constraints, so a linear divergence is reflected in a quadratic increase in score. If small changes are not to be disfavored, these constraints can be modified with 

	-relax:coord_cst_width    <width> 

which replaces the normal harmonic constraints with flat constraints out to a distance of *width*, such that any changes that leave the output within *width* of the input see no change in score at all.
To demonstrate this, run

	$> ../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -relax:constrain_relax_to_start_coords @general_relax_flags 

and compare to the results of the original, unconstrained relax run.


###Restricting the sequence it can sample
By default, Relax will not change the input sequence. It can be allowed to do so in a controlled way via [resfiles] and the options

	-relax:respect_resfile -packing:resfile *resfile*

which will set its internal packer to respect the provided resfile. This only controls packing behavior, not the minimizer; it can also be used to increase [rotamer] sampling around critical residues.
To demonstrate this, run 

	$> ../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -relax:respect_resfile -packing:resfile 1ubq.resfile @general_relax_flags

and compare to 1ubq.pdb. Note that the two structures are different; these differences arise during the minimization step.

###Changing the behavior of the repulsive algorithm
By default, Relax employs a ramping repulsive energy; this increases the structural space the algorithm is permitted to explore by gradually decreasing the weight of the repulsive term as the run continues. This behavior may be disabled by 
	-relax:ramp_constraints false
if the particular positions of the atoms in a given structure are believed to be accurate.

##FastRelax

There exists an updated version of relax called FastRelax that is capable of operating via script. The construction of these scripts is covered [here](https://www.rosettacommons.org/docs/wiki/application_documentation/structure_prediction/relax)

