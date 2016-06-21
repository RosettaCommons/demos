[[_TOC_]]
Demos
=====

These demos are designed to guide users through sample procedures in computational modeling from the point of view of solving a specific problem.
Each one includes:

1. An introduction to the task at hand.
2. Detailed step-by-step instructions on how to run the demo.
3. All the input data needed to run the demo.
4. Scripts containing the exact commands needed to run the demo.

Downloading the demos
=====================

All the demos documented here are included in the `demos` directory at the root of the weekly releases of Rosetta.
Weekly releases are available to academic users free of charge.
If you are a [[Rosetta Developer]] (i.e. you have signed the RosettaCommons developer's agreement and have access to the Rosetta source code repositories on GitHub), you can download all the demos by checking out the `demos.git` repository:

    git clone git@github.com:RosettaCommons/demos.git

Demos
=====

How-Tos
-------

How-to demos demonstrate the best way to solve a particular problem using the Rosetta framework.
Each demo defines a particular problem, and as new and improved ways to solve those problems are developed, the demos are updated by members of the Rosetta community.
These demos are located in the `demos/public` directory.

<<LinkDemos(public)>>

Protocol Captures
-----------------

Protocol captures are demos from published papers. 
They aren't meant to show the best way to solve problems in the current version of Rosetta, they are meant to show published solutions to problems that were addressed by members of the Rosetta community.
The purpose of these protocol captures is both to serve as a historical record and to assist those trying to reproduce past results.
These demos are located in the `demos/protocol_capture` directory.

<<LinkDemos(protocol_capture)>>

<!--- BEGIN_INTERNAL --->

Adding new demos
================

To add new demos, the first thing you need to do is become a [[Rosetta developer]].
Then you will be able to check out the `demos.git` repository:

    git clone git@github.com:RosettaCommons/demos.git

The demos are organized into three directories:

* `demos/public`:  
  For demos that are meant to show the best way to solve a particular problem.
These demos may be updated by the community as new ways to solve these problems are developed.

* `demos/protocol_capture`:  
  For demos that are associated with published papers and tht demostrate the specific algorithm described in that paper.
These demos are static, may only work with previous old of rosetta, and meant to serve more as a historical records.

* `demos/pilot`:  
  For demos that aren't meant to be included in the weekly releases yet.  I can't think of any good reasons for putting new demos here.

Each demo should go in its own directory within one of these three directories.
So to add a new demo, the first step is to create a descriptively named directory for it in the proper location.
For example, this is how you'd make a public demo called `my_demo`:

    cd demos/public
    mkdir my_demo

The only file that absolutely has to be in this directory is `README.md`.
A link to this file will automatically be added to the list of demos on this page.
Your readme should contain your name and email address, any relevant citations, a description of the problem your demo solves, links to all the scripts and input data your demo uses, and step-by-step instructions on how to run your demo.

Scripts and input files should of course be included in the directory as well.
How you organize these is up to you.
If you don't have many files, maybe just put everything in one directory.
If you have lots of files, maybe organize them into subdirectories.
Whatever makes the most sense for your demo.

Once you have finished writing your demo and have made sure that it runs properly, commit your changes and push them like usual: `git commit && git push`.
At this point, you can now edit your readme using the online Gollum interface available at the [[internal documentation site]].
However, only your readme can be maintained in this way.
You have to use git if you want to make any changes to your demo scripts or input data.

A few days after you push you demo to the `demos` repository, your demo will 
become available from this website.  A link to it will automatically be added 
to this page.  However, it will be hard to find (especially for new users) with 
just that.  To make your demo easier to find, spend a few minutes browsing the 
[[documentation wiki|https://www.rosettacommons.org/docs/wiki]] and adding 
links to any relevant pages you find.  The application section in particular 
would benefit from having lots of links to demos.  Note that (for technical 
reasons — Gollum gets really slow when there are too many pages in the wiki) 
the demos wiki is actually a whole different website than the documentation 
wiki.  So you have to use external links to link between the two wikis.

Under Development
-----------------

If you want to prevent a demo from being published to the static wiki until you've finished developing it or written a paper on it or something, put it in `demos/pilot`.
It will be visible here, but not on the public site.

<<LinkDemos(pilot)>>

<!--- END_INTERNAL --->
Demos
=====

These demos are designed to guide users through sample procedures in computational modeling from the point of view of solving a specific problem.
Each one includes:

1. An introduction to the task at hand.
2. Detailed step-by-step instructions on how to run the demo.
3. All the input data needed to run the demo.
4. Scripts containing the exact commands needed to run the demo.

Downloading the demos
=====================

All the demos documented here are included in the `demos` directory at the root of the weekly releases of Rosetta.
Weekly releases are available to academic users free of charge.
If you are a [[Rosetta Developer]] (i.e. you have signed the RosettaCommons developer's agreement and have access to the Rosetta source code repositories on GitHub), you can download all the demos by checking out the `demos.git` repository:

    git clone git@github.com:RosettaCommons/demos.git

Demos
=====

How-Tos
-------

How-to demos demonstrate the best way to solve a particular problem using the Rosetta framework.
Each demo defines a particular problem, and as new and improved ways to solve those problems are developed, the demos are updated by members of the Rosetta community.
These demos are located in the `demos/public` directory.

<<LinkDemos(public)>>

Protocol Captures
-----------------

Protocol captures are demos from published papers. 
They aren't meant to show the best way to solve problems in the current version of Rosetta, they are meant to show published solutions to problems that were addressed by members of the Rosetta community.
The purpose of these protocol captures is both to serve as a historical record and to assist those trying to reproduce past results.
These demos are located in the `demos/protocol_capture` directory.

<<LinkDemos(protocol_capture)>>

<!--- BEGIN_INTERNAL --->

Adding new demos
================

To add new demos, the first thing you need to do is become a [[Rosetta developer]].
Then you will be able to check out the `demos.git` repository:

    git clone git@github.com:RosettaCommons/demos.git

The demos are organized into three directories:

* `demos/public`:  
  For demos that are meant to show the best way to solve a particular problem.
These demos may be updated by the community as new ways to solve these problems are developed.

* `demos/protocol_capture`:  
  For demos that are associated with published papers and tht demostrate the specific algorithm described in that paper.
These demos are static, may only work with previous old of rosetta, and meant to serve more as a historical records.

* `demos/pilot`:  
  For demos that aren't meant to be included in the weekly releases yet.  I can't think of any good reasons for putting new demos here.

Each demo should go in its own directory within one of these three directories.
So to add a new demo, the first step is to create a descriptively named directory for it in the proper location.
For example, this is how you'd make a public demo called `my_demo`:

    cd demos/public
    mkdir my_demo

The only file that absolutely has to be in this directory is `README.md`.
A link to this file will automatically be added to the list of demos on this page.
Your readme should contain your name and email address, any relevant citations, a description of the problem your demo solves, links to all the scripts and input data your demo uses, and step-by-step instructions on how to run your demo.

Scripts and input files should of course be included in the directory as well.
How you organize these is up to you.
If you don't have many files, maybe just put everything in one directory.
If you have lots of files, maybe organize them into subdirectories.
Whatever makes the most sense for your demo.

Once you have finished writing your demo and have made sure that it runs properly, commit your changes and push them like usual: `git commit && git push`.
At this point, you can now edit your readme using the online Gollum interface available at the [[internal documentation site]].
However, only your readme can be maintained in this way.
You have to use git if you want to make any changes to your demo scripts or input data.

A few days after you push you demo to the `demos` repository, your demo will 
become available from this website.  A link to it will automatically be added 
to this page.  However, it will be hard to find (especially for new users) with 
just that.  To make your demo easier to find, spend a few minutes browsing the 
[[documentation wiki|https://www.rosettacommons.org/docs/wiki]] and adding 
links to any relevant pages you find.  The application section in particular 
would benefit from having lots of links to demos.  Note that (for technical 
reasons — Gollum gets really slow when there are too many pages in the wiki) 
the demos wiki is actually a whole different website than the documentation 
wiki.  So you have to use external links to link between the two wikis.

Under Development
-----------------

If you want to prevent a demo from being published to the static wiki until you've finished developing it or written a paper on it or something, put it in `demos/pilot`.
It will be visible here, but not on the public site.

<<LinkDemos(pilot)>>

<!--- END_INTERNAL --->
#How To Read These Tutorials
These tutorials were written such that a completely new user should be able to complete them in numeric order through Tutorial 10 and thereby gain an understanding of the basic mechanics of Rosetta. Prior to doing any of the later tutorials, run through Tutorials 1-10.
##Before Running Any of the Other Tutorials
Complete Tutorial 1 to install and compile rosetta; verify that the Rosetta/main/bin directory contains executables appropriate to your installation.
##Do The Following for Each Tutorial
In order for the hands-on portions of these tutorials to function correctly, you must make your current working directory that of the tutorial you want to run; IE for this tutorial, your current working directory must be (IE "you must be in") Rosetta/demos/tutorials/tutorial_0.
Therefore, cd into the directory for the tutorial you want to run; see Tutorial 2 for a brief description of how to use cd.
If you have previously run a tutorial and wish to re-run it, make sure you delete the contents of the tutorial_output folder; by default, Rosetta will abandon a job rather than overwrite an existing output file, so the output directory must be clean for the tutorials to properly execute.


# Working With Rosetta
##Required Software
Rosetta relies upon a number of additional programs to function. In order to productively use Rosetta, you must have:

1. A text editor capable of outputting plain text files; most people use [emacs](https://www.gnu.org/software/emacs/manual/html_node/emacs/index.html),[vi](https://www.washington.edu/computing/unix/vi.html),[nano](https://www.nano-editor.org/dist/v2.0/nano.html), or some more graphical way of interfacing with one of them.
2. A molecular viewing tool. Rosetta does not include a way to actually visualize the output files it creates; a tool such as [PyMOL](https://www.pymol.org/) is necessary to examine the output PDBs. Rosetta does include a PyMOLObserver for directly viewing its output in Pymol; instructions for setting up PyMOLObserver are found [here].
3. A terminal. Unix and Apple users have suitable terminals installed by default. Windows users attempting to run Rosetta remotely on a Mac or Unix machine will require a tool like PuTTY to provide the necessary interface.

It may also help to familiarize yourself with the basics of [command line interfaces.](https://bash.cyberciti.biz/guide/Main_Page) Particularly useful commands include:

	ls *foo* 	displays a list of files in the directory *foo*
	cd *foo* 	moves the current working directory to the directory *foo*
	mv *foo* *bar* 	moves all of *foo* into *bar* if *bar* is a directory or into a file called *bar* if not
	cp *foo* *bar* 	moves all of *foo* into *bar* if *bar* is a directory or into a file called *bar* if not
	ln -s *foo* 	*path_to_bar* establishes a link called a *symbolic link* or *symlink* from *foo* to *bar* 
	top 		lists the currently running processes, together with their process IDs. This is useful for:
	kill *foo_pid* stops the process with process ID *foo_pid*
	nohup		when prepended to a command, runs that command in the background, preventing it from writing output to the terminal
	

#The Three Main Interfaces to Rosetta
Rosetta is not a single program; it is instead built with a large number of executable files built to be run either individually, sequentially, or in concert. These files may be run directly from the command line, through an XML interface called RosettaScripts, or through Python via Pyrosetta.
##Running the Rosetta Executables directly
The executables that comprise Rosetta can be run directly by calling them through the terminal directly. These commands are constructed by first indicating the path to the executable you want to run, followed by the flags, options, and input files that dictate how you want the program to execute. 
As an example, consider the command 

	$>../../../main/source/bin/fixbb.default.linuxgccrelease -s 1ubq.pdb -out:suffix _my_first_rosetta_output @general_fixbb_flags

The first part of this command is the path to and name of the executable we want to run; here, we want to run fixbb, located in main/source/bin, which will attempt to design residues onto a peptide backbone provided to it. On your platform, the executable may not be named fixbb.linuxgccrelease; to find the version of the executable that will run, navigate to main/source/bin, type "ls fixbb", and then hit Tab to display the versions of fixbb that have been compiled on your system. This works for any other executable as well.
The second part, -s 1ubq.pdb, is an option, or *flag*, that modifies the execution of the program. In this case, the "-s" option indicates that we want to run fixbb on a single structure; the following "1ubq.pdb" is simply the filename of that structure. The counterpart to the -s option is -l, which takes as its argument a newline-delimited list of input files and is broadly conceptually equivalent to running multiple -s commands in sequence. Other common options include:

	-nstruct: This option dictates the number of output structures Rosetta is to generate; this is important for [Monte Carlo simulations].
	-database: This option should not be required for a successful Rosetta install, but indicates the path to the Rosetta database (which is now detected automatically.)
	-out:file:path [...] This option indicates the relative location to which Rosetta is to write any output files.
	-out:suffix This option adds a specified string to the end of all the output PDB files.
	-ignore_unrecognized_residue This option should not be required for a successful Rosetta install, but indicates that residues that Rosetta cannot process should be left out of the execution.

In order to demonstrate how Rosetta exectuables may be run and their execution controlled via flags, create a folder in your main Rosetta directory named *tutorials* or something similar, and there put the files 1ubq.pdb and 1qys.pdb. 

Now run 

	> ../../../main/source/bin/fixbb.default.linuxgccrelease -s 1ubq.pdb -out:suffix _my_second_rosetta_output @general_fixbb_flags

from within the tutorials directory you just created and observe that the program generates an output pdb file(1ubq_0001.pdb), a score file (score.sc) and a log file (output).
Try running the following commands:

	> ../../../main/source/bin/fixbb.default.linuxgccrelease -s 1ubq.pdb -nstruct 10 -out:suffix _my_third_rosetta_output @general_fixbb_flags
	> ../../../main/source/bin/fixbb.default.linuxgccrelease -s 1ubq.pdb -out:suffix _my_fourth_rosetta_output -ignore_unrecognized_residue @general_fixbb_flags

The first command should produce ten output files; the second should run similarly, albeit not identically, to the command provided at the start of this tutorial.
Also, try making a list of PDBs named *pdblist* including 1ubq and 1qys in the following format:
	
	1ubq.pdb
	1qys.pdb
	
and running

	$> ../../../main/source/bin/fixbb.default.linuxgccrelease -l pdblist

##Running Rosetta via RosettaScripts
Running RosettaScripts is similar to running an executable directly from the command line: execute rosetta_scripts and pass in the relevant options. The key difference is that the executable run is always rosetta_scripts, which requires the option "parser:protocol" indicating the XML script that it is to run. RosettaScripts is covered in more detail [here]. Try creating a file named "fixbb_script.xml" containing the following:
	<ROSETTASCRIPTS>
	script goes here
	</ROSETTASCRIPTS>
saving it, and running it with 

	>../../../main/source/bin/rosetta_scripts.default.linuxgccrelease -parser:protocol fixbb_script.xml'. You should see output similar to the command run at the start of this tutorial.

##Running Rosetta via PyRosetta
PyRosetta functions are run within Python scripts; it requires that rosetta be imported with the command
	from rosetta import /*
or similar, after which the Pyrosetta functions may be called within the script itself. PyRosetta tutorials may be called [here].

#Constructing a Rosetta Job

##Considerations of Scale
Rosetta is capable of running processes at many scales. The deterministic scoring functions may be run on a single processor in a matter of seconds; running an ab initio folding job can require a national-scale supercomputer running for thousands of CPU-hours. It may be worthwhile to review the CPU time required to run Rosetta for different tasks [here].

##Post-Processing
Because Rosetta is stochastic software, it is often necessary to perform statistical analysis on the structures it generates, called *decoys*. This may be done with a statistics package such as R, documentation for which is [here](https://cran.r-project.org/doc/manuals/r-release/R-intro.pdf).

# Where to Find Things In Rosetta
###Binaries
The Rosetta executables are all symlinked into Rosetta/source/bin and may be run from there; the actual executables are compiled into Rosetta/source/build, but should be run from source/bin to ensure robustness to updates.
###Source Code
The Rosetta source code exists in Rosetta/source/src, but does not contain content required for non-developer users of Rosetta.
###Weights files
Rosetta weights files (.wts files), used to parameterize the scorefunction weights, are in Rosetta/main/database/scoring/weights .
##Params files
Parameter files are in the corresponding subdirectories of /Rosetta/main/database/chemical/residue_type_sets .
###Tools
The tools directory, Rosetta/tools, contains a number of useful tools for manipulating Rosetta inputs and outputs.
###Scripts
RosettaScripts and useful python scripts may be found in Rosetta/main/source/scripts
#Relax
##The Relax protocol
Relax is the main protocol for relaxing a structure in Rosetta; that is, it samples conformations of a given structure close to it in 3d space to find the lowest-scoring variant, running both the [packer] and the [minimizer]. This is usually done to enable an apples-to-apples comparison between disparate structures, including crystal structures and the output of Rosetta's sampling protocols, by first minimizing them in local space according to the same score function. It is therefore advisable to run relax on any structures you intend to compare to each other.

To demonstrate this, run 

	$>../../../main/source/bin/score_jd2.default.linuxclangrelease -s 1ubq.pdb -out:suffix _crystal @crystal_score_flags

and note the score. Now, run

	$>../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -out:file:name relaxed.pdb -out:suffix _relaxed @general_relax_flags

and note the dramatic difference in score compared to the relatively minor difference in structure; this may be 300-500 REU depending on the success of the relaxation. In the first case, 1ubq.pdb was pulled directly from the Protein Databank and optimized according to some crystallographic energy function. In the second, Rosetta has moved the protein until it is optimal according to its internal score function.  
To further explore this, run

	>../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -out:suffix _N_relax_cycles -relax:cycles N -nstruct 10

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

	$>../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -relax:bb_move false -out:suffix _no_sidechain_motion @general_relax_flags

and align it to the original 1ubq.pdb. You should see a close alignment between the backbones before and after the run -- and a correspondingly higher final energy.
To further explore this, run

	$>../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb @general_relax_flags -out:suffix _lever_arm -in:file:movemap lever_arm_movemap

and observe that everything c-terminal to the region allowed to move has also moved. This is endemic to [movemaps] using [internal coordinates], and is called the [lever-arm effect]; while some protocols in Rosetta are written with this in mind, Relax allows lever-arm effects if not specifically prohibited from doing so within its MoveMap.

It can also be useful to disfavor dramatic movements in Relax without completely disallowing them. This may be done by adding constraints via

	-relax:constrain_relax_to_start_coords
	-relax:constrain_relax_to_native_coords -in:file:native

The former option disfavors output that is structurally dissimilar to the input; the latter similarly disfavors divergence from a provided input file. As the name implies, this is particularly useful for ensuring fidelity to some kind of native structure. These are implemented as harmonic constraints, so a linear divergence is reflected in a quadratic increase in score. If small changes are not to be disfavored, these constraints can be modified with 

	-relax:coord_cst_width    <width> 

which replaces the normal harmonic constraints with flat harmonic constraints out to a distance of *width*, such that any changes that leave the output within *width* of the input see no change in score at all.
To demonstrate this, run

	$>../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -out:suffix _constrained_relax -relax:constrain_relax_to_start_coords @general_relax_flags 

and compare to the results of the original, unconstrained relax run.


###Restricting the sequence it can sample
By default, Relax will not change the input sequence. It can be allowed to do so in a controlled way via [resfiles] and the options

	-relax:respect_resfile -packing:resfile *resfile*

which will set its internal packer to respect the provided resfile. This only controls packing behavior, not the minimizer; it can also be used to increase [rotamer] sampling around critical residues, but will not by itself ensure that particular rotamers are preserved.
To demonstrate this, run 

	$> ../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -out:suffix _relaxed_with_resfile -relax:respect_resfile -packing:resfile 1ubq.resfile @general_relax_flags

and compare to 1ubq.pdb. Note that the two structures are different; these differences arise during the minimization step.

###Changing the behavior of constraints
By default, Relax gradually decreases the weight of any constraints given as the run progresses in order to explore more space in which the constraints hold while also optimizing the final structure. This behavior may be deactivate with
	-relax:ramp_constraints false
if the constraints must be absolutely maintained. 

##FastRelax

There exists an updated version of relax called FastRelax that is capable of operating via script. The construction of these scripts is covered [here](https://www.rosettacommons.org/docs/wiki/application_documentation/structure_prediction/relax)

#Constraints
Many of the biological problems users wish to solve with Rosetta involve some biological or functional considerations that may not be reflected within a PDB file or evaluated by normal score functions. Constraints are a general way of scoring how well a structure adheres to these additional considerations; for example, one might wish to relax a structure with constraints in place to ensure that suspected disulfides are maintained.

Constraints are written like so: a geometrical function is written that, in a perfect sequence, will return some value N. This value, together with the output of the function in the structure under consideration, are compared by some function, and the output of that function multiplied by the some scalar weight and added to the score.

For example, a simple constraint might measure the distance between two atoms, subtract the ideal distance, and subtract the difference from the score. That constraint might look like this:
	
	AtomPair CA 20 CA 6 LINEAR_PENALTY 9.0 0 0 1.0 1.0 

The constraint begins with the definition of some geometrical ideal.  In this case, since we want to constrain two atoms to be a specific distance away from each other, we want an AtomPair constraint. The two atoms are defined as the alpha-carbons of residues 6 and 20 by the next four fields, while the seventh field indicates that in ideality they would be 9.0 Angstroms apart. This, then, is the geometrical function: if this structure were perfect, the alpha-carbons of residues 6 and 20 would be 9A apart.

The remaining fields dictate how the score of the structure should be penalized according to that ideal. LINEAR_PENALTY indicates that the actual result is to be compared linearly to the ideal result; the following fields indicate some parameters of that function. LINEAR_PENALTY includes the ability to add a flat zone in which the function returns a constant value for a number of distances. Here, that value is 0, and the width of the flat zone is also 0, so any deviation from ideality will be reflected in the score, but perfect ideality will not impact the score at all. LINEAR_PENALTY also allows for varying the slope of the line. 

To demonstrate this, begin by relaxing ubiquitin via the following command:

	$> ../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb -out:suffix _unconstrained @general_relax_flags

as in Tutorial 4 and compare the output to that of 

	$> ../../../main/source/bin/relax.default.linuxclangrelease -s 1ubq.pdb @general_relax_flags -out:suffix _unreasonably_constrained @unreasonable_constraint_flags

while varying the weight of the constraint (the -constraint:cst_fa_weight in unreasonable_constraints.) High (>1000) weights should produce demonstrably aphysical structures, unless you set the constraint ideal to something close to its native value. It is important to select a weight in proportion to the expected score value and how much you want your results to fit the constraint; a preparatory Relax run should give you an idea of the expected score range. It is also important that the constraint function chosen reflect any uncertainty about the biological system, for which reason there exist "flat" versions of several constraint functions (LINEAR_PENALTY and HARMONI among them) which allow for any value within a range to count equivalently.

##Commonly Used Constraints
AtomPair and AtomAngle constraints are parameterized as described above, although AtomAngle constraints require three atoms (with the vertex atom listed second.) AtomPairs and AtomAngles are *not robust to changes in Rosetta numbering*; if the length of the structure is expected to change, NamedAtomPair and NamedAtomAngle constraints exist which will preserve the specific atoms identified, rather than the positions of those atoms in primary sequence space. AtomPair values are returned in Angstroms, while AtomAngle results are returned in radians.

CoordinateConstraints work like AtomPair constraints, except the second "atom" is a point in 3D space rather than an actual atom. Since Rosetta uses relative coordinates, CoordinateConstraints require a second atom to define the coordinate frame; this second atom should not be one expected to move in concert with the first, but need not be any particular distance away from it. The related LocalCoordinateConstraint uses three more atoms to define the coordinate frame. These may be used straightforwardly in Relax runs via constrain_to_native_coords.

## Commonly Used Constraint Functions
HARMONIC constraints square the distance between the ideal and actual value, and are commonly used for various types of distance constraints. CIRCULARHARMONIC is the angular equivalent.

##How to Use Constraints
As mentioned above, constraints are a way to make Rosetta's scores reflect some experimental data about the system being scored (or designed) and disfavor structures that would conflict with that data. As an example, suppose we wish to relax one half of a protein-protein interaction but we know, perhaps from mutational studies, that certain residues on each subunit interact. It may make sense to include that information via constraints, which might look like this:

AtomPair CA 356A CA 423B HARMONIC 4.3 0.25 1
AtomPair CA 432A CA 356B HARMONIC 4.3 0.25 1

To demonstrate this, run

	$>../../../main/source/bin/relax.default.linuxclangrelease -s 4eq1.pdb -out:suffix _unconstrained @general_relax_flags

and 

	$>../../../main/source/bin/relax.default.linuxclangrelease -s 4eq1.pdb -out:suffix _constrained @general_relax_flags @constraint_flags

You should see that, while the rest of each subunit moves, the N terminus of each subunit moves very little relative to residue 423, as per the constraints we entered. (You may verify this more rigorously by measuring the distance from residue 356 to residue 423 on the other subunit.) As mentioned in the Relax tutorial, if we wished to prevent any of the amino acids from moving particularly far from their starting position, we could use the option
	-relax:constrain_to_starting_coords

This demo is still under development.
The entire workflow for this demo should be described in a file
named README.dox.  It should describe an entire work flow, with
command lines, tested if possible.

The contents of each demo directory should be:

starting_files/
  -- directory in which the raw input files are given - these
     are provided for you and serve as the base for your
     tutorial
rosetta_inputs/
  -- directory in which the modified starting files should
     be placed which will be used as inputs to Rosetta.
     You may need to make modifications like stripping
     extra chains from the input PDB; store the modified
     PDB here and leave the unaltered one in starting_files 
scripts/
  -- python scripts, shell scripts, used in the workflow
  -- awk, grep, sed, etc. command lines

README.dox
  -- A prose or list description of how to perform the protocol

FOR_AUTHORS.txt
  -- A description for the demo creators of what their demo
     should achieve.
  -- Most of what starts in this file should end up in the
     README file as well.
The goal of this tutorial is to do de novo structure prediction for an
small RNA. The target structure is the SARCIN/RICIN LOOP FROM RAT 28S
R-RNA (PDB 430D). The PDB and it's fasta file are included in the
starting_files directory. 

A good starting point would be the manual:
http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/app_RNA_denovo.html
There are also prep scripts in mini/demo/rna/ which help make the
structure Rosetta friendly. They should probably be included again
here. Other good references are the integration and scientific tests:
mini/test/scientific/tests/rna_denovo
mini/test/integration/tests/rna_denovo

[TROUBLE SHOOTING]

##################################################################
Databases 
##################################################################
minirosetta_database v. rosetta_database 

Someone should double check the public release for 3.3. It appears as though the /rosetta_database/chemical directory is not present 

[!!PUBLIC RELEASE FAILS!!]

##################################################################
Residue Names
##################################################################

The PDB uses capital letters for the residue names in a .fasta file. rna_denovo requires lower case letters



##################################################################
Directory Structure of this Demo
##################################################################


rosetta_inputs/
	run_flags
			-native ./rosetta_inputs/native.pdb
			-fasta ./rosetta_inputs/4d30_.fasta
			-params_file run.prm
			-nstruct 1 [Number of predictions made]
			-out::file::silent 4ds0.out
			-cycles 1000  [30,000 cycles is necessary for larger structures]
			-minimize_rna
			-filter_lores_base_pairs
			-output_lores_silent_file
			-dump [This outputs the PDB of the final predictions]
			-mute core.io.database

	The run_flags file contains command line options

	430D_.fasta
	This is the fasta file that conatins the sequence of the RNA molecule. MAKE SURE YOU USE LOWERCASE RESIDUE NAMES

	native.pdb 
	If hte structure is known it can be used to calculate RMSDs from decoys at the end

	run_params.prm 
	Defines an chain cuts and explicit base pairs

scripts/
	run.sh 
		command used to run rna_denovo


##################################################################
Output created
##################################################################

The output generated in this demo is the 4ds0.out file contianing the score for each prediction genertated and the corresponding PDB files S_XXXXXX.pdb that contain the all atom predictions.



(c) Copyright Rosetta Commons Member Institutions.
(c) This file is part of the Rosetta software suite and is made available under license.
(c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
(c) For more information, see http://www.rosettacommons.org. Questions about this can be
(c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

Running the fixbb protocol

Design in Centroid Mode:

../../bin/fixbb.macosgccrelease -database ~/minirosetta_database/ -l list.txt -ignore_unrecognized_res -centroid_input -score:weights score3 -mute core.io core.conformation -nstruct 2

This will produce 2 output sequences for each of the pdb files listed in list.txt.

Design in Fullatom Mode:

../../bin/fixbb.macosgccrelease -database ~/minirosetta_database/ -l list.txt -ignore_unrecognized_res  -mute core.io core.conformation -nstruct 2

By default, uses standard weights with score12 patch.

Fixbb is not smart enough, yet, to reuse precomputed pair energies for multiple
designs on the same backbone.

Dock protein complex with ligand
================================

We are predicting the conformation of the complex of FKBP12, FRAP, and 
rapamycin.  Rapamycin is a dimerizer that allows FK506-binding-protein (FKBP12) 
to form an interface with FKBP-rapamycin-associated protein (FRAP). The first 
section of this tutorial demonstrates how to prepare input files for the 
proteins and small molecule. The second section describes docking rapamycin to 
FKBP12.  In the third section we will dock the result from section 2 with FRAP. 

To complete this quest you must have at least a 3 member party, including a 
cleric (level 43, must have blessing spell), a warrior proficient with hammer, 
and a thief (unlock skill level 10).

This demo was written by Gordon Lemmon, Sergey Lyskov, and Loren Looger.

Part 1: Preparing the ligand
----------------------------

Unzip the PDB file 1FAP.pdb.gz. In order for the ligand docking to work 
correctly the 1-letter chain identifier for the ligand must be different from 
the protein chain ids.  Look inside the file 1FAP.pdb.  On line 307 and 308 we 
find that chain A is FKBP12 and chain B is FRAP.  Toward the bottom of the file 
Rapamycin is specified by the residue id RAP (2375-2442).  Make a new file with 
just the RAP lines.

    grep HETATM 1FAP.pdb | grep RAP > rap.pdb

Using your favorite text editor change the chain id found in rap.pdb from 
A to X.  Use clean_pdb.py (where to find?) to prepare 1FAP.pdb for Rosetta.

    clean_pdb.py 1FAP.pdb A # this should output a file with only atom records from chain A, 1FAP_A.pdb
    clean_pdb.py 1FAP.pdb B # this should output a file with only atom records from chain B, 1FAP_B.pdb

We now have a separate PDB file for both proteins and an additional PDB file 
for the ligand. We now must add hydrogens to our rapamycin molecule and save it 
in MOL format.  Simply open the file rap.pdb in Pymol and add hydrogens using 
the action menu `all->A->hydrogens->add` (or type `h_add` on the command line). 
Then and save it as a MOL file by using the file menu (`file->save 
molecule->ok`, select type as MOL file, and change the extension to 
.mol).  You should now find the file rap.mol in your directory. 

Rosetta requires a PARAMS file for each ligand.  These files describe the 
atoms, bonds and bond angles within the ligand.  To make a params file for 
rapamycin use the script `molfile_to_params.py` found here:

    /rosetta_source/src/python/apps/public/molfile_to_params.py -c -n RAP rap.mol

We use the -c option to produce centroid mode params used in Part 3 of this 
demo.

Notice the warnings that are produced by the script.  These are informing us 
that the ligand we are using is large and flexible, which means we will 
struggle to sample all of its flexibility during docking. Since we are starting 
with the correct conformation of Rapamycin we can ignore these warnings.

mol_to_params.py should have created a file called RAP_0001.pdb which has the 
same coordinates as rap.pdb but has been prepared for use with Rosetta.  
Combine 1FAP_A.pdb and RAP_0001.pdb into a new file:

    cat 1FAP_A.pdb RAP_0001.pdb > FKBP+RAP.pdb

Part 2: Docking of proteins and ligands
---------------------------------------

Copy the `flags` and `ligand_dock.xml` files from 
`rosetta_source/src/test/integration/tests/ligand_dock_scripts` to your 
directory.  We will use these as a starting point for our docking script.

We have modified the flags file to be specific for our study. 
We now modify the options in this file to be specific for our study. First we 
comment out the start_from lines, since our ligand is already in the correct 
starting position.  Other important options to consider optimizing include the 
following.  The `angstroms` option of `Translate` should represent 
the size of your binding pocket (your ligand will move within a sphere with a 
radius of this size).

Now we are ready to run our ligand docking protocol:

    Rosetta_scripts.linuxgccrelease @flags

This should produce a file with a model of rapamycin docked to FKBP: 
`FKBP+RAP_0001.pdb`.  This file serves as an input to protein docking.

Part 3: Docking of FKBP/RAP to FRAP
-----------------------------------

Combine 1FAP_B.pdb with FKBP+RAP_0001.pdb.  Put ATOM lines from 1FAP_B first, 
followed by ATOM lines from FKBP+RAP.pdb, and then HETATM lines from 
FKBP+RAP.pdb.

    egrep 'ATOM|HETATM' 1FAP_B.pdb FKBP+RAP_0001.pdb > combined.pdb

Prepare a flag file that specifies the centroid and full-atom PARAMS files for 
rapamycin.  Also specify combined.pdb as the input file.  Run the docking 
protocol:

    docking_protocol.linuxgccrelease @flags

This should produce an output file, `combined_0001.pdb`.  Using pymol you can 
see that the FKBP/RAP complex has moved relative to FRAP.

For a production run you will want to run this protocol 10,000 or more 
times.  Then find your best scoring models. An alternative strategy would be to 
produce thousands of models with Part 1 of this tutorial, then filter for the 
top few models of FKBP with RAP.  Use each of the top models as inputs for part 
2, producing several thousand models for each of these inputs.
AbInitio Structure Prediction Using Chemical-Shift Generated Fragments and NOE Distance Restraints
==================================================================================================

Written by Lei Shi.
Nikolas Sgourakis drafted the previous version.

---

We will use the chemical shifts to improve the fragments from which Rosetta builds up structures, and the NOEs to guide the Rosetta calculations towards the native structure.

Please see references at:
* rosetta abinitio: Bradley, P et al Science 2005
* chemical shift fragments: Shen Y et al. PNAS 2008;105:4685-4690
* chemical shift+NOE+RDC: Raman S, et al Science 2010

These Rosetta calculation steps are also described separately:
* Sgourakis NG et al JACS,2011,133(16):6288-98:

In this demo, we will use PDB 2JY7, which is a small protein (for demo purpose) and has experimental data deposited. Several scripts are provided in the scripts folder for formatting purposes:

	bmrb2talos.com
	cst_map_toCB.py
	upl2mini.csh
	scores.score.cfg

If you are from David Baker lab, there are scripts available to make setup easier without going through public servers. The following instructions should work just fine without having direct access to any Baker lab cluster.

Running the demo
----------------
1. Create following folders:  
    ```
    mkdir starting_inputs
    mkdir rosetta_inputs
    mkdir rosetta_inputs/talos_output
    mkdir rosetta_inputs/pick_cs_fragments
    ```

2. Download protein fasta and experimental data
Download fasta from http://www.pdb.org/pdb/explore/explore.do?structureId=2JY7  
    ```
    wget http://www.pdb.org/pdb/files/fasta.txt?structureIdList=2JY7 -O starting_inputs/t000_.fasta
    ```
Download chemical shift data from http://www.bmrb.wisc.edu/data_library/summary/index.php?bmrbId=15591  
    ```
    wget http://rest.bmrb.wisc.edu/bmrb/NMR-STAR2/15591 -O starting_inputs/raw.cs.bmrb
    ```
Download NOE data from http://restraintsgrid.bmrb.wisc.edu/NRG/MRGridServlet?pdb_id=2JY7&show_blocks=true&min_items=0:
    ```
    wget http://restraintsgrid.bmrb.wisc.edu/NRG/MRGridServlet?db_username=wattos1&format=ambi&mrblock_id=434910&pdb_id=2jy7&program=DYANA%2FDIANA&request_type=block&subtype=general+distance&type=distance
    echo "save file as starting_inputs/NOE_data.upl"
    ```

3. Format data for Rosetta use  
Formatting NOE: (Note only residues separated by more than 3 are kept in constraint)
The script `scripts/upl2mini.csh` only works with cyana format NOE:
    ```
    scripts/upl2mini.csh starting_inputs/NOE_data.upl > rosetta_inputs/NOE.cst
    scripts/cst_map_toCB.py rosetta_inputs/NOE.cst > rosetta_inputs/NOE.centroid.cst
    ```
Formmatting chemical shift data for TALOS:
    ```
    scripts/bmrb2talos.com starting_inputs/raw.cs.bmrb > rosetta_inputs/cs.talos
    ```

4. Generating talos predictions using http://spin.niddk.nih.gov/bax/nmrserver/talosn/ using rosetta_inputs/cs.talos
Save/copy pred.tab and predSS.tab to rosetta_inputs/talos_output

5. Generate fragment/profile from RobettaServer http://www.robetta.org/fragmentqueue.jsp using starting_inputs/t000_.fasta
Save/copy t000_.checkpoint to rosetta_inputs/

6. Pick fragments using secondary structure profile and chemical shift data:
```
Rosetta/main/source/bin/fragment_picker -database Rosetta/main/database/ -in::file::vall Rosetta//tools/fragment_tools/vall.apr24.2008.extended.gz -frags::n_frags 200 -frags::frag_sizes 3 9 -frags::sigmoid_cs_A 2 -frags::sigmoid_cs_B 4 -out::file::frag_prefix rosetta_inputs/pick_cs_fragments/frags.score -frags::describe_fragments rosetta_inputs/pick_cs_fragments/frags.fsc.score -frags::scoring::config scripts/scores.score.cfg -in:file:fasta starting_inputs/t000_.fasta -in:file:checkpoint rosetta_inputs/t000_.checkpoint -in:file:talos_cs rosetta_inputs/cs.talos -frags::ss_pred rosetta_inputs/talos_output/predSS.tab talos -in::file::talos_phi_psi rosetta_inputs/talos_output/pred.tab
```

7. Run Rosetta with the fragments made above and use NOEs to guide search
```
Rosetta/main/source/bin/minirosetta -database Rosetta/main/database/ -cst_fa_file rosetta_inputs/NOE.cst -cst_file rosetta_inputs/NOE.centroid.cst -abinitio:stage1_patch scripts/patch_atom_pair_constraint -abinitio:stage2_patch scripts/patch_atom_pair_constraint -abinitio:stage3a_patch scripts/patch_atom_pair_constraint -abinitio:stage3b_patch scripts/patch_atom_pair_constraint -abinitio:stage4_patch scripts/patch_atom_pair_constraint -score:patch scripts/patch_atom_pair_constraint -in:file:fasta starting_inputs/t000_.fasta -file:frag3 rosetta_inputs/pick_cs_fragments/frags.score.200.3mers -file:frag9 rosetta_inputs/pick_cs_fragments/frags.score.200.9mers -nstruct 1 -out:file:silent csrosetta_noe.out -run:protocol abrelax -abinitio::relax -overwrite
```
You can/should adjust the weights of NOE constraints in `scripts/patch_atom_pair_constraint`.
You should also change nstruct to generate desired number of models.
Larger is better depending on your available computer time, etc.
Note that the demo in abinitio_w_chemicalshift_only, you can add flags such as:
```
    -abinitio::rg_reweight 0.5
    -abinitio::rsd_wt_helix 0.5
    -abinitio::rsd_wt_loop 0.5
    -disable_co_filter true
    -abinitio::increase_cycles 10
```
to help sampling in centroid stage.
They are not used probably NOE constraint helps guided the search.

Processing the output
---------------------
1. Extract the low energy models:
    ```
    grep SCORE csrosetta_noe.out | sort –nk2 | head
    ```
The second column contains the energies of the lowest energy 10 models.
Select as the cutoff the energy on the last line.
You should also use NOE constraint energy as a criteria to select structures.
Example is only provided for total score.

2.  This command
    ```
    cull_silent.pl csrosetta_noe.out “score < cutoff”
    ```
will produce csrosetta.select.silent which contains the lowest energy 10 models.

3. Extract pdbs from selected silent file
    ```
    Rosetta/main/source/bin/extract_pdbs -database Rosetta/main/database/ -in::file::silent csrosetta.select.silent
    ```

4. Check convergence by superimposing the ten low energy models in pymol or your favorite molecular graphics.

5. Check convergence by clustering the lowest energy models (see clustering demo for instructions).

6. To see how NOE constraints are satisfied by a model:
    ```
    Rosetta/main/source/bin/r_cst_tool.linuxgccrelease -database Rosetta/main/database/ -in:file:s lowscore_1.pdb -cst_file rosetta_inputs/NOE.cst
    ```
`r_cst_tool` is a pilot program by Oliver Lange in `Rosetta/main/source/src/apps/pilot/olli/`
# Optimize Ligand Hydroxyl Hydrogens

This tutorial assumes a unix-style command line (or cygwin on windows).

## Build your ligand
In general the ligand atoms are marked by "HETATM", but check the ligand with pymol after you've grepped out the HETATM lines:
```
grep HETATM starting_inputs/cel5A_glucan.pdb > starting_inputs/cel5A_lig_noH.pdb
```
For this step you need to add on hydrogens. We use avogadro because its open source, but you can choose other software to place hydrogens:
http://avogadro.openmolecules.net
(this example is with 1.0.3 on mac)

Open in avogadro, choose Build--> add hydrogens, then save the molecule in MDL SDfile format cel5A_lig.mol

Alternatively, you can save the molecule in PDB format and convert the PDB file into a mol file using babel (http://openbabel.org/) as follows:
```
babel -ipdb starting_inputs/cel5A_lig.pdb -omol > starting_inputs/cel5A_lig.mol
```
Or the third alternative is to open the PDB file with pymol and save the molecule with a .mol extension to force pymol to save as a mol format. 

## Make a params file for the ligand
The next step is to create a params file for rosetta. Params file contains the internal coordinates of atoms, connectivity, charge of each atom, rosetta atom type. Most importantly for this demo it contains PROTON_CHI lines which specify the proton atoms that rosetta will simple around. 
-n specifies the name of the ligand in rosetta
```
python src/python/apps/public/molfile_to_params.py starting_inputs/cel5A_lig.mol -n cel
```
The output params file is called cel.params. 
The code will sample proton chi's at the explicit values stated in the params file, and then perform a minimization on those chi's.
You shouldn't have to change the params, but if you want to add sampling explicitly you can add more angles to sample. 
A sample proton chi is:
```
CHI 1  C1   C2   O1   H8
PROTON_CHI 1 SAMPLES 3 60 -60 180 EXTRA 0
```

## Setup files for rosetta
Pull out the protein without ligand:
```
grep ATOM starting_inputs/cel5A_glucan.pdb > rosetta_inputs/cel5A_input.pdb
```
Add back in the ligand pdb from molfile to params:
```
cat cel_0001.pdb >> rosetta_inputs/cel5A_input.pdb 
```

## Run rosetta enzdes
Now we can run the enzdes app in rosetta with minimal flags.
This optimizes proton chis on the ligand while also repacking sidechains:
```
path/to/EnzdesFixBB.[platform][compiler][mode] -s rosetta_inputs/cel5A_input.pdb -extra_res_fa cel.params -database path/to/minirosetta_database/ -out:file:o cel5A_score.out -nstruct 1 -detect_design_interface -cut1 0.0 -cut2 0.0 -cut3 10.0 -cut4 12.0 -minimize_ligand true
```
This optimizes proton chis on the ligand without repacking sidechains:
```
path/to/EnzdesFixBB.[platform][compiler][mode] -s rosetta_inputs/cel5A_input.pdb -extra_res_fa cel.params -database path/to/minirosetta_database/ -out:file:o cel5A_score.out -nstruct 1 -detect_design_interface -cut1 0.0 -cut2 0.0 -cut3 0 -cut4 0 -minimize_ligand true
```
Both runs should produce the PDB file cel5A_input__DE_1.pdb and the score file cel5A_score.out, which can be placed into the output directory under different names, e.g. cel5A_output_nopack.pdb or cel5A_output_w_repack.pdb

Flag descriptions:
```
-extra_res_fa specifies the params file for new residues types (the glucan in this case)
-out:file:o specifies the file name for the enzdes-style score output. This contains extra information about the output design, like packing, interface energy, and many more
-nstruct 1 specifies one run and one output pdb; the packing is stochastic so for more sampling use a higher nstruct. 10-100 is recommended for most ligands.
-detect_design_interface tells rosetta to set up the designable and packable residues (in a packer task) based on distance from the ligand. Distances are calculated from every ligand heavy atom to the CA of amino acids.
[default values in brackets]
-cut1: CA less than cut1 is designed [6]
-cut2: CA between cut1 and cut2, with CA --> CB vector pointing towards ligand, is designed [8]
-cut3: CA less than cut3 is re-packed [10]
-cut4: CA between cut3 and cut4 with CA --> CB vector pointing towards ligand, is re-packed [12]
-minimize_ligand true  allow ligand torsions to minimize
```

## More Complete Energy Function / Sampling
If desired, use a more complete energy function and more sampling as in this example flags file
Recommended full flags for a more careful run:
```
enzdes_flags
```

## Full Enzyme Design (?)
This setup can also be used for full enzyme design
This run is very close to a full design of the active site. For a full design just change cut1 and cut2, e.g.
```
-cut1 6 -cut2 8
```
AbInitio fold-and-dock of peptides using FlexPepDock
----------------------------------------------------
This demo illustrates how to run FlexPepDock ab-initio folding and docking of a peptide onto its receptor. The FlexPepDock ab-initio protocol is designed to generate high-resolution models of complexes between flexible peptides and globular proteins, given the approximate location of the peptide binding site. The ab-initio procol samples both rigid-body orientation and torsional space of the peptide extensively. No prior knowledge about the peptide backbone is necessary, as the protocol uses fragments to sample peptide backbone conformational space rigorously.

Protocol overview
-----------------
The input to the ab-initio protocol is a model of the peptide-protein complex in PDB format,starting from arbitrary (e.g., extended) peptide backbone conformation. It is required that the peptide is initially positioned in some proximity to the true binding pocket, but the exact starting orientation may vary.
Preliminary steps: (1) Generation of fragment libraries for the peptide sequence, including 3-mer, 5-mer and 9-mer fragments. (2) Pre-packing of the receptor and peptide to remove internal clashes that might confuse ranking.
Main protocol: Step 1: Monte-Carlo simulation for de-novo folding and docking of the peptide over the protein surface in low-resolution (centroid) mode, using a combination of fragment insertions, random backbone perturbations and rigid-body transformation moves. Step 2: The resulting low-resolution model is refined with FlexPepDock Refinement. As in the independent refinement protocol, the output models are then ranked based on their energy score, after their clustering for improved coverage of distinct conformations.

Refinement vs. ab-initio protocol
---------------------------------
The Refinement protocol is intended for cases where an approximate, coarse-grain model of the interaction is available that is close to the correct solution both in Cartesian and dihedral (phi, psi) space. The protocol iteratively optimizes the peptide backbone and its rigid-body orientation relative to the receptor protein including on-the-fly side-chain optimization (look at /demos/refinement_of_protein_peptide_complex_using_FlexPepDock/ to learn how to run refinement of protein-peptide complexes).
The ab-initio protocol extends the refinement protocol considerably, and is intended for cases where no information is available about the peptide backbone conformation. It simultaneously folds and docks the peptide over the receptor surface, starting from any arbitrary (e.g., extended) backbone conformation. It is assumed that the peptide is initially positioned close to the correct binding site, but the protocol is robust to the exact starting orientation. The resulting low-resolution models are refined using the FlexPepDock Refinement protocol.

Running the FlexPepDock ab-initio protocol
------------------------------------------
1. Generate an initial complex structure: An initial model can be built by placing the peptide in close promity to the binding site in an arbitary conformation. In this demo we have provided a starting structure with a peptide in extended conformation (2A3I.ex.pdb). Our goal is to optimize this structure using ab-initio FlexPepDock, towards a near-native model with a helical peptide conformation. Both the native structure (2A3I.pdb), as well as the starting structure (2A3I.ex.pdb) are provided in the input directory.

2. Prepack the input model: This step involves the packing of the side-chains in each monomer to remove internal clashes that are not related to inter-molecular interactions. The prepacking guarantees a uniform conformational background in non-interface regions prior to refinement. The prepack_flags file contains the flags for running the prepacking job. The run_prepack script will run prepacking of the input structure 2A3I.ex.pdb located in the input directory.

You need to change the paths of the Rosetta executables and database directories in the run_prepack script (also for run_refine; see below).

  ROSETTA_BIN="rosetta/main/source/bin"
  ROSETTA_DB="rosetta/main/database/"

After changing the paths, run the run_prepack script as:
   $./run_prepack

The output will be a prepacked structure, 2A3I.ex.ppk.pdb, located in the input directory; a scorefile named ppk.score.sc and a log file named prepack.log file located in the output directory. This prepacked structure will be used as the input for the ab-initio modeling step.

3. Create 3mer, 5mer & 9mer (peptide lingth >=9) fragment libraries: The scripts necessary for creating fragments are provided in the fragment_picking directory.
    a. Go to the fragment_picking directory.
    b. Save the peptide sequence in the xxxxx.fasta file.
    c. Run the make_fragments.pl script to generate the PSIPred secondary structure and PSI-Blast sequence profiles. You need to chnage the paths in the upper section of the make_fragments.pl file.
       Run as $perl make_fragments.pl -verbose -id xxxxx xxxxx.fasta
       This will create xxxxx.psipred_ss2, xxxxx.checkpoint along with other files.
    d. Run the executable fragment_picker.linuxgccrelease to create the frags. The flags are provided in the flags file and fragment scoring weights are provided in the psi_L1.cfg file.
    Run as $ROSETTA_BIN/fragment_picker.linuxgccrelease -database $ROSETTA_DB @flags >log
    e. Change the fragment numbering using shift.sh script.
    Run as $bash shift.sh frags.500.3mer X >frags.3mers.offset ; where X is the number of residues in the receptor. Do the same for 5mer and 9mer frags
    The offset fragment files will be used as input to the FlexPepDock ab-initio protocol. Put them in the input/frags directory.

4. Ab-initio folding and docking of the prepacked model: This is the main part of the protocol. In this step, the peptide backbone and its rigid-body orientation are optimized relative to the receptor protein using the Monte-Carlo with Minimization approach, including periodic on-the-fly side-chain optimization. The peptide backbone conformational space is extensively sampled using fragments derived from solved structures. The file abinitio_flags contains flags for running the ab-initio job. The run_abinitio script will run ab-initio modeling of the prepacked structure generated in the prepacking step located in the input directory.

After changing the Rosetta related paths run the run_abinitio script as:
    $./run_abinitio

The output will be an optimized structure (2A3I.ex.ppk_0001.pdb) located in the output directory; a scorefile named abintio.score.sc and a log file named abinitio.log, located in the output directory. This script has to be modified to run on a cluster during a production run (see below).


Specific changes needed for a production run
--------------------------------------------
For a production run it is recommended to generate large number of decoys (~10,000 to 50,000). In such a case you can run the job on a cluster. It is advided to use silent output format in such scenario to save space (See https://www.rosettacommons.org/manuals/rosetta3.1_user_guide/app_silentfile.html for details). Include the following lines to the abinitio_flags file:

-out:file:silent_struct_type binary
-out:file:silent decoys.silent

This will create the decoys.silent file containing data related to all the decoys in a compressed format. You can extract speicific decoy using the extract_pdbs.linuxgccrelease executable.
For example:
  $ROSETTA_BIN/extract_pdbs.linuxgccrelease -database $ROSETTA_DB -in:file:silent decoys.silent -in:file:tags 2A3I.ex.ppk_1234.pdb


Along with changes in the flags file you need to modify the run_abinitio file to run on a cluster. The file run_abinitio_slurm is the modified version of run_abinitio adapted to run on a slurm cluster. You should ask your cluster manager for relevent changes required.


Post Processing after a production run
--------------------------------------
In order to diversify our prediction, we cluster the results and select representative models. A clustering scripts is provided in the clustering directory. It will cluster the top 500 decoys based on a cutoff radius of 2.0 Angstrom, and select for each the top-scoring member (according to reweighted score,  reweighted_sc). A top scoring member from each cluster is reported in the file  cluster_list_reweighted_sc_sorted.

Runs as
$bash cluster.sh 2.0 ../input/2A3I.ex.pdb ../output/decoys.silent

Further information
-------------------
Detailed documentation on ab initio FlexPepDock is available under: https://www.rosettacommons.org/docs/latest/application_documentation/docking/flex-pep-dock.
Please cite: Raveh B, London N, Zimmerman L, Schueler-Furman O (2011) Rosetta FlexPepDock ab-initio: Simultaneous Folding, Docking and Refinement of Peptides onto Their Receptors. PLoS ONE 6(4): e18934. doi: 10.1371/journal.pone.0018934

The scripts and input files that accompany this demo can be found in the 
`demos/` directory of the Rosetta weekly releases.
# Predict eglinC DDGs
The entire workflow for this demo should be described in a file
named README.dox.  It should describe an entire work flow, with
command lines, tested if possible.

```
starting_files/
  -- directory in which the raw input files are given - these
     are provided for you and serve as the base for your
     tutorial
rosetta_inputs/
  -- directory in which the modified starting files should
     be placed which will be used as inputs to Rosetta.
     You may need to make modifications like stripping
     extra chains from the input PDB; store the modified
     PDB here and leave the unaltered one in starting_files 

scripts/
  -- python scripts, shell scripts, used in the workflow
  -- awk, grep, sed, etc. command lines
  extract_chains.pl
	-extracts specified chain of pdb and prints to standard out. execute without arguments for usage.
  sequentialPdbResSeq.pl
	-renumbers specified pdb and prints to standard out. execute without arguments for usage.
  (we didn't use the scripts topN_average.scr or compute_top3avg_energies.scr)

README.dox
  -- A prose or list description of how to perform the protocol

FOR_AUTHORS.txt
  -- A description for the demo creators of what their demo
     should achieve.
  -- Most of what starts in this file should end up in the
     README file as well.
```

## Running
The given starting structure was called starting_files/1CSE.pdb.

We started by:
1. removed the chain we didn't want to do ddG calculations on. the command for doing so is:
    ```
	./scripts/extract_chains.pl starting_files/1CSE.pdb I  > starting_files/1CSEi.pdb
    ```
    this outputs chain I to the file starting_files/1CSEi.pdb

2. renumbered the crystal structure starting from 1.
    ```
   ./scripts/sequentialPdbResSeq.pl  -pdbfile 1CSEi.pdb -res1 1 > 1CSEi.ren.pdb
    ```
    this script takes chain I of 1CSE and renumbers starting from 1, then outputs to 1CSEi.ren.pdb
    **WARNING:** if your pdb has a chainbreak (missing part of the poly-peptide chain), then your numbering will be inconsistent. For example, if you have a chain-break between 12 and 23, it will be renumbered as : 12 13 and so on..

3. minimized the starting input file with harmonic constraints on all C-alpha atoms within 9 Angstrom rmsd.
   this must be run from the directory which contains the pdb-file, otherwise you might get an error. 
he minimization protocol only takes in lists of files, so you need to do the following:
    ```
    cd starting_files/
    ls 1CSEi.ren.pdb > lst
    /rosetta_release_3/rosetta-3.3/rosetta_source/bin/minimize_with_cst.default.macosgccrelease -in:file:l lst -database ~/minirosetta_database/ -in:file:fullatom -ddg::out_pdb_prefix minimize_with_cst        
    ```

4. if you want to double-check that the minimization worked, you can score the structures as follows:
    ```
    ls 1CSEi.ren.pdb > test.lst 
    ls minimize_with_cst.1CSEi.ren_0001.pdb >> test.lst
    ~/rosetta_release_3/rosetta-3.3/rosetta_source/bin/score.default.macosgccrelease -in:file:l test.lst -database ~/minirosetta_database/ -in:file:fullatom -out:file:scorefile score.chk.fsc 
    ```
    and the score for minimize_with_cst.1CSEi.ren_0001.pdb should be lower than 1CSEi.ren.pdb.
    In this case, 1CSEi.ren.pdb has a score of 35.294 and minimize_with_cst.1CSEi.ren_0001.pdb has a score of -64.426. And for a given input structure you should always converge on a score (you should get the same score for each minimized-input structure).

5. prepare the mutation file:
    mutations should be in the form:
    (Wild-type-residue)(residue-position)(mutant-residue)
    with no spaces in between.
    For example, we have the file: mutations.multiples.txt which has the contents:
    ```
    V13A
    V14G
    V18G
    
    A21F
    E23A
    
    F25A
    ```
    output will be as follows:
    ```
    total 6
    3
    V 6 A
    V 7 G
    V 11 G
    2
    A 14 F
    E 16 A
    1
    F 18 A
    ```

    the new-lines mean that this is a complete batch of mutations. ( In this case we would make a triple mutant: V 13 -> A, V 14 -> G, and V 18 -> G. We would also make a double mutant: A 21 -> F, and E 23 -> A. And finally, we make the single mutant: F 25 -> A).
    **REMEMBER:** keep track of the offset between your initial pdb and renumbered pdb.
    To format the mutations.multiples.txt for input into ddgs run the script mutation_format_ddgs.pl as follows:
    ```
    perl mutation_format_ddgs.pl mutations.multiples.txt offset output_path.mut
    ```
    Where offest is the offset used to prepare the pdb for input into ddgs and output_path is the title you want for the .mut file. The output path is an optional variable. If no output path is provided, the script will print to standard output.

    The full explanation for the mutation format is here:
    http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/app_ddg_prediction.html
 
    Basically, the total keyword specifies how many *total* mutations you want to make. For example, if you make a triple mutant, and double mutant, and a single mutant, the value for total will be 6. Then in order to specify how many mutations you make in each 'batch' you specify with a number followed by a newline. 


6. run the ddg prediction application as follows: 
    ```
    ~/rosetta_release_3/rosetta-3.3/rosetta_source/bin/ddg_monomer.default.macosgccrelease -in:file:s minimize_with_cst.1CSEi.ren_0001.pdb -ddg::weight_file soft_rep -ddg::iterations 5 -ddg::dump_pdbs true -ddg::mut_file mutations.multiples.txt -database ~/minirosetta_database/ -ddg::local_opt_only false -ddg::min_cst false -ddg::mean true -ddg::min -ignore_unrecognized_res 
    ```
    This repacks the wild-type and the mutant structures 5 times, and in order to compute the ddG (which is Emutant - Ewt ) it averages the scores of the mutant and wild-type ensembles of structures.  It uses the soft_rep_design scoring function and it is a fixed backbone protocol. Normally I would run this protocol for 20 iterations , but for the demo purposes I'm only testing with 5 iterations.
    The full explanation of all options is here:
    http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/app_ddg_prediction.html
	The reason we choose these set of options is that it gives pretty reliable results in terms of correlation and stability predictions (in terms of destabilizing, neutral, or stabilizing) and it is relatively quick to run.     

    It dumps out the pdbs in the format:
     repacked_wt_round_X.pdb (where X is from 1-5)
     mut_(mutation_label).pdb where the mutation is listed for the triple mutant (for example) as mut_V6AV7GV11G_round_X.pdb (where X is from 1-5)

    Most importantly, it dumps out ddg-predictions in the file: ddg_predictions.out
     If you take a look at the file you will see a header:
    ```
    ddG: description total fa_atr fa_rep fa_sol fa_intra_rep pro_close fa_pair hbond_sr_bb hbond_lr_bb hbond_bb_sc hbond_sc dslf_ss_dst dslf_cs_ang dslf_ss_dih dslf_ca_dih fa_dun p_aa_pp ref 
    ddG: V6AV7GV11G     8.842    16.575    -2.875    -1.825    -1.496     0.000     0.029     0.000     0.000    -0.089     0.446     0.000     0.000     0.000     0.000    -0.557     0.564    -1.930 
    ```

    The description corresponds to the mutation information, the total is the predicted overall ddg. Positive means the mutation was destabilizing, and negative means the mutation was stabilizing. The numbers after the total correspond to the components that contribute to the overall score. Adding them, should sum to the total. 

     
# Stepwise Monte Carlo RNA Multiloop

## Authors
Rhiju Das, rhiju@stanford.edu

## StepWise Monte Carlo (examples for RNA)

### Brief Description

Solve structure of an RNA internal loop or multi-helix junction.

### Abstract

Ab initio and comparative modeling of biopolymers (RNA, protein, protein/RNA) often involves solving well-defined small puzzles (4 to 20 residues), like RNA aptamers, RNA tertiary contacts, and RNA/protein interactions. If these problems have torsional combinations that have not been seen previously or are not captured by coarse-grained potentials, most Rosetta approaches will fail to recover their structures.  This app implements a stepwise ansatz, originally developed as a 'stepwise assembly' enumeration that was not reliant on fragments or coarse-grained modeling stages, but was computationally expensive. The new mode is a stepwise monte carlo, a stochastic version of stepwise assembly. 


## Running

Following is for an internal loop ('two-way junction') drawn from the most conserved domain of the signal recognition particle, a core component of the machinery that translates membrane proteins in all kingdoms of life.

If you do not know the rigid body orientations of two helices (typical use case), run:

```
stepwise -in:file:fasta rosetta_inputs/1lnt.fasta -s rosetta_inputs/gu_gc_helix.pdb  rosetta_inputs/uc_ga_helix.pdb -out:file:silent swm_rebuild.out -extra_min_res 2 15 7 10 -terminal_res 1 8 9 16 -nstruct 20  -cycles 1000  -score:rna_torsion_potential RNA11_based_new  -native rosetta_inputs/native_1lnt_RNA.pdb
```

If you have starting coordinates for the two helix endpoints, you can start with that single PDB ('start_native_1lnt_RNA.pdb') instead:

```
stepwise -in:file:fasta rosetta_inputs/1lnt.fasta -s rosetta_inputs/start_native_1lnt_RNA.pdb -out:file:silent swm_rebuild.out -extra_min_res 2 15 7 10 -terminal_res 1 8 9 16 -nstruct 20  -cycles 1000  -score:rna_torsion_potential RNA11_based_new  -native rosetta_inputs/native_1lnt_RNA.pdb
```

To get out models:

```
extract_pdbs -silent swm_rebuild.out 
```

(Or use extract_lowscore_decoys.py which can be installed via tools/rna_tools/.)

Domain Insertion Demo
=====================

In this tutorial, we will demonstrate how to perform domain insertion with 
Rosetta.  Domain insertion is when you have two well-folded domains, and you 
insert one (B) into a flexible loop of the other (A), such that domain A is 
split into halves in primary sequence, but the whole protein still folds into A 
and B.  We will also allow for redesign of the remodeled loop accepting the 
insertion.

The sample PDBs are 1EMA (GFP) and 2LCT (an SH2 domain).

A dirty little secret is that there is no proper way to do domain insertion in 
Rosetta.  Nobody's gotten around to writing a real mode for it.  We are instead 
performing what is technically known as an "epic hack" to use a totally 
different suite, AnchoredDesign, for the purpose.  Unfortunately, there are 
many strange nomenclature issues forced by this: the "anchor" is the inserted 
domain, and the "scaffold" is the domain receiving the insertion, and the 
"target" is non-existent for domain insertion but must pretend-exist for the 
purpose of the code.  `AnchoredDesign` and `AnchoredPDBCreator` are 
[[extensively documented|public/anchored_design/README]] with a 
protocol capture released in 3.3, please refer to that documentation for more 
details.

Where the insertion will occur
------------------------------

Open the structures in the Pymol. The goal is to insert SH2 domain into GFP. 
Presumably you will have a scientific problem where you are either interested 
in a particular insertion isomer/position, or you only care that it's a good 
structure but not exactly where the insertion occurs.

We identify the loop region 207-220 (chain A) in GFP arbitrarily as the target 
for insertion.  It would be good to sample many insertion positions and loop 
lengths to find the MOST stable insertion; that is beyond the scope of this 
tutorial but should be easy by extension.

Preparing input structures for the first step, AnchoredPDBCreator
-----------------------------------------------------------------

We will have to do some manual editing to get the input PDBs ready.  Some of 
these tweaks are just Rosetta idiosyncrasies, some are AnchoredDesign issues.

- Prepare the SH2 pdb
    - Delete all the PDB file head matter, up to the first ATOM record.
    - Delete all lines including and after the first ENDMDL card.  This happens 
      to be an NMR model, and we only want the first NMR sub-model in the PDB, 
      and it's best to just delete them now.
    - Delete the peptide chain B out of the pdb.  It is irrelevant to this 
      problem.
    - To make the loop closures more plausible, we chose to delete residues 
      661-668 out of the SH2 domain.  This brings its termini closer together, 
      so that the loop accepting the insertion will not have to stretch to 
      accommodate it.
    - Compare 2lct.pdb to 2lct_prepared.pdb to see these changes.

- Prepare the GFP pdb
    - We will use it as is.  It has a fluorophore (duh), we'll use the flag 
      `-ignore_unrecognized_res` to silently edit it out of the input.  It is 
      not near the surface and won't affect the modeling.

- Prepare the “target” pdb
    - The target has no meaning for domain insertion, but is a requirement of 
      the AnchoredDesign suite.  Here, take a glycine from the SH2 domain and 
      move it far from the SH2 domain.  The easy way to do this is to copy a 
      single GLY residue into a new file and translate it by manually adding 
      many (900) angstroms to its x/y/z coordinates.  You can compare 
      `rosetta_inputs/AnchoredPDBCreator/pseudotarget.pdb` to residue 691 of 
      the SH2 domain to see what we did.

Performing step 1, AnchoredPDBCreator
-------------------------------------

AnchoredPDBCreator will perform the mechanical part of the insertion operation 
(actually suturing the sequences together), but it does not attempt to model 
the new interface much.  We will run it briefly to get a rough inserted 
structure which we will refine later.

To run AnchoredPDBCreator, use its executeable (of the same name):

    cd <demo directory>/rosetta_inputs/AnchoredPDBCreator
    <path/to/rosetta>/bin/AnchoredPDBCreator.<yoursystemsettings> @options

This command will create two new files in that directory.  The output structure 
will be named `S_0001.pdb`; there will also be a scorefile `score.sc`.

In realistic usage, you would generate several hundred models and choose a 
subset to subject to more processing by analyzing the LAM score (reported in 
both the PDB file and score.sc).  Also note that you will increase the value of 
the APDBC_cycles argument to result in longer trajectories.  LAM score 
(LoopAnalyzerMover) attempts to capture how well-closed and formed the subject 
loop is; it is further explained in the AnchoredDesign documentation.  
AnchoredPDBCreator results need only be judged on that criterion; the 
AnchoredDesign protocol will refine it anyway.

To sort through many models, run scripts/sort_by_LAM.sh in the result folder.  
It will sort the scorefile to put the best (lowest) LAM scores at the top.  
Manually examine the best handful and pick your favorite.  It will be used as 
the input to the next step, AnchoredDesign.  Don't stress over your choice here 
– it doesn't have to be a great structure, it's an input not an output.

Creating inputs for AnchoredDesign
----------------------------------

The primary input to AnchoredDesign is the result from AnchoredPDBCreator, 
`S_0001.pdb`.  You will want to load it up in a viewer of your choice for the 
next step.  Note that all numbering from this point on is relative to 
`S_0001.pdb`, NOT the original PDBs.  Also note it is chain B, not chain A; 
chain A is the pseudotarget.

- Creating an anchor file

  The anchor file tells AnchoredDesign what the rigid inserted region is (the 
  insert domain).  Here, it is residues 213-305 in chain B of S_0001.pdb; that 
  is what used to be 2lct_prepared.pdb.  The anchor file is formatted B 213 
  305, see AnchoredDesign documentation for more details.

- Creating a loops file.

  The loops file tells AnchoredDesign what regions are flexible loops.  It will 
  actually treat what used to be one loop plus the insertion as one huge loop, 
  but leave the insertion rigid.  So, our loops file will specify a loop 
  running from the N-terminus of the insert loop to the C-terminus, going 
  through the whole insert domain.  The file is at 
  rosetta_inputs/AnchoredDesign/loopfile; the loop file format documentation is 
  in the manual.

- Creating a resfile

  The resfile is optional; it will allow you to mutate residues in the loop to 
  design a loop that best accepts the insertion.  (If you do not use a resfile, 
  use the flag -packing:repack_only instead to preclude design).  Resfile 
  format documentation is available in the manual.  In our resfile, we have 
  specified that the flexible positions in the loop (the loop, but not the 
  inserted domain) can be designed to any residue.

Running AnchoredDesign to refine the insertion
----------------------------------------------

To run AnchoredDesign, use its executable (of the same name):

    cd <demo directory>/rosetta_inputs/AnchoredDesign
    <path/to/rosetta>/bin/AnchoredDesign.<yoursystemsettings> @options

AnchoredDesign will remodel the loop containing the insertion and sample the 
pseudo-rigid-body degree of freedom between the SH2 and GFP, while leaving the 
cores of each domain rigid.  It will also (optionally) design the loop region 
to create a loop that best accepts the insertion.

In this tutorial, the settings nstruct, refine_cycles, and perturb_cycles are 
set fairly low for speed.  In production, you will want to turn these flags up 
higher for better results; please see the options file for more details.

Interpreting results
--------------------

(The AnchoredDesign results will contain a meaningless chain A from the 
pseudotarget – delete or ignore it at your leisure.  Also, AnchoredDesign's 
interface metrics refer to the chain A – chain B interface, which won't exist; 
you should ignore those too.)

AnchoredDesign will create PDB files (here of the form S_0001_*.pdb) and a 
scorefile, score.sc.  Interpreting the results requires all your scientific 
intuition.  For a first pass, you can sort the models by total score with the 
command sort_by_score.sh in the scripts directory.  As before, low scores are 
better.  You can examine the other score terms (reported in the score file), 
and even per-residue scores (reported at the end of each PDB), to help you 
decide which model you think is most physically plausible.  You will want to 
examine your models individually in a viewer to pick the best.  Look for 
well-formed interfaces between the two domains and well-closed loops with good 
geometry.
Dock protein complex with ligand
================================

We are predicting the conformation of the complex of FKBP12, FRAP, and 
rapamycin.  Rapamycin is a dimerizer that allows FK506-binding-protein (FKBP12) 
to form an interface with FKBP-rapamycin-associated protein (FRAP). The first 
section of this tutorial demonstrates how to prepare input files for the 
proteins and small molecule. The second section describes docking rapamycin to 
FKBP12.  In the third section we will dock the result from section 2 with FRAP. 

To complete this quest you must have at least a 3 member party, including a 
cleric (level 43, must have blessing spell), a warrior proficient with hammer, 
and a thief (unlock skill level 10).

This demo was written by Gordon Lemmon, Sergey Lyskov, and Loren Looger.

Part 1: Preparing the ligand
----------------------------

Unzip the PDB file 1FAP.pdb.gz. In order for the ligand docking to work 
correctly the 1-letter chain identifier for the ligand must be different from 
the protein chain ids.  Look inside the file 1FAP.pdb.  On line 307 and 308 we 
find that chain A is FKBP12 and chain B is FRAP.  Toward the bottom of the file 
Rapamycin is specified by the residue id RAP (2375-2442).  Make a new file with 
just the RAP lines.

    grep HETATM 1FAP.pdb | grep RAP > rap.pdb

Using your favorite text editor change the chain id found in rap.pdb from 
A to X.  Use clean_pdb.py (where to find?) to prepare 1FAP.pdb for Rosetta.

    clean_pdb.py 1FAP.pdb A # this should output a file with only atom records from chain A, 1FAP_A.pdb
    clean_pdb.py 1FAP.pdb B # this should output a file with only atom records from chain B, 1FAP_B.pdb

We now have a separate PDB file for both proteins and an additional PDB file 
for the ligand. We now must add hydrogens to our rapamycin molecule and save it 
in MOL format.  Simply open the file rap.pdb in Pymol and add hydrogens using 
the action menu `all->A->hydrogens->add` (or type `h_add` on the command line). 
Then and save it as a MOL file by using the file menu (`file->save 
molecule->ok`, select type as MOL file, and change the extension to 
.mol).  You should now find the file rap.mol in your directory. 

Rosetta requires a PARAMS file for each ligand.  These files describe the 
atoms, bonds and bond angles within the ligand.  To make a params file for 
rapamycin use the script `molfile_to_params.py` found here:

    /rosetta_source/src/python/apps/public/molfile_to_params.py -c -n RAP rap.mol

We use the -c option to produce centroid mode params used in Part 3 of this 
demo.

Notice the warnings that are produced by the script.  These are informing us 
that the ligand we are using is large and flexible, which means we will 
struggle to sample all of its flexibility during docking. Since we are starting 
with the correct conformation of Rapamycin we can ignore these warnings.

mol_to_params.py should have created a file called RAP_0001.pdb which has the 
same coordinates as rap.pdb but has been prepared for use with Rosetta.  
Combine 1FAP_A.pdb and RAP_0001.pdb into a new file:

    cat 1FAP_A.pdb RAP_0001.pdb > FKBP+RAP.pdb

Part 2: Docking of proteins and ligands
---------------------------------------

Copy the `flags` and `ligand_dock.xml` files from 
`rosetta_source/src/test/integration/tests/ligand_dock_scripts` to your 
directory.  We will use these as a starting point for our docking script.

We have modified the flags file to be specific for our study. 
We now modify the options in this file to be specific for our study. First we 
comment out the start_from lines, since our ligand is already in the correct 
starting position.  Other important options to consider optimizing include the 
following.  The `angstroms` option of `Translate` should represent 
the size of your binding pocket (your ligand will move within a sphere with a 
radius of this size).

Now we are ready to run our ligand docking protocol:

    Rosetta_scripts.linuxgccrelease @flags

This should produce a file with a model of rapamycin docked to FKBP: 
`FKBP+RAP_0001.pdb`.  This file serves as an input to protein docking.

Part 3: Docking of FKBP/RAP to FRAP
-----------------------------------

Combine 1FAP_B.pdb with FKBP+RAP_0001.pdb.  Put ATOM lines from 1FAP_B first, 
followed by ATOM lines from FKBP+RAP.pdb, and then HETATM lines from 
FKBP+RAP.pdb.

    egrep 'ATOM|HETATM' 1FAP_B.pdb FKBP+RAP_0001.pdb > combined.pdb

Prepare a flag file that specifies the centroid and full-atom PARAMS files for 
rapamycin.  Also specify combined.pdb as the input file.  Run the docking 
protocol:

    docking_protocol.linuxgccrelease @flags

This should produce an output file, `combined_0001.pdb`.  Using pymol you can 
see that the FKBP/RAP complex has moved relative to FRAP.

For a production run you will want to run this protocol 10,000 or more 
times.  Then find your best scoring models. An alternative strategy would be to 
produce thousands of models with Part 1 of this tutorial, then filter for the 
top few models of FKBP with RAP.  Use each of the top models as inputs for part 
2, producing several thousand models for each of these inputs.
# SWA Protein Main

## Author
Rhiju Das, rhiju@stanford.edu

# Loop Remodeling by Enumeration: the Core Step of Protein 'Stepwise Assembly'

## Brief Description

Build a loop denovo by enumerating through phi,psi angles, and closing the chain by CCD. Should give a complete enumeration for a loop up to 5 residues in length.

## Abstract

Consistently predicting protein structure at atomic resolution from sequence alone remains an unsolved problem in computational biophysics. Practical challenges involving protein loops arise frequently in ab initio modeling, comparative modeling, and protein design, but even these cases can become intractable as loop lengths exceed 10 residues and if surrounding side-chain conformations are erased. This demo illustrates a novel approach to protein modeling that is more powerful than prior methods that strive for atomic resolution. The central innovation is a ‘stepwise ansatz’ inspired by recent ab initio RNA algorithms, which resolves a conformational sampling bottleneck through residue-by-residue conformer enumeration and dynamic programming.


Reference: R. Das (2013) "Atomic-accuracy prediction of protein loop structures enabled by an RNA-inspired ansatz", under review.
More info: http://arxiv.org/abs/1208.2680

## Example Rosetta Command Line

This rebuilds residues 5-8 on the knottin scaffold 2it7 which has had the loop and all the protein's sidechains removed:

```
swa_protein_main -rebuild  -s1 rosetta_inputs/noloop5-8_2it7_stripsidechain.pdb   -input_res1 1-4 9-28   -sample_res 5 6  -bridge_res 7 8  -cutpoint_closed 7    -superimpose_res 1-4 9-28  -fixed_res 1-4 9-28   -calc_rms_res 5-8  -jump_res 1 28  -ccd_close  -out:file:silent_struct_type binary  -fasta rosetta_inputs/2it7.fasta  -n_sample 18  -nstruct 400  -cluster:radius 0.100  -extrachi_cutoff 0  -ex1  -ex2  -score:weights score12.wts  -pack_weights pack_no_hb_env_dep.wts  -in:detect_disulf false  -add_peptide_plane  -native rosetta_inputs/2it7.pdb  -mute all   -out:file:silent 2it7_rebuild.out -disulfide_file rosetta_inputs/2it7.disulf
```

If you plot column 26 against 1 in 2it7_rebuild.out you should see that the lowest energy models have backbone rmsds of well under 1 Angstrom.

```
%(bin)s/extract_pdbs.%(binext)s -in:file:silent 2it7_rebuild.out   -in:file:silent_struct_type binary  -in:file:tags S_0 -database %(database)s -run:constant_seed -nodelay  2>&1 \
```

**NOTE:** Running 'StepWise Assembly' on a longer loop requires a more complex workflow that carries out buildup of the loop across all possible residue-by-residue build paths. This requires a master python script to setup the job, and another master python script to queue up the resulting computation, which is described as a directed acyclic graph. This full workflow is being presented in a separate demo [swa_protein_long_loop].

## Versions
This should work directly out of trunk for any version of Rosetta after June 2012; however for versions before March 2013, the name of the "swa_protein_main" application was "stepwise_protein_test".


Chemically Conjugated Docking
=============================

Included are demos for three applications:

* [[UBQ_E2_thioester|public/chemically_conjugated_docking/UBQ_E2_thioester/readme]]
* [[UBQ_Gp_CYD-CYD|public/chemically_conjugated_docking/UBQ_Gp_series/readme]]
* [[UBQ_Gp_LYX-Cterm|public/chemically_conjugated_docking/UBQ_Gp_series/readme]]

All are closely related.
These demos are copies of their integration tests.
See the linked READMEs for further details.
This demo contains the starting structures for the original published use of UBQ_E2_thioester, and is a copy of UBQ_E2_thioester's integration test.  To run the demo (in the inputs subfolder):

    UBQ_E2_thioester.linuxgccrelease @options

Edit the options file first to set paths as needed.

Please refer to UBQ_E2_thioester's documentation, online or at rosetta_source/doc/apps/public/scenarios/UBQ_conjugated.dox, and also the publication, Saha A, Kleiger G, Lewis S, Kuhlman B, Deshaies RJ. Essential role for ubiquitin-ubiquitin-conjugating enzyme interaction in ubiquitin discharge from Cdc34 to substrate. Molecular Cell. 2011 Apr 8;42(1):75-83.

Note that the provided outputs are from the integration test, which runs in ~30 s.  You will need to run the code for much longer (both longer individual runs, and many trajectories) to get scientifically useful results.  The options file and documentation provide details on how to do that.
This demo contains the starting structures for the original use of the UBQ_Gp series of executables, and is a copy of their integration tests.  To run the demos (in the inputs subfolder):

    UBQ_Gp_CYD-CYD.linuxgccrelease @options
    UBQ_Gp_LYX-Cterm.linuxgccrelease @options

Each uses the same input files.  Read the options file first to set local paths as needed.

Please refer to the documentation, online or at rosetta_source/doc/apps/public/scenarios/UBQ_conjugated.dox, and also the publication, Baker R, Lewis SM, Wilkerson EM, Sasaki AT, Cantley LC, Kuhlman B, Dohlman HG, Campbell SL.  Site-Specific Monoubiquitination Activates Ras by Impeding GTPase Activating Protein Function.  Submitted.

Note that the provided outputs are from the integration test, which runs in ~30 s.  You will need to run the code for much longer (both longer individual runs, and many trajectories) to get scientifically useful results.  The options file and documentation provide details on how to do that.
# RNA Design

This code is intended to carry out fixed backbone design of RNA sequences given an input backbone.

This demo redesigns a 'UUCG' tetraloop on a single-base pair RNA 'helix', as a small 6-nucleotide test case. As illustration, only 3 designs are output. It takes about 15 seconds to run. Run:

```
  rna_design.linuxgccrelease @flags -database <PATH TO ROSETTA DATABASE>
```

The output will show up in:

```
 chunk001_uucg_RNA.pack.txt 
```

with scores in

```
 chunk001_uucg_RNA.pack.out
```

and PDBs in S_0001.pdb, S_0002.pdb, etc.  
The typical sequence output is cuuggg (native is cuucgg). 
# Tutorial for using modified residue types in centroid-level applications.

In this tutorial, you will describe how to wrangle modified residue types through combined centroid/fullatom protocols in Rosetta.

##  Background

1. Many Rosetta protocols use a centroid phase followed by a fullatom phase.  However, the centroid residue set is compatible with the canonical 20 amino acids, and very little else.  For example, ligands and post-translational modifications are incompatible with centroid mode; Rosetta generally crashes when trying to switch these residues to centroid (or silently drops the post-translational modification).

2. You have been provided with a PDB file containing a two phosphoresidues as input (one p-SER and one p-TYR). This tutorial will show you how run the docking\_protocol application using this PDB model as input without Rosetta exiting with an error.

## Clean up the PDB file

1. This is an NMR model. Rosetta will read this as a giant complex of superimposed structures, which is bad. You need to find the line labeled "ENDMDL" and delete all lines below it.
    ```
    gunzip -c starting_files/2lax.pdb.gz > rosetta_inputs/2lax_edited.pdb
    <your favourite text editor> rosetta_inputs/2lax_edited.pdb
    ```

2. Rosetta wants the phosphorylated residues to be named the same as their canonical countertypes. Open the pdb file rosetta\_inputs/2lax\_edited.pdb for editing. For residue 202 and residue 206, rename the three-letter code "TPO" to "TYR" and "SEP" to "SER". Rosetta will identify these as phosphorylated based on their atom names.

    - note: If at this point you were to run any protocol that utilizes centroid mode, you would see this error when Rosetta crashed:
        ```
        can not find a residue type that matches the residue TYR_p:phosphorylatedat position 38

        ERROR: core::util::switch_to_residue_type_set fails
        ```

## Create centroid residue parameter patch files

For modified residues (e.g. phosphorylation or acetylation), Rosetta uses "patch" files to modify the pre-existing residue parameter files. Unfortunately, a centroid-level patch file for phosphorylated residues does not exist. You need to create two, one for P-TYR and one for P-SER.
1. Go to the database directory ```<my_rosetta_directory>/rosetta_database/chemical/residue_type_sets/centroid/patches/``` and copy the file tyr_sulfated.txt to tyr_phosphorylated.txt. Open that file for editing. Change these lines:
```
NAME sulfated
TYPES SULFATION
```
to this:
```
NAME phosphorylated
TYPES PHOSPHORYLATION

NOT VARIANT_TYPE PHOSPHORYLATION 
```
Now copy this new file tyr_phosphorylated.txt to ser_phosphorylated.txt. Open ser_phosphorylated.txt for editing. Change this line:
```
AA TYR
```
to this:
```
AA SER
```
Now, we have patch files for our centroid-level phosphorylated residues. All we have to do now is point Rosetta to these files. Open the file  ```<my_rosetta_directory>/rosetta_database/chemical/residue_type_sets/centroid/patches.txt```. Add these two lines to the bottom of the file:
```
patches/tyr_phosphorylated.txt
patches/ser_phosphorylated.txt
```
 Now run the docking_protocol application like this:
```
<my_rosetta_directory>/rosetta_source/bin/docking_protocol.linuxgccrelease -s rosetta_inputs/2lax_edited.pdb
```

This application should now work correctly, converting the input pdb coordinates to a centroid model, performing rigid-body docking, then converting back into a full-atom model ( which should now *correctly convert the phosphorylated residue types* ) before performing a final docking rotamer packing optimization (see the docking_protocol documentation for more information). 
# RNA Assembly

# Author

This README was written in Sep. 2011, by Rhiju Das (rhiju@stanford.edu); updated in Feb. 2012 after directory restructuring. Thanks to P. Kerpedjiev for suggestions.

# Brief Info
This demo illustrates a protocol to assemble models of large RNAs by first building their helical stems and inter-helical motifs, and then putting them together.

It is being published in a (primarily experimental) paper "A two-dimensional mutate-and-map strategy for non-coding RNA structure" by W. Kladwang, C. VanLang, P. Cordero, and R. Das (2011), Nature Chemistry.

# Running the demo

The example input files are in rosetta_input; you may wish to copy them locally with the command:

```
    cp rosetta_inputs/* .
```

Everything needed to run the job is created by the command:

```
    python scripts/setup_rna_assembly_jobs.py  add.fasta add_secstruct.txt 1y26_RNA.pdb  add_mutate_map_threetertiarycontacts.cst
```

The first two arguments are required -- the sequence_file and the secondary structure file [either in dot/bracket notation, or specifying Watson/Crick base pairs as pairs of numbers]. 

The last two arguments are optional; they supply the native pdb and any constraints, here derived from a high-throughput "mutate-and-map" strategy for RNA structure determination.

You may need to change the path to your rosetta executable ('EXE_DIR') in setup_rna_assembly_jobs.py [in which case you should get a warning!]. The script currently assumes that you are in the rosetta_demos/public/RNA_Assembly directory within the rosetta codebase.

Then run the Rosetta commands in :

```
    README_STEMS
    README_MOTIFS
    README_ASSEMBLE
```

You can see examples of these files and their output in example_output/. Please note that for these files I changed the 'nstruct' commands to create 100 models per motif. In reality you will want to make 2000-4000 MOTIF models, and then several thousand ASSEMBLE models. [You can just use one STEM model per helix, as that is supposed to be an ideal helix.] For some scripts to generate lots of models on a computer cluster, see note below.

The final 'outfile' is  add_assemble.out. We can extract models from it using:

```
    extract_pdbs<.exe> -in:file:silent add_assemble.out -tags S_000001 -database <path to your database> -out:file:residue_type_set rna
```

or using scripts like my `extract_lowscore_decoys.py`.

Caveat: The above protocol is a bit inflexible in that the motifs are modeled separately from each other -- if a loop/loop interaction occurs in the final global model it will not really be modeled correctly by the isolated loops. We are working on iterative methods to tackled this global de novo assembly question. For now the protocol seems to work well if there are experimental constraints that e.g., connect the loops.

# Appendix: Generating lots of models on a cluster

We often use either condor or LSF ('bsub') queuing systems on clusters, and have some handy scripts to set up these jobs, which are included in the Rosetta distribution.

The syntax (e.g., for README_MOTIFS) is:

```
    python ../../../rosetta_tools/rna/rosetta_submit.py  README_MOTIFS MOTIF_OUTPUT 50
```

this will create the directory MOTIF_OUTPUT with 50 subdirectories 0/, 1/, etc. for 50 parallel jobs. You can then run:

```
    source bsubMINI
```

or

```
    condor_submit condorMINI
```

to kick off jobs. The rosetta_submit.py script is pretty straightforward and should be easy to edit for any kind of cluster queuing system.

After the jobs have been running for a while, or have finished, you can type:

```
    python ../../../rosetta_tools/rna/easy_cat.py MOTIF_OUTPUT 
```

to concatenate your files.  [Do a similar thing with README_ASSEMBLE after getting the results of MOTIF].
Do the best picker option for the starter demo: BestFragmentsProtocol 

Submit the following sequence (in the FASTA format) to 
[[psipred|http://bioinf.cs.ucl.ac.uk/psipred]]:

    > 2JSV X
    MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE

After psipred finishes get the "download the machine learning scores results in 
plain text format" and give it a meaningful name, e.g. `2jsv.psipred.ss2`. This 
should have a header with the string "VFORMAT"

* Put this vall in the flags:

        [path]/rosetta/rosetta_database/sampling/filtered.vall.dat.2006-05-05.gz

  Update the flags with your correct database location

* Run the picker like this:

        [path]/rosetta/rosetta_source/bin/fragment_picker.linuxgccrelease @best-frags.flags 

  Output:
  * `frags.200.3mers`: for use with `-in:file:frag3`
  * `frags.200.9mers`: for use with `-in:file:frag9`
  * fsc files for fragment information.

* Now get the quality of these fragments:

        [path]/rosetta/rosetta_source/bin/r_frag_quality.linuxgccrelease -database [path]/rosetta/rosetta_database -in:file:native input_files/2jsvX.pdb -f output_files/frags.200.9mers

  Output:
  * frag_qual.dat

  Use gnuplot to visualize the quality data:

        plot "frag_qual.dat" u 2:4 w p

  This command simply makes a scatter plot with column 2 on the x-axis and 
  column 4 on the y-axis. 

More detailed information in:  
http://www.rosettacommons.org/manuals/archive/rosetta3.4_user_guide/dc/d10/app_fragment_picker.html

# Presenting author: Dominik Gront (dgront@chem.uw.edu.pl)
# Protocol Name: fragment picker : nnmake style
# Brief Description: The protocol substitutes nnmake
   
# Source code location:
Check out the mini SVN:
https://svn.rosettacommons.org/source/trunk/mini/

Fragment picker is located in:
https://svn.rosettacommons.org/source/trunk/mini/src/core/fragment/picking

Applications are in:
https://svn.rosettacommons.org/source/trunk/mini/src/apps/pilot/dgront/fragmentpicker

# running:
1) set up the path to minirosetta database
2) set up the path to vall database
3) run the picker:
picker.linuxgccrelease @quota-protocol.flags# Template Demo Directory

The entire workflow for this demo should be described in a file named README.md, formatted in gollum markdown.  
It should describe an entire work flow, with command lines, tested if possible.

The contents of each demo directory should be:

```
starting_files/
  -- directory in which the raw input files are given - these
     are provided for you and serve as the base for your
     tutorial
rosetta_inputs/
  -- directory in which the modified starting files should
     be placed which will be used as inputs to Rosetta.
     You may need to make modifications like stripping
     extra chains from the input PDB; store the modified
     PDB here and leave the unaltered one in starting_files 
scripts/
  -- python scripts, shell scripts, used in the workflow
  -- awk, grep, sed, etc. command lines

README.dox
  -- A prose or list description of how to perform the protocol

FOR_AUTHORS.txt
  -- A description for the demo creators of what their demo
     should achieve.
  -- Most of what starts in this file should end up in the
     README file as well.
```
This document briefly demos how to use favor_native_residue for Rosetta 3.4 in 
a fixed backbone context.  favor_native_residue is a mode (which debuted in 
Rosetta 2) that pushes the packer to prefer NOT making a mutation when the 
mutation is only weakly beneficial.  This is useful in design contexts where 
few mutations are desired.

The fixbb executable does not support favor_native_residue because the 
implementation is a little tricky.  (It could, but it currently does not).  The 
problem is that favor_native_residue uses constraints to do the energy 
favoring.  Constraints must be added to the Pose and ScoreFunction.  JD2 
doesn't instantiate the Pose until after the outer setup of fixbb has already 
completed, so there is no chance to load the constraints into fixbb.

To run the demo:

    /path/to/rosetta_scripts.default.<OS><compiler><mode> @options -database /path/to/rosetta_database/

The options file is annotated.

The XML file (inputs/favor_native_residue.xml) contains two movers.  
FavorNativeResidue, applied before the packing mover, causes the 
favor-native-residue behavior.  "bonus" is the energy bonus to assign to native 
residues.  1.5 will likely be overwhelming; 0.5 will be useful; 0.05 will be 
very weak.  You'll have to tune it to your application.
# Relax a Large Structure

The FastRelax protocol minimizes the input structure according to Rosetta's force field.
The attached XML file (demo.xml) is an example rosetta script that specifies how to relax a part of a large protein.
The script defines one mover (FastRelax) with a MoveMap, which specified the residues to relax.
The mover is referenced in the PROTOCOLS section as the single operation.

An empty MoveMap selects all residues for relaxation.
The attached script forces the positions 1-800 and 1000-1200 to be fixed and allows the rest of the protein to be moved.
On 3E0C.pdb for example (a protein that has 1011 positions), this will only relax positions between 800-1000 (watch out: the numbering is according to the pose and might not correspond to pdb numbering if there are missing numbers in the pdb).
The bb and chi parameters specify if the backbone and/or side chains should be relaxed.
The repeats parameter control the length of the relax simulation.

Multiple residue spans can be specified and will be parsed in the specified order.
For example, the snippet below turns off relax for all residues, and then specified that only residue 100-200 should be relaxed:

```
	<Span begin=1 end=1011 chi=0 bb=0/>
	<Span begin=100 end=200 chi=1 bb=1/>
```

If not explicitly specified, the following parameter values are used for FastRelax:

```
	- scorefxn: score12
	- repeats: 8
	- task_operations: InitializeFromCommandline, IncludeCurrent, and RestrictToRepacking
```

Command line arguments to run this script:

```
	./rosetta_scripts.linuxgccrelease -database ~/minirosetta_database/ -s starting_files/3E0C.pdb -parser:protocol demo.xml -ignore_unrecognized_res
```

The flag ignore_unrecognized_res asks Rosetta to ignore any residue types it doesn't recognize (HOH in this example).


Additional flags that can help:

```
	-linmem_ig <int>
	might reduce the memory usage of packing (althought it didn't seem to help on this example). 

	-relax::min_type lbfgs_armijo_nonmonotone
	This uses lbfgs minimizer instead of the bfgs minimizer that uses less memory. You can control this amount of this saving by specifying (default=64): -optimization::lbfgs_M. If the protein has 200 positions, lbfgs with lbfgs_M set to 200 will behave exactly like bfgs. You can trade-off memory usage for running time by trying to minimize for a longer time. Use (default=200) -optimization::lbfgs_max_cycles
```

This is a very simple design on a fixed backbone demo. If you have never run 
Rosetta before this is a good first demo to run, because it is very simple and 
has few options.

Use the files from the integration test, for example copying them to a new 
working directory:

    rosetta/rosetta_tests/integration/tests/fixbb

* Run like this:

        rosetta/rosetta_source/bin/fixbb.linuxgccrelease @flags_fullatom_dun10 -database ~/rosetta/rosetta_database/ > log.txt &

  The following files should be produced:

        1l2y_0001.pdb
        log.txt
        score.sc

* Open the structure and the input structure in pymol to observe sequence 
  changes from design.

* Systematically list sequence changes in the form of a sequence profile:

        ls 1l2y_0001.pdb > list.txt  # this would typically be many designed structures all in a list
        python rosetta/rosetta/tools/protein_tools/scripts/SequenceProfile.py -l list.txt -t 1l2y.pdb

* To control which residues are allowed at each sequence position you would add 
  a resfile (included in this demo) like so:

        rosetta/rosetta_source/bin/fixbb.linuxgccrelease @flags_fullatom_dun10 -database ~/rosetta/rosetta_database/ -resfile resfile.txt -out:suffix _resout > log_resout.txt &

  Open up the resfile.txt file to see its format. Briefly, NATRO leaves the 
  natural rotamer (and amino acid). NATAA leaves the amino acid at a position 
  but allows rotamer to change. ALLAA allows full design with any amino acid. 
  PIKAA followed by a list of single-letter-code amino acids restricts design 
  to just those amino acids.  So, for example:

        1 A PIKAA NT

  indicates that residue 1 can be either N or T.
# SWA Protein Loop Long

## Author
Rhiju Das, rhiju@stanford.edu

# Loop Remodeling by Enumeration: The Core Step of Protein 'Stepwise Assembly'

## Brief Description

This demo is an expansion of the demo swa_protein_main/. This demo involves a more detailed workflow for loops longer than 4-5 residues. It consists of buildup of the loop in 1-2 residue segments from either end, followed by chain closure. The workflow is described as a directed acyclic graph (DAG) of Rosetta jobs with well-defined dependencies.

All possible buildup-paths are followed through a dynamic-programming-like recursion, and up to a 1000 models are retained at each intermediate buildup path. The results are high in accuracy and low in energy, but requires significant computational expense (1000s of CPU-hours) and the ability to run a complex DAG on a cluster.

## Abstract

Consistently predicting protein structure at atomic resolution from sequence alone remains an unsolved problem in computational biophysics. Practical challenges involving protein loops arise frequently in ab initio modeling, comparative modeling, and protein design, but even these cases can become intractable as loop lengths exceed 10 residues and if surrounding side-chain conformations are erased. This demo illustrates a novel approach to protein modeling that is more powerful than prior methods that strive for atomic resolution. The central innovation is a ‘stepwise ansatz’ inspired by recent ab initio RNA algorithms, which resolves a conformational sampling bottleneck through residue-by-residue conformer enumeration and dynamic programming.


Reference: R. Das (2013) "Atomic-accuracy prediction of protein loop structures enabled by an RNA-inspired ansatz", under review.
More info: http://arxiv.org/abs/1208.2680

## Running

### Setup
You need to define an environment variable $ROSETTA with your Rosetta directory. Add to your .bashrc or .bash_profile a line like:

```
export ROSETTA='/Users/rhiju/src/rosetta/'   [change to your Rosetta directory]
```
 
You also need your system to know where the python scripts are for generating the DAG and running the jobs:

```
PATH=$PATH:$ROSETTA/rosetta_tools/SWA_protein_python/generate_dag/
PATH=$PATH:$ROSETTA/rosetta_tools/SWA_protein_python/run_dag_on_cluster/
```

### Example Python Command Line to Generate DAG
**Note: needs to be run in a copy of rosetta_inputs/**

This rebuilds residues 3-8 on the knottin scaffold 2it7 which has had the loop and all the protein's sidechains removed:

```
generate_swa_protein_dag.py  -loop_start_pdb noloop_2it7_stripsidechain.pdb  -native 2it7.pdb -fasta 2it7.fasta -cluster_radius 0.25 -final_number 1000   -denovo 1   -disulfide_file 2it7.disulf  -loop_res 3 4 5 6 7 8
```

If you want to do a quick run to test the overall workflow, you can keep the models within RMSD of 1.0 A to the experimental loop, use coarser backbone sampling, and save only 10 structures per run:

```
generate_swa_protein_dag.py  -loop_start_pdb noloop_2it7_stripsidechain.pdb  -native 2it7.pdb -fasta 2it7.fasta -cluster_radius 0.25 -final_number 1000   -denovo 1   -disulfide_file 2it7.disulf  -loop_res 3 4 5 6 7 8 -n_sample 9 -rmsd_screen 1.0 -nstruct 10
``` 

The outputs are:

```
  protein_build.dag
    text file outlining all the jobs, preprocessing and preprocessing script commands, and their dependencies. In condor DAGMAN format.

  CONDOR/
    directory with job definitions, in format similar to condor job definition format.

  REGION_3_2/, REGION_4_2/, ...
    directories that will hold outputs of each rosetta job. The numbers correspond to the N-terminal residue of the loop fragment reaching from the C-terminus endpoint of the loop, and the C-terminal residue of the loop fragment reaching from the N-terminal takeoff of the loop.  REGION_3_2 corresponds to models in which the loop fragments have met at the boundary between residues 2 and 3.
```

## How to Run the DAG

On condor clusters:

```
condor_submit_dag protein_build.dag
```

We have found condor_dagman can be slow due to latency in queuing rosetta jobs via Condor, unfortunately. A further problem is that there is currently no good universal solution to running DAGs on clusters, although packages like Pegasus and newer versions of Hadoop appear promising. 

For our own purposes, we have developed in-house python scripts (available in `rosetta_tools/SWA_protein_python/run_dag_on_cluster/`) to run the jobs by kicking off a master node that can queue jobs to slave nodes. Most recently, we have been using PBS/torque ('qsub') queueing, and you can run the jobs using 100 cores with the command:

```
 qsub README_QSUB 
```

which will queue from a master node:

```
 SWA_pseudo_dagman_continuous.py -j 100 protein_build.dag  > SWA_pseudo_dagman_continuous.out 2> SWA_pseudo_dagman_continuous.err
```

Further development for LSF clusters and MPI queuing is also under way. Please contact rhiju [at] stanford.edu with questions, or suggestions for supporting more general queuing systems.

**Note:** These scripts should work directly out of trunk for any version of Rosetta after March 2013.


AbInitio Structure Prediction Using Chemical-Shift Generated Fragments and NOE Distance Restraints
==================================================================================================

Written by Lei Shi.
Nikolas Sgourakis drafted the previous version.

---

We will use the chemical shifts to improve the fragments from which Rosetta builds up structures, and the NOEs to guide the Rosetta calculations towards the native structure.

Please see references at:
* rosetta abinitio: Bradley, P et al Science 2005
* chemical shift fragments: Shen Y et al. PNAS 2008;105:4685-4690
* chemical shift+NOE+RDC: Raman S, et al Science 2010

These Rosetta calculation steps are also described separately:
* Sgourakis NG et al JACS,2011,133(16):6288-98:

In this demo, we will use PDB 2JY7, which is a small protein (for demo purpose) and has experimental data deposited. Several scripts are provided in the scripts folder for formatting purposes:

	bmrb2talos.com
	cst_map_toCB.py
	upl2mini.csh
	scores.score.cfg

If you are from David Baker lab, there are scripts available to make setup easier without going through public servers. The following instructions should work just fine without having direct access to any Baker lab cluster.

Running the demo
----------------
1. Create following folders:  
    ```
    mkdir starting_inputs
    mkdir rosetta_inputs
    mkdir rosetta_inputs/talos_output
    mkdir rosetta_inputs/pick_cs_fragments
    ```

2. Download protein fasta and experimental data
Download fasta from http://www.pdb.org/pdb/explore/explore.do?structureId=2JY7  
    ```
    wget http://www.pdb.org/pdb/files/fasta.txt?structureIdList=2JY7 -O starting_inputs/t000_.fasta
    ```
Download chemical shift data from http://www.bmrb.wisc.edu/data_library/summary/index.php?bmrbId=15591  
    ```
    wget http://rest.bmrb.wisc.edu/bmrb/NMR-STAR2/15591 -O starting_inputs/raw.cs.bmrb
    ```
Download NOE data from http://restraintsgrid.bmrb.wisc.edu/NRG/MRGridServlet?pdb_id=2JY7&show_blocks=true&min_items=0:
    ```
    wget http://restraintsgrid.bmrb.wisc.edu/NRG/MRGridServlet?db_username=wattos1&format=ambi&mrblock_id=434910&pdb_id=2jy7&program=DYANA%2FDIANA&request_type=block&subtype=general+distance&type=distance
    echo "save file as starting_inputs/NOE_data.upl"
    ```

3. Format data for Rosetta use  
Formatting NOE: (Note only residues separated by more than 3 are kept in constraint)
The script `scripts/upl2mini.csh` only works with cyana format NOE:
    ```
    scripts/upl2mini.csh starting_inputs/NOE_data.upl > rosetta_inputs/NOE.cst
    scripts/cst_map_toCB.py rosetta_inputs/NOE.cst > rosetta_inputs/NOE.centroid.cst
    ```
Formmatting chemical shift data for TALOS:
    ```
    scripts/bmrb2talos.com starting_inputs/raw.cs.bmrb > rosetta_inputs/cs.talos
    ```

4. Generating talos predictions using http://spin.niddk.nih.gov/bax/nmrserver/talosn/ using rosetta_inputs/cs.talos
Save/copy pred.tab and predSS.tab to rosetta_inputs/talos_output

5. Generate fragment/profile from RobettaServer http://www.robetta.org/fragmentqueue.jsp using starting_inputs/t000_.fasta
Save/copy t000_.checkpoint to rosetta_inputs/

6. Pick fragments using secondary structure profile and chemical shift data:
```
Rosetta/main/source/bin/fragment_picker -database Rosetta/main/database/ -in::file::vall Rosetta//tools/fragment_tools/vall.apr24.2008.extended.gz -frags::n_frags 200 -frags::frag_sizes 3 9 -frags::sigmoid_cs_A 2 -frags::sigmoid_cs_B 4 -out::file::frag_prefix rosetta_inputs/pick_cs_fragments/frags.score -frags::describe_fragments rosetta_inputs/pick_cs_fragments/frags.fsc.score -frags::scoring::config scripts/scores.score.cfg -in:file:fasta starting_inputs/t000_.fasta -in:file:checkpoint rosetta_inputs/t000_.checkpoint -in:file:talos_cs rosetta_inputs/cs.talos -frags::ss_pred rosetta_inputs/talos_output/predSS.tab talos -in::file::talos_phi_psi rosetta_inputs/talos_output/pred.tab
```

7. Run Rosetta with the fragments made above and use NOEs to guide search
```
Rosetta/main/source/bin/minirosetta -database Rosetta/main/database/ -cst_fa_file rosetta_inputs/NOE.cst -cst_file rosetta_inputs/NOE.centroid.cst -abinitio:stage1_patch scripts/patch_atom_pair_constraint -abinitio:stage2_patch scripts/patch_atom_pair_constraint -abinitio:stage3a_patch scripts/patch_atom_pair_constraint -abinitio:stage3b_patch scripts/patch_atom_pair_constraint -abinitio:stage4_patch scripts/patch_atom_pair_constraint -score:patch scripts/patch_atom_pair_constraint -in:file:fasta starting_inputs/t000_.fasta -file:frag3 rosetta_inputs/pick_cs_fragments/frags.score.200.3mers -file:frag9 rosetta_inputs/pick_cs_fragments/frags.score.200.9mers -nstruct 1 -out:file:silent csrosetta_noe.out -run:protocol abrelax -abinitio::relax -overwrite
```
You can/should adjust the weights of NOE constraints in `scripts/patch_atom_pair_constraint`.
You should also change nstruct to generate desired number of models.
Larger is better depending on your available computer time, etc.
Note that the demo in abinitio_w_chemicalshift_only, you can add flags such as:
```
    -abinitio::rg_reweight 0.5
    -abinitio::rsd_wt_helix 0.5
    -abinitio::rsd_wt_loop 0.5
    -disable_co_filter true
    -abinitio::increase_cycles 10
```
to help sampling in centroid stage.
They are not used probably NOE constraint helps guided the search.

Processing the output
---------------------
1. Extract the low energy models:
    ```
    grep SCORE csrosetta_noe.out | sort –nk2 | head
    ```
The second column contains the energies of the lowest energy 10 models.
Select as the cutoff the energy on the last line.
You should also use NOE constraint energy as a criteria to select structures.
Example is only provided for total score.

2.  This command
    ```
    cull_silent.pl csrosetta_noe.out “score < cutoff”
    ```
will produce csrosetta.select.silent which contains the lowest energy 10 models.

3. Extract pdbs from selected silent file
    ```
    Rosetta/main/source/bin/extract_pdbs -database Rosetta/main/database/ -in::file::silent csrosetta.select.silent
    ```

4. Check convergence by superimposing the ten low energy models in pymol or your favorite molecular graphics.

5. Check convergence by clustering the lowest energy models (see clustering demo for instructions).

6. To see how NOE constraints are satisfied by a model:
    ```
    Rosetta/main/source/bin/r_cst_tool.linuxgccrelease -database Rosetta/main/database/ -in:file:s lowscore_1.pdb -cst_file rosetta_inputs/NOE.cst
    ```
`r_cst_tool` is a pilot program by Oliver Lange in `Rosetta/main/source/src/apps/pilot/olli/`
Model a Missing Loop
====================

Authors: Roland Pache, Michal Sperber, Steven Combs, George Rosenberger  
Last updated: August 2011 (RosettaCon9)

---

This demo shows how missing electron densities of several consecutive residues 
can be modeled using the loop modeling application (loopmodel) and the 
KInematic Closure algorithm (KIC).

The starting structure (1tr2_missing_density.pdb) is based on vinculin (1TR2). 
5 loop residues have been removed (32-36) and should be replaced by the 
sequence VDGKA for loop modeling (simulating missing electron density). 
Afterwards, this PDB structure can be used to model the loop. For this demo, 
the water molecules (HOH) have been removed and the structure was truncated to 
the first 132 residues.

Running the demo
----------------

1.  Insert the new residues into the structure file

    Open the file 1TR2_missing_density.pdb in the text editor of your choice. 
    Search the first gap line (residue 32). Search for the first residue Valine 
    in the file and copy all atoms to the new line. Repeat this step for all 
    other residues (DGKA) and insert the coordinates below Valine. The file 
    should then look like 1TR2_manually_added_dummy_residues.pdb.   Renumber 
    the residues you copied from another place to 32-36 and remove all 
    eventually inserted new lines. Save this file as 
    1TR2_manually_added_dummy_residues_renumbered.pdb.

2.  Create the loop file.

    Create a new file, called 1TR2.loop, and open it in your text editor. 
    Insert the following line:

        LOOP 31 37 37 0 1

    This excerpt from the loopmodel documentation describes the meaning of the 
    6 columns in that line:

        column1  "LOOP":     Literally the string LOOP, identifying this line as a loop
                             In the future loop specification files may take other data.
        column2  "integer":  Loop start residue number
        column3  "integer":  Loop end residue number
        column4  "integer":  Cut point residue number, >=startRes, <=endRes.
        column5  "float":    Skip rate. default - never skip (0)
        column6  "boolean":  Extend loop. Set to 1

    For this example, we select the one residue before and after the loop to 
    have real coordinates that can be used as anchor points by the KIC loop 
    modeling algorithm. The cut point residue number is set to the last loop 
    residue, since it must be inside the loop. The skip rate is set to 0 for 
    this short example (since we want to model this loop) and the extend loop 
    setting is set to true to idealize all bond lengths, bond angles and 
    torsion angles of the loop residues before modeling.

3.  Execution of the algorithm and definition of the flags

    Assuming that Rosetta 3.3 is installed and all paths are set correctly, 
    open your shell and change the directory to the one where the demo files 
    are stored.

        loopmodel.linuxgccrelease -database /pathtoyourdb/ -loops:input_pdb 1TR2_manually_added_dummy_residues_renumbered.pdb -loops:loop_file 1TR2.loop -loops:remodel perturb_kic -loops:refine refine_kic -ex1 -ex2 -nstruct 1 -loops:max_kic_build_attempts 100 -in:file:fullatom

    Brief descriptions for all the components of this command-line:

        loopmodel.linuxgccrelease : loopmodel application (linuxgccrealease or macosgccrelease)
        -database : path to your Rosetta 3.3 DB
        -loops:input_pdb 1TR2_manually_added_dummy_residues_renumbered.pdb : name of your edited pdb
        -loops:loop_file 1TR2.loop : name of your loops file
        -loops:remodel perturb_kic : kinematic closure based loop modeling low resultion stage (side chains: centroids)
        -loops:refine refine_kic : kinematic closure based loop modeling high resultion stage (side chains: fullatom)
        -ex1 : extra chi rotamers for chi-1 angle (+/- 1 stddev from the optimal rotamer for better loop reconstruction)
        -ex2 : extra chi rotamers for chi-2 angle (+/- 1 stddev from the optimal rotamer for better loop reconstruction)
        -nstruct 1 : number of structures to generate (set to at least 1000 for real application; computationally expensive)
        -loops:max_kic_build_attempts 100 : the maximal number of trials the algorithm should do to find a closed confirmation for the loop (default: 100); can be increased for difficult problems.
        -in:file:fullatom : keep native amino acid side chain confirmations of the non-loop residues (residues within 10 Angstrom of the loop will be remodeled by default)

4.  Analysis of the results

    If you set -nstruct > 1, look at the Rosetta energy score in the standard 
    output and identify the model with the lowest energy. Compare them visually 
    using your favorite molecular visualization application (pymol, etc).
    Sample output files (incl. visualization using pymol) can be found in the 
    output directory.
    
Antibody Docking
================

The entire workflow for this demo should be described in this file.
It should describe an entire work flow, with command lines, tested if possible.

Authors:
* Jianqing Xu (xubest at gmail dot com)
* Christine Tinberg
* Jeff Gray
* Angela Loihl

Demo files
----------

`starting_files/`
* directory in which the raw input files are given - these are provided for you and serve as the base for your tutorial

`rosetta_inputs/`
* empty in this demo
* directory in which the modified starting files should be placed which will be used as inputs to Rosetta.  You may need to make modifications like stripping extra chains from the input PDB; store the modified PDB here and leave the unaltered one in starting_files 

`scripts/`
* empty in this demo
* python scripts, shell scripts, used in the workflow
* awk, grep, sed, etc. command lines

`README.md`
* A prose or list description of how to perform the protocol

`FOR_AUTHORS.txt`
* A description for the demo creators of what their demo should achieve.
* Most of what starts in this file should end up in the README file as well.

Running the demo
----------------

1.  The first thing a user should notice is that there’s an antibody 
    modeler in Rosetta 3, but still under development.  As of August 
    2011, Rosetta3 should be used for camelid antibody modeling, but 
    Rosetta2 should be used for other antibody modeling and for antibody 
    docking via SnugDock.  The current stable version of Rosetta2 
    (Rosetta++) is the released Rosetta-2.3.1

2.  Obviously, you need structures of both antibody and antigen in order 
    to do antibody-antigen docking. If you don't have antibody structures, 
    but have antibody sequences, you can use Gray lab antibody homology 
    modeling server (http://antibody.graylab.jhu.edu/) and input the sequence 
    of the light and heavy chain. You will get best 10 structures.

    If you want to manually run the scripts yourself, you can download 
    the scripts source code, example, and instructions from the link below: 
    https://svn.rosettacommons.org/source/trunk/antibody/.
    Again, please realize that the H3 loop modeling is still from Rosetta++.

3.  Download the released version of Rosetta++: 
    https://svn.rosettacommons.org/source/branches/releases/rosetta-2.3.1/ 

    For people outside of the rosetta community: 
    https://www.rosettacommons.org/software/academic/2.3.1/RosettaSnugDock-2.3.1.tgz

4.  Compile rosetta:
    ```
    tar –zxvf RosettaSnugDock-2.3.1.tgz (if you download the second link)
    cd rosetta++
    scons mode=release –j12    (assuming you can use 12 CPUs)
    ```

5.  The SnugDock example in Rosetta++ can be found at:
    https://svn.rosettacommons.org/source/branches/releases/rosetta-2.3.1/example/
    * Besides the original example shown above, we made a new example in the current directory.
      Please be careful with different flags used in the command line.
      The documentations of SnugDock options can be found at:
      http://www.rosettacommons.org/guide/SnugDock.
      Please also be careful with each "paths.txt" file.
    * Do Ensemble Prepack:
      ```
      cd ./PrePack_input # you need AB_model*.pdb, ABRM.fab, ABRM.pdb, ABRM.unbound.pdb, Antigen.pdb, pdblist1, pdblist2
      cd ../Prepack # you need `paths.txt` and `prepack.bash` file available
      ./prepack.bash
      ```
      After prepack, you will see `*.ppk` files in the PrePack_input directory.
      Your original pdblist1 and pdblist2 files will be modified as well, please see the README file inside that folder.
    * Do SnugDock+Ensemble:
      ```
      cd ../SnugDock # you need EnsembleDock_plus_SnugDock.bash and paths.txt file
      ./EnsembleDock_plus_SnugDock.bash # please realize that the example we used here is a little slow, due to the protein size
      ```

6.  Some extra information you may need, besides the example tutorial linked above:

    * Make fab file:
      The CDR loops of antibody should point to the antigen.
      By specifying the antibody loops in the fab file, one can reduce the computational cost for global docking.
      The scripts to make fab file can be found at:  
      https://svn.rosettacommons.org/source/branches/releases/rosetta-2.3.0/rosetta_scripts/docking/  
      Run makefab.pl on your pdb of choice.
      ```
      ./makefab.pl `input pdb` <heavy and/or light chain i.e. HL>
      ./makefab.pl AB_model1.pdb HL
      ```

    * Make `FR02.pdb` complex file:
      Use pymol to open both the antibody and antigen in one session and save both into one pdb file.
      In the example: `ABRM.pdb`.
      It's better to point the antibody CDRs to the antigen, and keep them at a certain distance.

1. *.ppk files are the resutls of prepacking
2. after prepack, the "pdblist1" and "pdblist2" files will be modified

    for example, in this case, the orignal "pdblist1" file is just:
	AB_model1.pdb
	AB_model2.pdb
	AB_model3.pdb
	AB_model4.pdb
	AB_model5.pdb
	AB_model6.pdb
	AB_model7.pdb
	AB_model8.pdb
	AB_model9.pdb
	AB_model10.pdb
 
β-Strand Homodimer Design
=========================

This outlines how to use the applications involved in finding exposed beta-strands and then designing a protein with an exposed beta strand to be a homodimer. 
Written by Ben Stranges (stranges at unc dot edu)

There are three applications associated with this demo:

    rosetta/rosetta_source/src/apps/public/scenarios/beta_strand_homodimer_design/homodimer_design.cc
    rosetta/rosetta_source/src/apps/public/scenarios/beta_strand_homodimer_design/homodimer_maker.cc
    rosetta/rosetta_source/src/apps/public/scenarios/beta_strand_homodimer_design/exposed_strand_finder.cc

It would also be a good idea to look at the doxygen for these apps at:

    rosetta/rosetta_source/doc/apps/public/scenarios/beta_strand_homodimer_design.dox

It explains how everything works. 

The idea is that you scan through a list of pdbs and find ones with exposed beta strands, you then make potential homodimers along these strands then design the residues at the interface to stabilize the homodimer. 

Running the Demo
----------------

1.  Find exposed beta-strands:
    ```
    rosetta/rosetta_source/bin/exposed_strand_finder.linuxgccrelease -s 2a7b_mpm.pdb.gz -database rosetta/rosetta_database @finder_options > exposed_strands.txt
    ```

    The contents of `finder_options` are:

    * general options  
        `-ignore_unrecognized_res true`  
        `-packing::pack_missing_sidechains`  
        `-out::nooutput`: this protocol manages its own output, prevent job distributor from helping  
        `-mute core basic protocols.jd2.PDBJobInputter`

    * app specific options  
        `-beta_length 5`: how long of an exposed strand do you look for  
        `-sat_allow 2`: how many satisfied bb atoms do you allow in this range

    * allow alignment of a found strand to some target protein  
        `-check_rmsd false`: setting this to false prevents the code from doing rmsd comparisons

    * `-native anti_model.pdb`

    * `-strand_span B 5 11`

    In reality you'll probably want to use `-l` instead of `-s` to pass a bigger list of pdbs to look for exposed strands in.
    This is just an example that will work.
    When this runs look at the output in exposed_strands.txt, the important line is this one:

        ExposedStrand: FILE:         2a7b_mpm   CHAIN:    A   START:   806   END:   812   H_BONDS:   0

    This means that there is an exposed strand in pdb 2a7b_mpm in chain A between residues 806 and 812 and there are 0 satisfied bb_bb H bonds in every other residue along this span. 
    This information is then used in the next step.

    This application has another mode that allows you to match a found exposed strand onto a beta-strand involved in the interaction with another protein. It is still in development and not really known to work. Use at your own risk. To activate it you need to set -check_rmsd to true and pass a structure with -native and use the -strand_span option.

2.  Make the potential homodimers:  
    There are two ways to do this.
    If you are only doing it for one structure is is easy just to use this command line:

    ```
    rosetta/rosetta_source/bin/homodimer_maker.linuxgccrelease -s 2a7b_mpm.pdb.gz -database rosetta/rosetta_database @finder_options > maker_tracers
    ```
    ```
    -run::chain A
    -sheet_start 806
    -sheet_stop 812
    -window_size 5 
    -ignore_unrecognized_res true
    -mute core protocols.moves.RigidBodyMover basic.io.database
    ```

    If you need to run a bunch of these from the output of the exposed strand finder I have provided a script that reads a file (runner_input).
    This file should be multiple lines instead of just the one here.
    It's structure is
    ```
    /path/to/pdb/1av3.pdb chainletter betastart betaend
    ```
    You will need to modify some of the paths in runner.sh. Then run this command:
    ```
    ./runner.sh runner_input
    ```

    This outputs a bunch of pdbs:
    ```
    2a7b_mpm_A806_anti_wind_1_step_-1.pdb
    2a7b_mpm_A808_parl_wind_2_step_1.pdb
    2a7b_mpm_A808_parl_wind_2_step_0.pdb
    2a7b_mpm_A808_parl_wind_2_step_-1.pdb
    2a7b_mpm_A806_parl_wind_1_step_1.pdb
    2a7b_mpm_A806_parl_wind_1_step_0.pdb
    2a7b_mpm_A806_parl_wind_1_step_-1.pdb
    2a7b_mpm_A806_anti_wind_1_step_1.pdb
    ```

    However for these purposes we are only interested in `2a7b_mpm_A806_anti_wind_1_step_1.pdb` so I removed the rest in the interest of saving space.

3.  Homodimer design:  
    The next step is to take the output from above and make the files you need for symmetry.
    To do this you will need to use the symmetry script as so:
    ```
    perl rosetta/rosetta_source/src/apps/public/symmetry/make_symmdef_file.pl -m NCS -a A -i B -p 2a7b_mpm_A806_anti_wind_1_step_1.pdb > symmdef
    ```

    This makes a bunch of files:
    ```
    2a7b_mpm_A806_anti_wind_1_step_1_symm.pdb
    2a7b_mpm_A806_anti_wind_1_step_1_model_AB.pdb
    2a7b_mpm_A806_anti_wind_1_step_1_INPUT.pdb
    2a7b_mpm_A806_anti_wind_1_step_1.kin
    symmdef
    ```

    However, the only one you really need is `2a7b_mpm_A806_anti_wind_1_step_1_INPUT.pdb` and `symmdef` so I removed the rest.
    Now you are ready for the full design runs. I suggest using mpi compiled executables but the command line below is general. 
    Run this command:
    ```
    rosetta/rosetta_source/bin/homodimer_design.linuxgccrelease -s 2a7b_mpm_A806_anti_wind_1_step_1_INPUT.pdb.gz -database rosetta/rosetta_database -symmetry:symmetry_definition symmdef @design_options
    ```
    See the comments in the design_options file for descriptions of what does what.

    This will output designed structures with the name: `2a7b_mpm_A806_anti_wind_1_step_1_INPUT_000x.pdb.gz` and a score file: `score.fasc`.
    From there it is up to you to to chose how you will determine which designs meet your needs.

    Then you should probably run the output through the InterfaceAnalyzer. See documentation for it here:

    rosetta_source/doc/apps/public/analysis/interface_analyzer.dox
Design the Rac/Raf interface
============================

This demo will walk through the steps of designing a protein-protein interface.
The goal of this protocol is to predict mutations on Raf that will enable it to 
bind Rac.

Demo files
----------

    README.md
    starting_files/
      -- 1c1y.pdb.gz
      -- 2ov2.pdb.gz
    rosetta_inputs/
      -- design_script.xml (dock/design script)
      -- raf-rac.pdb (raf-rac starting interface)
    scripts/
      -- score_vs_rmsd.R (creates score vs rmsd plot)
    output_files/
      -- directory to store output files

Running the Demo
----------------

1. Create model of Raf-Rac interaction.
    * Uncompress pdb files:
      ```
      gunzip 1c1y.pdb.gz
      gunzip 2ov2.pdb.gz
      ```
    * Open 1c1y and 2ov2 with Pymol.
    * Superimpose chain A of 1c1y with chain A of 2ov2.  Pymol command:
      ```
      super 2ov2 and chain A, 1c1y and chain A
      ```
    * Select modeled Raf-Rac complex. Pymol command:
      ```
      select raf-rac, (1c1y and chain B)+(2ov2 and chain A)
      ```
    * Save Molecule raf-rac.  In pymol:
      ```
      File->Save Molecule->raf-rac.pdb
      ```

2. Run dock/design protocol using rosetta_scripts.
    * Change working directory to the output directory.
      ```
      cd output_files
      ```
    * Run design_script.xml using rosetta_scripts executable.
      ```
      ~/mini/bin/rosetta_scripts.macosgccrelease -s ../rosetta_inputs/raf-rac.pdb -database ~/minirosetta_database/ -parser:protocol ../rosetta_inputs/design_script.xml -in:file:native ../rosetta_inputs/raf-rac.pdb -ex1 -ex2 -ignore_unrecognized_res -nstruct 1 -overwrite
      ```
    * Options:
      * `-ex1 -ex2`: expand rotamer library for chi1 and chi2 angles used in repacking/design
      * `-ignore_unrecognized_res`: ignores HETATM lines in input PDB file
      * `-nstruct 1`: specifies how many times the protocol is run (one decoy is output for each run)
      * `-overwrite`: overwrites decoys from previous runs

    Each run should take about 90 seconds to complete.
    For production runs, nstruct should be set to 1000 or greater.
    This protocol returns decoys named `raf-rac_####.pdb` and a file named 
    `score.sc` that contains scores for each decoy.

3. Postprocessing output from dock/design run.
    * To see the sequence changes in chain B due to design grep out chain B: 
      grep " B " raf-rac.pdb > chainBin.pdb then use 
            rosetta_tools/protein_tools/scripts/SequenceProfile.py
      to diff the sequences

    * Run R script to generate Interface Score vs. Interface RMSD plot. R is 
      easy to install from [[here|http://cran.r-project.org/bin/macosx/]] as a 
      pkg that self installs (do "which R" to check install first). The above 
      link is fro MacOSX, but R is available for many other operating systems 
      as well.

            CMD BATCH scripts/score_vs_rmsd.R

      Output: `score_vs_rmsd.pdf`
# Optimize Ligand Hydroxyl Hydrogens

This tutorial assumes a unix-style command line (or cygwin on windows).

## Build your ligand
In general the ligand atoms are marked by "HETATM", but check the ligand with pymol after you've grepped out the HETATM lines:
```
grep HETATM starting_inputs/cel5A_glucan.pdb > starting_inputs/cel5A_lig_noH.pdb
```
For this step you need to add on hydrogens. We use avogadro because its open source, but you can choose other software to place hydrogens:
http://avogadro.openmolecules.net
(this example is with 1.0.3 on mac)

Open in avogadro, choose Build--> add hydrogens, then save the molecule in MDL SDfile format cel5A_lig.mol

Alternatively, you can save the molecule in PDB format and convert the PDB file into a mol file using babel (http://openbabel.org/) as follows:
```
babel -ipdb starting_inputs/cel5A_lig.pdb -omol > starting_inputs/cel5A_lig.mol
```
Or the third alternative is to open the PDB file with pymol and save the molecule with a .mol extension to force pymol to save as a mol format. 

## Make a params file for the ligand
The next step is to create a params file for rosetta. Params file contains the internal coordinates of atoms, connectivity, charge of each atom, rosetta atom type. Most importantly for this demo it contains PROTON_CHI lines which specify the proton atoms that rosetta will simple around. 
-n specifies the name of the ligand in rosetta
```
python src/python/apps/public/molfile_to_params.py starting_inputs/cel5A_lig.mol -n cel
```
The output params file is called cel.params. 
The code will sample proton chi's at the explicit values stated in the params file, and then perform a minimization on those chi's.
You shouldn't have to change the params, but if you want to add sampling explicitly you can add more angles to sample. 
A sample proton chi is:
```
CHI 1  C1   C2   O1   H8
PROTON_CHI 1 SAMPLES 3 60 -60 180 EXTRA 0
```

## Setup files for rosetta
Pull out the protein without ligand:
```
grep ATOM starting_inputs/cel5A_glucan.pdb > rosetta_inputs/cel5A_input.pdb
```
Add back in the ligand pdb from molfile to params:
```
cat cel_0001.pdb >> rosetta_inputs/cel5A_input.pdb 
```

## Run rosetta enzdes
Now we can run the enzdes app in rosetta with minimal flags.
This optimizes proton chis on the ligand while also repacking sidechains:
```
path/to/EnzdesFixBB.[platform][compiler][mode] -s rosetta_inputs/cel5A_input.pdb -extra_res_fa cel.params -database path/to/minirosetta_database/ -out:file:o cel5A_score.out -nstruct 1 -detect_design_interface -cut1 0.0 -cut2 0.0 -cut3 10.0 -cut4 12.0 -minimize_ligand true
```
This optimizes proton chis on the ligand without repacking sidechains:
```
path/to/EnzdesFixBB.[platform][compiler][mode] -s rosetta_inputs/cel5A_input.pdb -extra_res_fa cel.params -database path/to/minirosetta_database/ -out:file:o cel5A_score.out -nstruct 1 -detect_design_interface -cut1 0.0 -cut2 0.0 -cut3 0 -cut4 0 -minimize_ligand true
```
Both runs should produce the PDB file cel5A_input__DE_1.pdb and the score file cel5A_score.out, which can be placed into the output directory under different names, e.g. cel5A_output_nopack.pdb or cel5A_output_w_repack.pdb

Flag descriptions:
```
-extra_res_fa specifies the params file for new residues types (the glucan in this case)
-out:file:o specifies the file name for the enzdes-style score output. This contains extra information about the output design, like packing, interface energy, and many more
-nstruct 1 specifies one run and one output pdb; the packing is stochastic so for more sampling use a higher nstruct. 10-100 is recommended for most ligands.
-detect_design_interface tells rosetta to set up the designable and packable residues (in a packer task) based on distance from the ligand. Distances are calculated from every ligand heavy atom to the CA of amino acids.
[default values in brackets]
-cut1: CA less than cut1 is designed [6]
-cut2: CA between cut1 and cut2, with CA --> CB vector pointing towards ligand, is designed [8]
-cut3: CA less than cut3 is re-packed [10]
-cut4: CA between cut3 and cut4 with CA --> CB vector pointing towards ligand, is re-packed [12]
-minimize_ligand true  allow ligand torsions to minimize
```

## More Complete Energy Function / Sampling
If desired, use a more complete energy function and more sampling as in this example flags file
Recommended full flags for a more careful run:
```
enzdes_flags
```

## Full Enzyme Design (?)
This setup can also be used for full enzyme design
This run is very close to a full design of the active site. For a full design just change cut1 and cut2, e.g.
```
-cut1 6 -cut2 8
```
ERRASER Demo
============


This demo illustrates the ERRASER (Enumerative Real-Space Refinement ASsitted 
by Electron density under Rosetta) protocol, which improves an RNA 
crystallographic model using Rosetta under the constraint of experimental 
electron density map. It was written in Mar. 2012, by Fang-Chieh Chou (fcchou 
at stanford dot edu) and based on a paper to be published:

* Chou, F.C., Sripakdeevong, P., Dibrov, S.M., Hermann, T., and Das, R. Correcting pervasive errors in RNA crystallography with Rosetta, arXiv:1110.0276. 

A preprint is available at:

* http://arxiv.org/abs/1110.0276

Setting up the demo
-------------------

The example input files are in rosetta_input; you may wish to copy them locally 
with the command:

    cp rosetta_inputs/* ./

Python codes needed to run the job are located at 
rosetta/rosetta_tools/ERRASER/

The following setup steps are required prior to running this demo.

1. Download and install PHENIX from http://www.phenix-online.org/. PHENIX is 
   free for academic users.

2. Ensure you have correctly setup PHENIX. As a check, run the following 
   command:
   
       phenix.rna_validate 

3. Check if you have the latest python (v2.7) installed. If not, go to the 
   rosetta/rosetta_tools/ERRASER/ folder and run 

        ./convert_to_phenix.python

    This will change the default python used by the code to phenix-built-in 
    python, instead of using system python.

4. Set up the environmental variable "$ROSETTA", point it to the Rosetta 
   folder. If you use bash, append the following lines to ~/.bashrc:

        ROSETTA=YOUR_ROSETTA_PATH; export ROSETTA" # Change YOUR_ROSETTA_PATH to the path in your machine!

    Also add the ERRASER script folder to $PATH. Here is a bash example:

        PATH=$PATH:YOUR_ROSETTA_PATH/rosetta_tools/ERRASER/" # Change YOUR_ROSETTA_PATH to the path in your machine!

Now you are ready to go!

Running the demo
----------------

As a fast quick run, run the following command:

    erraser.py -pdb 1U8D_cut.pdb -map 1U8D_cell.ccp4 -map_reso 1.95 -fixed_res A33-37 A61 A65

In this example, we specify the input pdb file, ccp4 map file and the map 
resolution. The -pdb and -map are required option, and -map_reso is optional but 
recommended (default is 2.0 if no input is given). In this example, the 
"1U8D_cut.pdb" file is a segment cutting from a deposited PDB file. The input 
pdb file should follow the standard PDB format, and no pre-processing is 
needed. The input map must be a CCP4 2mFo-DFc map. To avoid overfitting, Rfree 
reflection should be removed during the creation of the map file.

Note that we also manually fixed the position of residue 33-37, 61 and 65 in 
chain A, therefore we will only optimize residue 62-64. The -fixed_res argument 
is optional.

After the job finished successfully, you should see the output pdb file 
"1U8D_cut_erraser.pdb" in the current folder. A sample output is in the 
example_output folder for comparsion.

Note that the output file is in the standard PDB format and inherits all the 
ligands, metals and waters from the input pdb file (these atoms are not 
optimized in ERRASER). The user can then refine the output model using PHENIX 
or other refinement packages without any post-processing.

By inspecting the structure, you should be able to see a backbone conformation 
change at residue 63-64 after ERRASER.

Arguments for erraser.py (for your reference)
---------------------------------------------

#### Required:

* -pdb  
  Format: -pdb \<input pdb>  

  The starting structure in standard pdb format

* -map  
  Format: -map <map file>  

  2mFo-DFc map file in CCP4 format. Rfree should be excluded.

#### Commonly used:

* -map_reso  
  Format: -map_reso <float>  
  Default: 2.0  

  The resolution of the input density map. It is highly recommended to input 
  the map resolution whenever possible for better result.

* -out_pdb  
  Format: -out_pdb <string>  
  Default: \<input pdb name>\_erraser.pdb.  

  The user can output to other name using this option.

* -n_iterate  
  Format: -n_iterate <int>  
  Default: 1  

  The number of rebuild-minimization iteration in ERRASER. The user can 
  increase the number to achieve best performance. Usually 2-3 rounds will be 
  enough. Alternatively, the user can also take a ERRASER-refined model as the 
  input for a next ERRASER run to achieve mannual iteration.

* -fixed_res  
  Format: -fixed_res <list>  
  Default: <empty>  
  Example: A1 A14-19 B9 B10-13  (chain ID followed by residue numbers)  

  This allows users ton fix selected RNA residues during ERRASER. For example, 
  because protein and ligands are not modeled in ERRASER, we recommand to fix 
  RNA residues that interacts strongly with these unmodeled atoms. ERRASER will 
  automatically detect residues covalently bonded to removed atoms and hold 
  them fixed during the rebuild, but users need to specify residues having 
  non-covalent interaction with removed atoms mannually.

* -kept_temp_folder  
  Format: -kept_temp_folder <True/False>  
  Default: False  

  Enable this option allows user to examine intermediate output files storing 
  in the temp folder. The default is to remove the temp folder after job 
  completion.

#### Other:

* -rebuild_extra_res  
  Format/Default: Same as -fixed_res  

  This allows users to specify extra residues and force ERRASER to rebuild 
  them. ERRASER will automatically pick out incorrect residues, but the user 
  may be able to find some particular residues that was not fixed after one 
  ERRASER run. The user can then re-run ERRASER with -rebuild_extra_res 
  argument, and force ERRASER to remodel these residues.

* -cutpoint_open  
  Format/Default: Same as -fixed_res  

  This allows users to specify cutpoints (where the nucleotide next to it is 
  not connected to itself) in the starting model. Since ERRASER will detect 
  cutpoints in the model automatically, the users usually do not need to 
  specify this option.

* -use_existing_temp_folder  
  Format: -use_existing_temp_folder <True/False>  
  Default: True  

  When is True, ERRASER will use any previous data stored in the existing temp 
  folder and skip steps that has been done. Useful when the job stopped 
  abnormally and the user try to re-run the same job. Disable it for a fresh 
  run without using previously computed data.

* -rebuild_all  
  Format: -rebuild_all <True/False>  
  Default: False  

  When is True, ERRASER will rebuild all the residues instead of just 
  rebuilding errorenous ones. Residues in "-fixed_res" (see below) are still 
  kept fixed during rebuilding. It is more time consuming but not necessary 
  leads to better result. Standard rebuilding with more iteration cycles is 
  usually prefered.

* -native_screen_RMSD  
  Format: -native_screen_RMSD <float>  
  Default: 2.0  

  In ERRASER default rebuilding, we only samples conformations that are within 
  2.0 A to the starting model (which is the "native" here). The user can modify 
  the RMSD cutoff. If the value of native_screen_RMSD is larger than 10.0, the 
  RMSD screening will be turned off.
# RECCES 
## (Reweighting of Energy-function Collection with Conformational Ensemble Sampling)

## Author
This README was written in May 2015, by Fang-Chieh Chou (fcchou@stanford.edu).

## Brief
This demo illustrates the RECCES pipeline for computing the folding free energy of RNA helices. To keep it simple, we just show how to run RECCES simulated tempering (ST) on one construct.

Python codes needed to run the job are located at tools/recces/. You need to include this path into your PYTHONPATH to run the following demo. The python codes have only been tested with Python v2.7.

## Detail
We first do a quick a prerun. First create a separate folder:
```
    mkdir prerun
    cd prerun/
```
Now run the following commands:
```
    recces_turner -score:weights stepwise/rna/turner -n_cycle 300000 -seq1 gu -seq2 ac -temps 0.8 -out_prefix prerun
    recces_turner -score:weights stepwise/rna/turner -n_cycle 300000 -seq1 gu -seq2 ac -temps 1.0 -out_prefix prerun
    recces_turner -score:weights stepwise/rna/turner -n_cycle 300000 -seq1 gu -seq2 ac -temps 1.4 -out_prefix prerun
    recces_turner -score:weights stepwise/rna/turner -n_cycle 300000 -seq1 gu -seq2 ac -temps 1.8 -out_prefix prerun
    recces_turner -score:weights stepwise/rna/turner -n_cycle 300000 -seq1 gu -seq2 ac -temps 3.0 -out_prefix prerun
    recces_turner -score:weights stepwise/rna/turner -n_cycle 300000 -seq1 gu -seq2 ac -temps 7.0 -out_prefix prerun
    recces_turner -score:weights stepwise/rna/turner -n_cycle 300000 -seq1 gu -seq2 ac -temps 30  -out_prefix prerun
```
These preruns generate data for computing the ST weights in the following run. They may each take about 1 hour.

The ST weights can then be determined using the following code snippet:
```
    from recces.util import weight_evaluate
    weight_evaluate('./', 'prerun_hist_scores.gz')
```
Here it outputs a list of temperatures and the corresponding ST weights.

Now we create a new folder `ST` and run simulated tempering:
```
    mkdir ST
    recces_turner -score:weights stepwise/rna/turner -seq1 gu -seq2 ac -n_cycle 9000000 -temps 0.8 1 1.4 1.8 3 7 30 -st_weights 0 7.33 14.6 17.32 18.87 18.34 17.09 -out_prefix ST -save_score_terms
```
We also need to run a simulation at infinite temperature:
```
    recces_turner -score:weights stepwise/rna/turner -seq1 gu -seq2 ac -n_cycle 300000 -temps -1 -out_prefix kT_inf -save_score_terms
```
Note that here we use the "-save_score_terms" option to cache the contributions of each score term, so we may easily reweight the score function.

After the end of the run, the free energy can be computed using python codes:
```
    from recces.data import SingleSimulation, KT_IN_KCAL, N_SCORE_TERMS
    curr_wt = [0.73, 0.1, 0.0071, 0, 4.26, 2.46, 0.25, 0, 1.54, 4.54]
    sim = SingleSimulation('ST/', curr_wt)
    print sim.value, sim.value * KT_IN_KCAL
```
Here sim.value gives the free energy of the molecule in the unit of kT (T = 37C). To convert it to kcal/mol, we multiply the number by KT_IN_KCAL. In the `example_output` directory, the values in kT and in kcal/mol were: `1.35766728357 0.836772916126`. In an independent rerun, we got: `1.25722532755 0.774867389305`.

We may also reweight the score function and obtain the new free energy:
```
    import numpy as np
    sim.reweight(np.ones(N_SCORE_TERMS))
    print sim.value, sim.value * KT_IN_KCAL
```
In the `example_output` directory, we get: `31.8988957998 19.660289662`. In an independent rerun, we got: `31.9550773046 19.694916085`.

To determine the nearest neighbor energies, you need to add/subtract several free energies based on simulations of several single-stranded and double-stranded constructs.Homology Modeling with End Extension
====================================

Created: Aug 5, 2011  
Last Updated: Aug 5, 2011

This demo will illustrate a protocol for modeling the structure of a protein 
given its sequence, provided:

- There are sequence homologs in the PDB.
- But an n-terminal region of around 140 residues has no identifiable sequence 
  similarity to a protein structure in the pdb.

The input to the protocol is a sequence of interest.  The output will be a set 
of Rosetta-generated structural models.

Outline
-------

We will break down homology_modeling_with_end_extension into 2 parts:

1. A protocol for homology modeling of the protein (the region of the protein 
   that has a homolog in the pdb).

2. A protocol for *ab initio* modeling of the n-terminal region while keeping 
   the previously modeled homologous region in part 1 rigid.

Part 1: Homology modeling
--------------------------

There are 2 ways to do the homology modeling:

1.  Use the rosetta_comparative_modeling_protocol

2.  Use the [[Robetta server|http://www.robetta.org]].  Robetta can 
    automatically generate models for the region that has a homolog in the pdb. 
    If using Robetta, use its output (the generated models for the homologous 
    region) and skip to Part 2 (Broker Rigid Chunk Claimer).

The following 7 steps are for the rosetta_comparative_modeling_protocol option, 
based on a protocol developed by James Thompson and TJ Brunette.

1.  Get homologs and alignments:

    Prior to homology modeling, PBD homologs for the protein of interest need 
    to be identified and alignments generated.  Several external programs can 
    be used for this step like HHpred, psi-blast, or hhsearch.

    In this demo, we will use HHpred. To submit an HHpred job go to 
    http://toolkit.tuebingen.mpg.de/hhpred, paste the sequence of your protein 
    of interest into the text field, and then click the "Submit job" button. 
    Getting results should take a few minutes. On the results page, you can 
    save the alignments to a file by clicking the Save link under the "Results" 
    tab. We include an example output file below:

        ./starting_files/query.hhr

    At this point, you may truncate the query so that the long (not well 
    aligned) n-terminal region is removed. Then resubmit the truncated fasta to 
    HHpred for more accurate alignments and continue with the comparative 
    modeling protocol using the truncated sequence. For this example, we do not 
    truncate the sequence for simplicity.  We will model the N-terminal region 
    in Part 2 using the the broker in context of the homology model.

2.  Get template PDBs:

        ./scripts/cm_scripts/bin/get_pdbs.pl ./starting_files/query.hhr

    This script will get the PDB for each alignment from the RCSB and clean the 
    PDB files so they are compatible with Rosetta. The PDB files are named by 
    the PDB code and chain id followed by a ".pdb" file extension.

3.  Convert alignment format to 'grishin' (Rosetta uses this format):

        ./scripts/cm_scripts/bin/convert_aln.pl -format_in hhsearch -format_out grishin ./starting_files/query.hhr > query.grishin

    This script will generate the following output file:

        query.grishin

    The 'grishin' format is the standard alignment format that is used by 
    Rosetta. An example of this format is shown below:

        ## t288_ 1be9A_4
        # hhsearch_3 33
        scores_from_program: 0.000000 0.998400
        2 VPGKVTLQKDAQNLIGISIGGGAQYCPCLYIVQVFDNTPAALDGTVAAGDEITGVNGRSIKGKTKVEVAKMIQEVKGEVTIHYNKLQ
        9 EPRRIVIHRGS-TGLGFNIIGGED-GEGIFISFILAGGPADLSGELRKGDQILSVNGVDLRNASHEQAAIALKNAGQTVTIIAQYKP
        --

    The first two lines represent the identifier of the query and template 
    sequences, The string "1be9A_4" should be a unique identifier for the 
    sequence alignment where the first five characters are the PDB code and 
    chain id of the template PDB. The starting positions are given at the 
    begining of each alignment line.

4.  Fix alignments to match the PDBs:

    HHpred alignments may not match the sequence and/or amino acid positions in 
    the PDBs.  This is common and should be fixed using the following 
    application. The example command below assumes that the Rosetta release is 
    installed in the directory where this demo exists.

        ../rosetta-3.3/rosetta_source/bin/fix_alignment_to_match_pdb.macosgccrelease -in:file:alignment query.grishin -cm:aln_format grishin -out:file:alignment query_fixed.grishin -database ../rosetta-3.3/rosetta_database -in:file:template_pdb *.pdb

    This script will generate the following output file:

        query_fixed.grishin

5.  Generate fragment files:

    The next step is to generate 3mer and 9mer fragment files, which can easily 
    be done using the Robetta fragment server at 
    http://robetta.bakerlab.org/fragmentsubmit.jsp.

    Save the fragment files along with the psipred secondary structure 
    prediction file locally. We include fragment files and the psipred 
    secondary structure prediction file:

    * 3mer fragment file: `./rosetta_inputs/aat000_03_05.200_v1_3`
    * 9mer fragment file: `./rosetta_inputs/aat000_09_05.200_v1_3`
    * psipred file: `./rosetta_inputs/t000_.psipred_ss2`

6.  Run the rosetta comparative modeling application

    The Rosetta comparative modeling application 1) generates an incomplete 
    model based on the template structure by copying coordinates over the 
    aligned regions, 2) rebuilds the missing parts using the Rosetta loop 
    modeling protocol, and then 3) runs full-atom refinement of the protein 
    model using the Rosetta full-atom energy function. For more information, 
    see the documentation for the "loop modeling" and "relax" applications.

    The example command below uses the following input files:

        ./starting_files/query.fasta                 (query fasta sequence file)
        ./rosetta_inputs/aat000_03_05.200_v1_3       (3mer fragment file)
        ./rosetta_inputs/aat000_09_05.200_v1_3       (9mer fragment file)
        ./rosetta_inputs/t000_.psipred_ss2           (psipred file)
        ./query_fixed.grishin                        (the alignment file)
        *.pdb                                        (the template PDB files)

    The command-line is:

        ../rosetta-3.3/rosetta_source/bin/minirosetta.macosgccrelease @./rosetta_inputs/comparative_modeling.args -loops:frag_files ./rosetta_inputs/aat000_09_05.200_v1_3 ./rosetta_inputs/aat000_03_05.200_v1_3 none -in:file:psipred_ss2 ./rosetta_inputs/t000_.psipred_ss2 -in:file:alignment ./rosetta_inputs/query_fixed.grishin -in:file:template_pdb *.pdb -database ../rosetta-3.3/rosetta_database -out:file:silent query.out -in:file:fasta ./starting_files/query.fasta -nstruct 1

    This application will generate the following output file which contains the 
    models in a compressed "silent file" format:

        query.out

    For simplicity we set `-nstruct` to 1 in this demo, but in real world 
    applications you will want to set it higher to incease the amount of 
    conformational sampling. The standard protocol is to generate around 1,000 
    to 10,000 separate models, select the lowest 10% of models by Rosetta 
    energy, and then choose clusters using the "Cluster" application.

7.  Extract pdbs from silent file:

        ../rosetta-3.3/rosetta_source/bin/extract_pdbs.macosgccrelease -in:file:fullatom -in:file:silent query.out -database ../rosetta-3.3/rosetta_database

    This script will generate extracted PDB files with the following name 
    format where 1N9LA would be replaced with the PDB code and chain id of the 
    template:

        S_1N9LA_1_0001.pdb

    Select a homology model to use in Part 2. In real world applications you 
    would base selection on the results of clustering as briefly described 
    above.

Part 2: Broker Rigid Chunk Claimer (Bruno Correia's protocol)
-------------------------------------------------------------

Given the homology model, Part 2 describes how to model the terminal region 
using ab initio modeling in the context of the homology model generated in Part 
I.  To do this, we will use the fragment insertion sampling protocol with the 
"broker rigid chunk claimer" in Rosetta.

1.  Generate a full length model if the n-terminal region was omitted

        ../rosetta-3.3/rosetta_source/bin/full_length_model.macgccrelease -in:file:fasta ./starting_files/query.seq -loops:frag_files ./rosetta_inputs/aat000_09_05.200_v1_3 ./rosetta_inputs/aat000_03_05.200_v1_3 none -loops:frag_sizes 9 3 1 -in:file:s S_1N9LA_1_0001.pdb

2.  Make a region file that specifies the region you'd like to keep rigid (this 
    is the region that was modeling in Part I using the homology modeling 
    application). The region file included with the demo is:

        rosetta_inputs/1N9LA_143_236.region

    The format of the region file consists of a single line that contains the 
    start and end positions that span the region that will be kept rigid.  In 
    the example below positions 143-236 will be kept rigid:

        RIGID 143 236 0 0 0

3.  Create a setup file for using the RigidChunkClaimer

    Example setup file:

        ./rosetta_inputs/1N9LA.tpb

    Format of the setup file:

        CLAIMER RigidChunkClaimer
        REGION_FILE 1N9LA_143_236.region
        PDB_FILE 1N9LA.pdb
        END_CLAIMER

4.  Run the broker protocol using Rosetta:

        ../rosetta-3.3/rosetta_source/bin/minirosetta.macosgccrelease @./rosetta_inputs/broker.args -in:file:fasta ./starting_files/query.fasta -broker:setup ./rosetta_inputs/1N9LA.tpb -database ../rosetta-3.3/rosetta_database -frag3 ./rosetta_inputs/aat000_03_05.200_v1_3 -frag9 ./rosetta_inputs/aat000_09_05.200_v1_3 -nstruct 1

    For the simplicity of this protocol demo, we set -nstruct to 1. In real 
    world applications, we suggest making many more, perhaps 10000 model 
    structures, dependent on the length of the extension.

5.  Evaluate the models via clustering.  Refer to the [[clustering 
    demo|public/clustering/readme]] for more information.





Snapshots:
ftp://snapshots.rcsb.org
﻿Anchored Design
===============
 This document describes how to use the AnchoredDesign protocol, both in benchmarking and design mode.  As the protocol's components [AnchorFinder](http://www.rosettacommons.org/docs/latest/anchor-finder.html), [AnchoredPDBCreator](http://www.rosettacommons.org/docs/latest/anchored-pdb-creator.html), and [AnchoredDesign](http://www.rosettacommons.org/docs/latest/anchored-design.html) are reasonably extensively documented elsewhere, this protocol capture is meant to be used alongside that online documentation. Presented at RosettaCon2010 (in poster form) was a description of the protocol itself, plus benchmarking results, plus some early design results.  The accompanying paper ([Lewis SM, Kuhlman BA. Anchored design of protein-protein interfaces. PLoS One. 2011;6(6):e20872. Epub 2011 Jun 17.](http://www.ncbi.nlm.nih.gov/pubmed/21698112) (pubmed link)) describes only benchmarking results, but the tools to do design are described here.  A paper on design results is forthcoming.

Note that this protocol capture is somewhat focused on just the AnchorFinder portion (the least important part of the process), because the other portions are documented elsewhere but AnchorFinder largely is not.  


Contained here:
* Instructions on choosing appropriate benchmarks – AnchorFinder or otherwise
* Instructions on preparing those structures – for benchmarking or via AnchoredPDBCreator
* How to benchmark AnchoredDesign against those structures
* How to use AnchoredDesign to design interfaces
* command lines/option files, and discussion of the options

Not contained here:
* A speck of code – that lives in your Rosetta release.
* Submission scripts for running jobs.  I don't know your architecture.  It's all MPI compatible so it's not hard.

Sort-of contained here:
* The raw_documentation directory includes copies of the doxygen-style manual documentation for AnchorFinder, AnchoredDesign and AnchoredPDBCreator.  These copies are guaranteed NOT to be up to date; look in your copy of the code or the web links above instead.

Overview
--------
The purpose of the AnchoredDesign protocol is to create new protein-protein interactions using information (the anchor) borrowed from a known interaction at the same interface of one partner.  Because this protocol is intended to design protein-protein interactions, the obvious test is to see whether it can recover known structures of such interactions.  This document is massively overconcerned with benchmarking because it accompanies the paper in which AnchoredDesign is benchmarked; unless you are actually trying to replicate the benchmarking you can ignore most of those details and skip to the design tools.

The protocol modifies a loop region around an anchor in designing binders.   Selection of structures for benchmarking therefore requires interface loops with anchor regions.

The ideal anchor has several qualities:
* one or a few contiguous residues – the protocol can only have one anchor
* many interactions across an interface – this is the whole point of an anchor
* not embedded in the middle of a secondary structure element (helix/strand) – the anchor needs to be in a loop because the interface space will be sampled via loop remodeling of the anchor loop

Examples of anchors might be a phosphotyrosine inserting into an SH2 domain, a polyproline sequence binding an SH3 domain, etc.

For design, one would choose an anchor based on one's target.  For benchmarking, you are free to choose anything that has an interface loop with a good anchor.  To select benchmarking structures, I wrote the AnchorFinder protocol.  AnchorFinder has some value in highlighting which residues might make good anchors for a given target (although computational alanine scanning, not covered here, is more likely to be useful).

After finding suitable structures with the help of AnchorFinder, the next step is to pick anchors and loops out of those structures, and in general prepare them for Rosetta's use (removing solvent atoms, etc).  At this point you're ready to run AnchoredDesign.


Compiling AnchorFinder
----------------------
(This section applies to benchmarking only)

Your Rosetta code distribution should include an application called AnchorFinder. If you wish to search large numbers of PDBs for potential anchors (I searched a local copy of the entire PDB structure set), then you will wish to modify the code slightly before running it.  Running any part of Rosetta against huge numbers of unprepared, straight-from-the-PDB structures is challenging because the PDB reader in Rosetta is not robust against nonstandard file formats, etc.

To compile AnchorFinder such that it will be robust, examine the manual documentation on RobustRosetta (also included in this protocol capture).  Briefly, this documentation describes changes that A) make Rosetta slower (thus they aren't on by default) and B) cause it to throw C++ exceptions when it hits errors instead of crashing.  The job distributor catches the errors, skips the bad structures, and continues.  You must recompile after making these changes.

You do not want to use compiled executeables OTHER than AnchorFinder with these changes made – they will significantly slow the code down.  AnchorFinder is quite fast so it's not a problem.

When running AnchorFinder, watch your memory usage.  When I used it, there was a patch in the JobDistributor which deleted starting poses for PDBs that had already been processed.  This patch was rejected by the community and since been replaced by a different patch to do the same thing; AnchorFinder is a run-once sort of thing so it has not been tested against the new method.

Using AnchorFinder
------------------
(This section is minimally relevant if not benchmarking)

At this point you should have a compiled copy of AnchorFinder with the necessary changes to the code.  You can then list your PDBs in one or many -l files (or -s) for use in Rosetta.  The format for Rosetta's -l flag is one path per line:

    A.pdb
    B.pdb
    C.pdb
    ...

Depending on your available architecture, it may be better to split the run up into many -l on separate processors.  I don't know what's best for you.

If you want to do the whole PDB – it's a good idea to skip the largest PDB files ahead of time, particularly ribosome structures.  These take a very long time to process through the PDB reader, and due to heavy nucleic acid content are skipped anyway.  You can also either toss the NMR structures ahead of time or use the -obey_ENDMDL flag to only read the first model.

AnchorFinder will automatically remove nonprotein atoms from the Poses before examination.  It also skips anything that is monomeric, has no protein residues, or smaller than 20 residues after processing.

It will then look through the structures searching for regions with certain command-line-defined characteristics.  These characters are:
* length of windows for consideration - 4 or 5 or 6, etc, contiguous residues at a time?  This flag is -window_size.  I suggest 5 residue windows.
* What fraction of this window should have loop secondary structure as assigned by DSSP?  -loopness controls this.  It takes it as a decimal between 0-1, I suggest 0.6 (which translates to 3/5 residues for a 5 residue window) to 1 (all residues loop).
* How many cross-interface interactions per residue should the window have?  I suggest a minimum of 4.  This translates to 20 (redundancy included) cross-interface interactions for a 5 residue window.  By redundancy, I mean residues 43 and 44 on chain A can both interact with residue 234 on chain B and it will count as two interactions.  Specify this with -nbrs_per_residue.
* What file name should the good interactions be printed to?  I leave it as an exercise to the reader to pick their own file name.  Specified with -bestoutfile; defaults to goodfile.out.

Running AnchorFinder, while not particularly slow, is still something you only want to do once.  The defaults suggested above produce lots of output, which can then be further processed quickly without reloading PDBs.  To expedite this, AnchorFinder produces two levels of output.  All residues have their data printed to a file named (pdbname).data – you can reprocess this to get data for differing window lengths, loopnesses, etc.  Windows passing the loopness and interactions filters are printed to the specified output file.

A suggested options file for AnchorFinder is available with this document.

Interpreting AnchorFinder Results
---------------------------------
(This section is minimally relevant if not benchmarking)

After you've run AnchorFinder, you'll have a fairly large pile of output: pdbname.data for all pdbs, plus goodfile.out for the better windows.

`pdbname.data` looks like this:

Rows are residues, columns are chains, data are neighbors in that chain for each residue

    residue chain   PDBdata DSSP    1       2
    1       1       2 D     L       7       0
    2       1       3 D     L       10      0
    3       1       4 D     L       14      0
    ...

The columns are residue and chain in Rosetta numbering, residue/chain in PDB numbering, DSSP value, and then N columns for the N chains in the protein.  The number in those columns is the number of cross-interface neighbors on that chain for that position.

`goodfile.out` looks like this:

    PDB pdb2vk1 window 45 loopness 5 nbrs 0 28 0 0 start 46 A pymol select pdb2vk1 and chain A and resi 46-50
    PDB pdb2vk1 window 108 loopness 5 nbrs 0 25 0 0 start 109 A pymol select pdb2vk1 and chain A and resi 109-113
    PDB pdb2vk1 window 109 loopness 5 nbrs 0 36 0 0 start 110 A pymol select pdb2vk1 and chain A and resi 110-114
    PDB pdb2vk1 window 110 loopness 5 nbrs 0 46 0 0 start 111 A pymol select pdb2vk1 and chain A and resi 111-115
    PDB pdb2vk1 window 111 loopness 5 nbrs 0 46 0 0 start 112 A pymol select pdb2vk1 and chain A and resi 112-116
    PDB pdb2vk1 window 112 loopness 5 nbrs 0 47 0 0 start 113 A pymol select pdb2vk1 and chain A and resi 113-117

Each line identifies the PDB, the window number, its loopness, its number of neighbors on each chain in the PDB (variable # of columns), the starting residue PDB numbering for the window, and a Pymol selection for the window.

Inputs and outputs for this stage from a convenience sample (PDBs 3cy?) are included with this protocol capture.

At this point, the data is yours to play with.  I searched for windows with large numbers of neighbors on only one chain using sifter.py (included), then sorted for those with the largest number of neighbors (sort -n -k1 `input`).  After that it was all manual filtering to choose structures for the benchmarks.

Choosing Loop and Anchor — Benchmarking
---------------------------------------
(This section applies to benchmarking only)

OK, so you ran AnchorFinder, looked at the results, and/or picked what protein you want to run through AnchoredDesign.  How do you choose a loop/anchor?

If you ran AnchorFinder, look at the AnchorFinder result lines that came up as good:

    92 PDB pdb1zr0 window 526 loopness 5 nbrs 0 0 92 0 start 13 D pymol select pdb1zr0 and chain D and resi 13-17

Load this PDB into pymol (1zr0.pdb) and activate the suggested selection.  You'll see that it is in a surface loop of one partner which sticks an arginine straight into its binding partner – a perfect anchor.  (This is a chosen example; not all AnchorFinder hits are this nice.)

Choosing the anchor is entirely up to human effort; here the arginine 15 is an obvious choice.

For choosing loops, I just traveled up and down the chain in both directions until I hit secondary structure, significant backbone-backbone hbonding, or the protein core.  Here I'd choose a loop of D10 to L17 – more N-terminal than that affects the core, and more C-terminal affects a sheet.

Anchor and loop file specifications are included in the release documentation and the examples here.

Note that for the included example, the PDB has been renumbered from 1.  Scripts to do this are occasionally included with Rosetta distributions and not included here.  It will be convenient to also remove waters, ligands, etc.

If you are doing benchmarking, skip to the [Running AnchorDesign](#Running 
AnchorDesign) section.

Choosing a System and Anchor — Design
-------------------------------------
(This section applies to the design case only)

In the design case, you will be choosing your proteins based on what you want designed.  Your target is forced by what targets:
* have crystal structures, and
* are related to your biological problem.

Choosing an anchor then requires:
* a cocrystal of your target with some partner from which to source the anchor.

You can run this cocrystal through AnchorFinder and let it suggest anchors to you, but for one structure you can just look at it yourself.  Look for loops on the partner that insert into the target, or do computational alanine scanning, or examine the literature for mutations that disrupt the interface.

Choosing a Scaffold, Design Positions, and Loops — Design
---------------------------------------------------------
(This section applies to the design case only)

In the design case, you will be replacing your target's partner with some new scaffold to form a mostly de novo interface.  Your scaffold must meet a few requirements:
* flexible, mutateable surface loops (for AnchoredDesign to modify)
* experimentally tractable (hey, your funeral if it's not)
* whatever other functionality you need for your desired design

The protocol was written with the fibronectin monobody scaffold in mind.

Choosing which loops are flexible is dependent on biological knowledge of the scaffold.  In fibronectin's case, many papers have been published establishing the mutability of the BC and FG loops.

Choosing which positions are designable is similarly dependent on your scaffold.  AnchoredDesign carries the assumption that the non-anchor loop positions are designable, and non-loop positions are not, but nothing in the code enforces that.  Use a resfile (documented with the manual) to specify which positions are designable.  The code will automatically prevent design of the anchor (you can turn that off).  The code will automatically prevent design of positions that are not close to either the interface or a flexible loop (you cannot turn that off), so take care in specifying designable positions on opposite faces of your protein.  Proximity is redetermined at each design opportunity so positions peripheral to the interface may not be designed regularly.

Choosing Loop Lengths and Anchor Position
-----------------------------------------
(This section applies to the design case only)

OK, so you know which scaffold to use, and which anchor, and which target.  You are ready to create your starting structure for AnchoredDesign, in which the anchor will be inserted into the scaffold, and the anchor will be aligned properly to the target, dragging the scaffold and target together.  The protocol used for this is called AnchoredPDBCreator; further details are below.

One important part of conformational space that AnchoredDesign cannot search is the space of loop lengths and anchor positions.  You may want to try, for a loop of length N, all combinations of loops of length N-3 to N+3, or even more for long loops.  As you are designing the loop to form an interface, there is no reason to believe its native length is particularly relevant.  You will have to do this searching at this stage: create starting structures for all loop lengths, run them all through AnchoredDesign, and pick off the best ones later.

Loops can be shortened directly by just deleting residues mid-loop before handing the scaffold to AnchoredPDBCreator – it can insert a 3 residue anchor into a 6 residue window, and close the gap.  Loop lengthening must be done externally.  One way to lengthen loops is to manually modify a PDB to contain enough residues in the loop (copy-and-paste a residue, renumber as necessary), then use the loop_modeling executeable's build_initial mode to close the loop.  Further instructions are included in their own folder in this packet.

A paired space is anchor placement space.  Besides choosing which anchor to use (try several), exactly where it is placed within a loop can vary.  For a loop of length 7, and an anchor of length 2, (assuming a flexible residue on each side), you have the following 4 choices:

    X = scaffold
    - = loop
    A = anchor
    X1234567X
    X-AA----X
    X--AA---X
    X---AA--X
    X----AA-X

Again, this space is not searched by AnchoredDesign and must be searched by trying all the inputs.

Using AnchorPdbCreator
----------------------
(This section applies to the design case only)

AnchoredPDBCreator is the protocol which assembles an anchor, scaffold, and target into a starting structure for AnchoredDesign.  Its code documentation is included in this packet.

Briefly, AnchoredPDBCreator takes as input 4 files:
* The target structure, as a PDB, with the partner removed
* The anchor structure, drawn from the cocrystal with the target, containing ONLY the residues being used as an anchor, as a PDB
* The scaffold structure, as a PDB, with loop residues added/deleted as desired
* A scaffold_loop specification, which declares which residues in the scaffold are flexible and where the anchor insertion should occur.

It is ABSOLUTELY VITAL to recognize that AnchoredPDBCreator does NOT produce interfaces, it only produces starting structures for AnchoredDesign.  It is entirely plausible that its structures will have the target and scaffold totally eclipsed.  This is fine, AnchoredDesign will fix it.

AnchoredPDBCreator's results should be interpreted by analyzing ONLY the closure of the anchored loop.  Use the result with the best loop geometry.  Loop geometry can be measured by examining the LoopAnalyzerMover output tagged to the end of result PDBs:

LoopAnalyzerMover: unweighted bonded terms and angles (in degrees)

    position phi_angle psi_angle omega_angle peptide_bond_C-N_distance rama_score omega_score dunbrack_score peptide_bond_score chainbreak_score
     pos phi_ang psi_ang omega_ang pbnd_dst    rama  omega_sc dbrack pbnd_sc   cbreak
      17  -106.8   175.8     178.2    1.322   0.998    0.0342   7.01   -2.68   0.0182
      18  -82.33   64.67    -178.5    1.329   0.211    0.0217   3.11   -3.42   0.0203
      19  -83.63   149.4     177.2    1.329   -1.07    0.0795      0   -3.43    0.584
      20  -75.25   171.1    -178.7    1.329  -0.264    0.0161  0.348   -3.43   0.0151
      21  -58.53  -42.95     174.6    1.329   -0.58     0.294      0   -3.43      2.7
      22  -76.02   159.9    -179.8    1.326  -0.811  0.000404   0.97   -3.45   0.0424
      23  -72.63   130.1     179.4    1.325   -1.29   0.00372   0.24   -3.46   0.0281
      24  -94.91   116.5     179.8    1.323   -1.21   0.00028  0.721   -3.45   0.0694
      25  -65.42   150.7     179.4    1.335   -1.58     0.004      0   -3.32     1.38
      26  -64.68   147.9     179.1    1.323   -1.45    0.0079   1.61   -3.32    0.211
      27  -56.44  -66.68      -180    1.329    1.34  8.08e-30   7.87   -3.43 2.37e-05
      28  -124.4  -56.48     177.6    1.329    2.08    0.0568  0.608   -3.43   0.0533
      29  -124.1   28.78    -177.7    1.264   0.341    0.0542   2.39    2.65     2.07
      30   81.57  -134.3    -176.4    1.329      20     0.126   5.06    2.65    0.128
      31  -112.9   147.2     172.7    1.318  -0.744     0.538  0.534   -3.35     1.38
    total_rama 15.9674
    total_omega 1.23676
    total_peptide_bond -38.3223
    total_chainbreak 8.70689
    total rama+omega+peptide bond+chainbreak -12.4113

    LAM_total -12.4113

In this particular example, position 29 is clearly problematic: the peptide bond distance is too short, as reported by the pbnd_dst, pbnd_sc, and cbreak columns.

You should be running AnchoredPDBCreator for at least 100 trajectories before choosing a starting structure.

Running AnchorDesign
--------------------
If you are benchmarking, the crystal structure of the complex is the appropriate input for AnchoredDesign.  If you are designing, the best result from AnchoredPDBCreator is your starting structure.

The input files for AnchoredDesign provide an example with 1zr0 for running AnchoredDesign.  It is a heterodimer so you can pretend it was AnchoredPDBCreator sourced if you want.  (You can also look in the AnchoredDesign integration test at test/integration/tests/AnchoredDesign for such an input).

* The anchor file specifies which residues form anchor.
* The PDB file is the pdb.
* loopsfile_extended is for benchmarking – the extension column is true, which tells AnchoredDesign to forget the starting loop conformation before sampling
* loopsfile_native tells AnchoredDesign NOT to forget the starting loop conformation – this is probably what you would use for design, although there is no reason you can't reject the starting loop conformation for design.  (This is cheating for benchmarking, but useful for making relaxed natives)
* options is the command line options file.  The active options are for a simple “does it run” test; parameters for a longer running “real” test are included.
* A resfile is necessary if you wish to design, and is only active in that option set.  

Interpreting AnchorDesign — Benchmarking
----------------------------------------
(This section applies to the benchmarking case only)

If you are duplicating the benchmarking results, you passed the rmsd flag. AnchoredDesign will have output a lot of RMSD values allowing you to determine the performance of the protocol against the structures you chose to benchmark.  The paper describes the score versus RMSD metrics used to determine quality (including the I_sup_bb_RMSD, ch2_CA_RMSD, and loop_CA_sup_RMSD.  The structures themselves don't really matter; you are ensuring that the low-scoring structures have low RMSD.

Interpreting AnchorDesign — Design
----------------------------------
(This section applies to the design case only)

In the design case, the other fields of the AnchoredDesign output come in to play. There are three classes of output:
* scorefunction terms
* LoopAnalyzerMover output,
* InterfaceAnaylzerMover output.

Generally, you should rank your structures according to total_score (the Rosetta scorefunction).  This tells you what Rosetta thinks is best.

Next, you use the LoopAnalyzerMover output (described above) and InterfaceAnalyzerMover output to determine which structures have flaws not caught by total_score.  Toss structures that those filters think have problems.  Pick the ones you think are best, order the DNA, and pray.  When it works great, feel free to send me kudos, citations, or money!

InterfaceAnalyzerMover
----------------------
InterfaceAnalyzerMover output looks like this:

    Residues missing H-bonds:
    Residue 	 Chain 	 Atom 
    38 	 A 	 NE2
    101 	 A 	 OE1
    248 	 A 	 O
    250 	 A 	 O
    344 	 B 	 N
    384 	 B 	 O
    477 	 B 	 O

    pymol-style selection for unstat hbond res 
    select start_5411_unsat, /start_5411//A/38+101+248+250+ + /start_5411//B/344+384+477+

    pymol-style selection for interface res 
    select start_5411_interface, /start_5411//A/31+32+33+34+35+36+37+38+39+40+41+54+56+57+59+60+61+62+64+65+66+92+95+98+99+100+101+102+103+106+194+195+224+225+226+227+228+229+230+247+248+249+250+251+252+253+265+ + /start_5411//B/314+315+316+317+318+319+320+321+322+323+324+337+339+340+342+343+344+345+347+348+349+350+375+378+379+381+382+383+384+385+386+389+473+476+477+478+479+480+481+482+488+506+507+508+509+510+511+512+513+

The first section documents where Rosetta thinks there are unsatisfied hydrogen bonds at the interface.  This code is known to be oversensitive to missing bonds, but it's better than nothing.

The next sections print PyMOL selections for interface residues for easier visualization.

InterfaceAnalyzerMover also includes columns into the scorefile:

    dSASA_int 2396.33
    dG_separated -35.3379
    dG_separated/dSASAx100 -1.47467
    delta_unsatHbonds 7
    packstat 0
    dG_cross -27.6963
    dG_cross/dSASAx100 -1.15578
    AllGly_dG -2.83564
    cen_dG -10.3844
    nres_int 96
    per_residue_energy_int -1.15006
    side1_score -361.478
    side2_score -267.353
    nres_all 520
    side1_normalized -1.25513
    side2_normalized -1.15238
    complex_normalized -1.74616
    hbond_E_fraction 0.368537

Most of these are experimental and not useful (and not part of AnchoredDesign; InterfaceAnalyzerMover has other clients).  The useful ones are dG_separated/dSASAx100, which measures the Rosetta energy of binding per unit area of SASA (scaled by a factor of 100).  This ensures you pick an interface that is energetic for its size, not large but sloppy.
Beta-3-peptide design
=====================

This demo illustrates a protocol to redesign and model beta-3-peptides in Rosetta.
It was written in July 2012 by Fang-Chieh Chou (fcchou at stanford dot edu).
Detailed application of the method is described in the following paper:

* Molski, M.A., Goodman, J.L., Chou, F.C., Baker, D., Das, R., and Schepartz, A. (2012). Remodeling a beta-peptide bundle. Chemical Science.

Running the Demo
----------------

The example input files are in rosetta_input.
The beta-3-peptide residues is less commonly used so is not turned on by default. Therefore the following setup steps are required prior to running this demo.
Before running the command lines, edit the file `rosetta_database/chemical/residue_type_sets/fa_standard/residue_types.txt`, find the section "Beta-peptide Types" and uncomment all the items listed below:

* residue_types/beta-peptide/B3A.params
* etc.

The example command lines are given in rosetta_input/cmdline.
Cd into rosetta_input/ and type `sh cmdline` to run the commands.
After the job finished successfully, you should see the output pdb files in the corresponding folders and log files in the current folder.
Sample outputs are given in the example_output folder for comparison.

Command-lines
-------------

1. Redesigning
   ```
   beta_peptide_modeling.linuxgccrelease -database  ../../../../rosetta_database -force_field beta_peptide -native redesign/acdy_LLLL_LLLL.pdb -algorithm redesign -ex1 -ex2 -packing::pack_missing_sidechains false -packing::extrachi_cutoff 0 -repack_res 2 5 8 11 14 17 20 23 -n_repeat 4 -repeat_size 24
   ```

   `-native` gives the starting model.

   `-repack_res` specify the residue being redesigned/repacked.

   `-n_repeat` and `-repeat_size` set up the symmetry based repacking. In this example, each residue n listed in `repack_res` and n+24, n+48, n+72 (4 copies total) are repacked/redesigned. Also each set of 4 residues are enforced to be the same residue type and have the same rotameric conformation.

2. Repacking
   ```
   beta_peptide_modeling.linuxgccrelease -database ../../../../rosetta_database -force_field beta_peptide_soft_rep_design -native repack_and_minimize/acdy_LFFL_LFFL.pdb -algorithm repack -ex1 -ex2 -packing::pack_missing_sidechains false -packing::extrachi_cutoff 0 -repack_res 2 5 8 11 14 17 20 23 -n_repeat 4 -repeat_size 24
   ```
   Similar to above, but only repack the original residues.

3. Minimization
   ```
   beta_peptide_modeling.linuxgccrelease -database  ../../../../rosetta_database -force_field beta_peptide -native repack_and_minimize/acdy_LFFL_LFFL_repack.pdb -algorithm minimize -ex1 -ex2 -packing::pack_missing_sidechains false -packing::extrachi_cutoff 0 
   ```
   Minimize a given pdb, return the final model and Rosetta energies.

Extra flags and score files
---------------------------

`-score::no_smooth_etables` true Use the old Rosetta attractive/repulsive score without smoothing.

`-no_symmetry` true Disable the symmetry-based redesigning.

`-force_field beta_peptide_mm` / `-force_field beta_peptide_soft_rep_mm` Use the molecular mechanics-based Rosetta scoring files. Replace the corresponding `-force_field` flags with the MM ones in the above command lines. Also, extra flags `-apply_dihedral_cst false` should be used to turn off dihedral constraint when MM potential is used.

# Multistate Design of Antibodies
This demo includes input files for running multistate design using 10 states.
It can either be run on a single processor using the
```
   mpi_msd.default.linuxgccrelease (or, more generally the mpi_msd.default.{os}{compiler}{release/debug} executable)
```
or by distributing those states over multiple processors using the
```
   mpi_msd.mpi.linuxgccrelease (or, more generally the mpi_msd.mpi.{os}{compiler}{release/debug} executable)
```
and by launching the job with between 1 and 10 CPUs.  Using 5 CPUs would place two states
on each processor.

Files to look at first:
```
command :         contains an example command line which you will have to modify for your system
command_mpi:      contains an example command line for running the mpi executable which you will have to modify for your system
1USM_het.flags:   contains the set of flags read by the command line in the command/command_mpi files.
fitness.daf:      the file which declares the set of states included in this design task, and the fitness function itself
entity.resfile:   the resfile which declares the accessible regions of sequence space
```
# Peptide Specificity

## What is this?
This is the demo for the pepspec and pepspec_anchor_dock applications.

## What does it do?
First, a single peptide "anchor" residue will be docked onto the surface of a peptide binding protein (c-CRK SH3 domain) in which the original peptide has been deleted. This docking will be repeated to generate an ensemble of peptide anchor-docked structures. To do this, the pepspec_anchor_dock application uses the relative orientation of peptide anchor residues in three homologous peptide complex structures (1CKB, 1N5Z, 1OEB).
Next, these anchor residue docked structures are used by the pepspec application to design putative binding peptides on the surface of the SH3 domain. Peptide residues are added to the anchor residue while peptide sequences and structures are simultaneously explored. Structures and sequences are saved for later processing.

# Demo Summary 
1. generating an ensemble of anchor prolines docked to an SH3 scaffold, then

2. exploring the sequence specificity of peptides designed from those docked proline anchors

3. generate a specificity PWM from the designed peptides

# Demo Detail
1. run the anchor docking protocol to generate 10 structures of proline docked to 1CKA.align.nopep.pdb
    ```
    ~/mini/bin/pepspec_anchor_dock.linuxgccrelease /path/to/minirosetta_database @dock.args
    ```

    This will generate 10 pdbs, a pdblist file, and a cst (contraint) file. These files are used in the next stage.

    - 1cka.docked_[1-10].pdb are the anchor-docked structures of the input structure 1CKA.align.nopep.pdb
    - 1cka.docked.pdblist is simply a list of the pdb files generated
    - 1cka.docked.cst is a peptide constraint file

2. run the peptide design protocol:
    ```
    ~/mini/bin/pepspec.linuxgccrelease -database /path/to/minirosetta_database @spec.args
    ```

    This will generate a folder "1cka_spec.pdbs" with some designed peptide - protein complexes and
a 1cka.spec file with all the peptide sequences, energies, and some other values.

    This *.spec file can be used for deriving a PWM using the script

3. Genrate a normalized PWM (position-weight-matrix)
    ```
    ~/mini/analysis/apps/gen_pepspec_pwm.py 1cka_spec.spec 3 0.1 binding-prot_score /path/to/minirosetta_database/pepspec_background.binding-prot-0.1.pwm
    ```

    This will sort the peptide sequences in 1cka_spec.spec, in which there are 3 residues n-term to the anchor residue, by the "binding-prot_score" score term (one of many scores calculated in 1cka_spec.spec), then filter out the lowest-scoring 10% peptide sequences, and use these sequences to construct a peptide PWM based on the position-specific amino acid frequencies. This resultant PWM is then normalized by a background PWM (pepspec_background.binding-prot-0.1.pwm) in order to eliminate some residue frequency biases caused by artifacts of the rosetta score function.

    The script creates a PWM file (1cka_spec.binding-prot_score.0.1.pwm), a list of the sequences from which that PWM was constructed (1cka_spec.binding-prot_score.0.1.seq), and a normalized PWM file (1cka_spec.binding-prot_score.0.1.norm.pwm).

    For a real, production-level run, I would reset the values of these flags:
    ```
    <dock.args>
    -pepspec::n_peptides 100
    
    <spec.args>
    -pepspec::n_peptides 1000
    -pepspec::diversify_lvl 5
    ```

**Note:**
	In a real production-level run, the peptide backbone coordinate constraints (e.g. 1cka.docked.cst) may be generated from only one or two homologue complexes. If this is the case, or you have reason to believe that peptides may bind with backbone conformations not represented in the homologue set, then use of the constraints may prevent your simulation from ever generating the correct peptide backbone structures.
	To avoid use of constraints, simply do not include a reference to the constraint file in you command line arguments. (e.g. delete the line "-pepspec::homol_csts 1cka.docked.cst" from spec.args).
	If you're not using constraints, you need to do much, MUCH more sampling (i.e. an order of magnitude more).
The unconstrained sequence+structure space gets big fast.

HBS Design Demo
===============

Written by Kevin Drew (Bonneau Lab), kdrew at nyu dot edu

This demo shows how to run the hbs_design application.  A hydrogen bond 
surrogate (HBS) is a helical mimetic scaffold used for inhibiting protein 
interactions. The demo shows the design of an hbs inhibitor for the MDM2-P53 
protein interaction.

Algorithm
---------

1. Pertubation phase: rigid body movement of hbs wrt target
2. Design phase: design user specified residues on hbs scaffold and minimize
3. Repeat 10x

Command
-------

    hbs_design<.exe> -database <path to your database> @input/flags

Input Files
-----------

* `./input/flags`: User specified options.
* `./input/mdm2_hbs.pdb`: Input structure where target is chain 1 and hbs is 
  chain 2

Options
-------

* `-hbs_design_positions`: Residues on hbs to design (numbering is relative to 
  hbs, for example 3 is the third residue on hbs), default repacks with no 
  design.
* `-pert_num`: Number of pertubations during pertubation phase, default 10, 
  production 100.
* `-design_loop_num`: Number of pertubation + design cycles, default 10, 
  production 10.
* `-nstruct`: For production runs, use 1000.

Pre-processing
--------------

...

Post-processing
---------------

Similar to other multi chain design protocols, the ddG is computed and is a 
good indicator of a good design.  First sort by total score, take top 5 percent 
and then sort by REPACK_ENERGY_DIFF (ddG).

Limitations
-----------

This app is inflexible to adjusting Monte Carlo temperatures, score functions, 
degree of rigid body pertubations, designing noncanonical amino acids, etc. The 
app also requires the hbs is close to a plausible binding mode with respect to 
the target.

AbInitio fold-and-dock of peptides using FlexPepDock
----------------------------------------------------
This demo illustrates how to run FlexPepDock ab-initio folding and docking of a peptide onto its receptor. The FlexPepDock ab-initio protocol is designed to generate high-resolution models of complexes between flexible peptides and globular proteins, given the approximate location of the peptide binding site. The ab-initio procol samples both rigid-body orientation and torsional space of the peptide extensively. No prior knowledge about the peptide backbone is necessary, as the protocol uses fragments to sample peptide backbone conformational space rigorously.

Protocol overview
-----------------
The input to the ab-initio protocol is a model of the peptide-protein complex in PDB format,starting from arbitrary (e.g., extended) peptide backbone conformation. It is required that the peptide is initially positioned in some proximity to the true binding pocket, but the exact starting orientation may vary.
Preliminary steps: (1) Generation of fragment libraries for the peptide sequence, including 3-mer, 5-mer and 9-mer fragments. (2) Pre-packing of the receptor and peptide to remove internal clashes that might confuse ranking.
Main protocol: Step 1: Monte-Carlo simulation for de-novo folding and docking of the peptide over the protein surface in low-resolution (centroid) mode, using a combination of fragment insertions, random backbone perturbations and rigid-body transformation moves. Step 2: The resulting low-resolution model is refined with FlexPepDock Refinement. As in the independent refinement protocol, the output models are then ranked based on their energy score, after their clustering for improved coverage of distinct conformations.

Refinement vs. ab-initio protocol
---------------------------------
The Refinement protocol is intended for cases where an approximate, coarse-grain model of the interaction is available that is close to the correct solution both in Cartesian and dihedral (phi, psi) space. The protocol iteratively optimizes the peptide backbone and its rigid-body orientation relative to the receptor protein including on-the-fly side-chain optimization (look at /demos/refinement_of_protein_peptide_complex_using_FlexPepDock/ to learn how to run refinement of protein-peptide complexes).
The ab-initio protocol extends the refinement protocol considerably, and is intended for cases where no information is available about the peptide backbone conformation. It simultaneously folds and docks the peptide over the receptor surface, starting from any arbitrary (e.g., extended) backbone conformation. It is assumed that the peptide is initially positioned close to the correct binding site, but the protocol is robust to the exact starting orientation. The resulting low-resolution models are refined using the FlexPepDock Refinement protocol.

Running the FlexPepDock ab-initio protocol
------------------------------------------
1. Generate an initial complex structure: An initial model can be built by placing the peptide in close promity to the binding site in an arbitary conformation. In this demo we have provided a starting structure with a peptide in extended conformation (2A3I.ex.pdb). Our goal is to optimize this structure using ab-initio FlexPepDock, towards a near-native model with a helical peptide conformation. Both the native structure (2A3I.pdb), as well as the starting structure (2A3I.ex.pdb) are provided in the input directory.

2. Prepack the input model: This step involves the packing of the side-chains in each monomer to remove internal clashes that are not related to inter-molecular interactions. The prepacking guarantees a uniform conformational background in non-interface regions prior to refinement. The prepack_flags file contains the flags for running the prepacking job. The run_prepack script will run prepacking of the input structure 2A3I.ex.pdb located in the input directory.

You need to change the paths of the Rosetta executables and database directories in the run_prepack script (also for run_refine; see below).

  ROSETTA_BIN="rosetta/main/source/bin"
  ROSETTA_DB="rosetta/main/database/"

After changing the paths, run the run_prepack script as:
   $./run_prepack

The output will be a prepacked structure, 2A3I.ex.ppk.pdb, located in the input directory; a scorefile named ppk.score.sc and a log file named prepack.log file located in the output directory. This prepacked structure will be used as the input for the ab-initio modeling step.

3. Create 3mer, 5mer & 9mer (peptide lingth >=9) fragment libraries: The scripts necessary for creating fragments are provided in the fragment_picking directory.
    a. Go to the fragment_picking directory.
    b. Save the peptide sequence in the xxxxx.fasta file.
    c. Run the make_fragments.pl script to generate the PSIPred secondary structure and PSI-Blast sequence profiles. You need to chnage the paths in the upper section of the make_fragments.pl file.
       Run as $perl make_fragments.pl -verbose -id xxxxx xxxxx.fasta
       This will create xxxxx.psipred_ss2, xxxxx.checkpoint along with other files.
    d. Run the executable fragment_picker.linuxgccrelease to create the frags. The flags are provided in the flags file and fragment scoring weights are provided in the psi_L1.cfg file.
    Run as $ROSETTA_BIN/fragment_picker.linuxgccrelease -database $ROSETTA_DB @flags >log
    e. Change the fragment numbering using shift.sh script.
    Run as $bash shift.sh frags.500.3mer X >frags.3mers.offset ; where X is the number of residues in the receptor. Do the same for 5mer and 9mer frags
    The offset fragment files will be used as input to the FlexPepDock ab-initio protocol. Put them in the input/frags directory.

4. Ab-initio folding and docking of the prepacked model: This is the main part of the protocol. In this step, the peptide backbone and its rigid-body orientation are optimized relative to the receptor protein using the Monte-Carlo with Minimization approach, including periodic on-the-fly side-chain optimization. The peptide backbone conformational space is extensively sampled using fragments derived from solved structures. The file abinitio_flags contains flags for running the ab-initio job. The run_abinitio script will run ab-initio modeling of the prepacked structure generated in the prepacking step located in the input directory.

After changing the Rosetta related paths run the run_abinitio script as:
    $./run_abinitio

The output will be an optimized structure (2A3I.ex.ppk_0001.pdb) located in the output directory; a scorefile named abintio.score.sc and a log file named abinitio.log, located in the output directory. This script has to be modified to run on a cluster during a production run (see below).


Specific changes needed for a production run
--------------------------------------------
For a production run it is recommended to generate large number of decoys (~10,000 to 50,000). In such a case you can run the job on a cluster. It is advided to use silent output format in such scenario to save space (See https://www.rosettacommons.org/manuals/rosetta3.1_user_guide/app_silentfile.html for details). Include the following lines to the abinitio_flags file:

-out:file:silent_struct_type binary
-out:file:silent decoys.silent

This will create the decoys.silent file containing data related to all the decoys in a compressed format. You can extract speicific decoy using the extract_pdbs.linuxgccrelease executable.
For example:
  $ROSETTA_BIN/extract_pdbs.linuxgccrelease -database $ROSETTA_DB -in:file:silent decoys.silent -in:file:tags 2A3I.ex.ppk_1234.pdb


Along with changes in the flags file you need to modify the run_abinitio file to run on a cluster. The file run_abinitio_slurm is the modified version of run_abinitio adapted to run on a slurm cluster. You should ask your cluster manager for relevent changes required.


Post Processing after a production run
--------------------------------------
In order to diversify our prediction, we cluster the results and select representative models. A clustering scripts is provided in the clustering directory. It will cluster the top 500 decoys based on a cutoff radius of 2.0 Angstrom, and select for each the top-scoring member (according to reweighted score,  reweighted_sc). A top scoring member from each cluster is reported in the file  cluster_list_reweighted_sc_sorted.

Runs as
$bash cluster.sh 2.0 ../input/2A3I.ex.pdb ../output/decoys.silent

Further information
-------------------
Detailed documentation on ab initio FlexPepDock is available under: https://www.rosettacommons.org/docs/latest/application_documentation/docking/flex-pep-dock.
Please cite: Raveh B, London N, Zimmerman L, Schueler-Furman O (2011) Rosetta FlexPepDock ab-initio: Simultaneous Folding, Docking and Refinement of Peptides onto Their Receptors. PLoS ONE 6(4): e18934. doi: 10.1371/journal.pone.0018934

CS-Rosetta RNA Demo
===================

The cs_rosetta_rna application refines (minimizes) and scores a RNA structure under the hybrid CS-ROSETTA-RNA all-atom energy function:

    E(hybrid) = E(Rosetta) + E(shift)

where E(Rosetta) is the standard Rosetta all-energy function for RNA [1], and E(Shift) is the chemical shift pseudo-energy term [2].
Input RNA PDB structures can be generated by FARFAR [1] and/or Stepwise Assembly [3] structure modeling methods, or can be an experimental NMR or crystallographic structure.

Example inputs
--------------

1.  `rosetta_inputs/GA-AG_mismatch/1MIS_NMR.pdb`  
    NMR PDB structure of a tandem GA:AG mismatch internal loop.

2.  `rosetta_inputs/GA-AG_mismatch/1MIS_exp_1H_chem_shifts.str`  
    Experimental non-exchangeable 1H chemical shift data for the
    tandem GA:AG mismatch interna loop.

    Each line in this file represent one chemical shift data point and 
    contains the following nine space-delimited columns (based on the STAR
    v2.1 format):

        col 1: Atom_shift_assign_ID (INT)
        col 2: Residue_author_seq_code (INT)
        col 3: Residue_seq_code (INT)
        col 4: Residue_label (STRING)
        col 5: Atom_name (STRING)
        col 6: Atom_type (STRING)
        col 7: Chem_shift_value (FLOAT)
        col 8: Chem_shift_value_error (STRING)
        col 9: Chem_shift_ambiguity_code (STRING)

    Note that residue_seq_code (col 3), residue_label (col 4), and 
    atom_name (col 5) should be consistent with the data in the PDB file. 
    Also, col 8 and col 9 are currently not used internally by
    the cs_rosetta_rna application.

3. `rosetta_inputs/GA-AG_mismatch/1MIS_params`  
    Parameters file (in FARNA format) for the tandem GA:AG mismatch.

    Also, the input data files for all 23 RNA motifs benchmarked in ref. [2]
    are provided in the Supplemental Data Zip file, available at the following
    URI: http://dx.doi.org/10.1038/nmeth.2876

Example command-lines
---------------------

1.  Refine (minimize) and score a PDB structure under the hybrid CS-ROSETTA-RNA all-atom energy function:

        <path-to-rosetta-bin>/cs_rosetta_rna.<release> \
            -mode minimize_pdb \
            -pdb <input_pdb> \
            -score:rna_chemical_shift_exp_data <exp_cs_data_file> \
            -params_file <input_param_file> \
            -analytic_etable_evaluation false

2.  Score (but not refine) a PDB structure under the hybrid CS-ROSETTA-RNA all-atom energy function:

        <path-to-rosetta-bin>/cs_rosetta_rna.<release> \
            -mode score_pdb \
            -pdb <input_pdb> \
            -score:rna_chemical_shift_exp_data <exp_cs_data_file> \
            -params_file <input_param_file> \
            -analytic_etable_evaluation false

3.  Refine (minimize) and score the tandem GA:AG_mismatch NMR PDB structure under the CS-ROSETTA-RNA all-atom energy function:

        ~/Rosetta/rosetta_git/Rosetta/main/source/bin/cs_rosetta_rna.graphics.macosgccrelease \
            -mode minimize_pdb \
            -pdb rosetta_inputs/GA-AG_mismatch/1MIS_NMR.pdb \
            -score:rna_chemical_shift_exp_data rosetta_inputs/GA-AG_mismatch/1MIS_exp_1H_chem_shifts.str \
            -params_file rosetta_inputs/GA-AG_mismatch/1MIS_params \
            -analytic_etable_evaluation false

4.  Score (but not refine) the tandem GA:AG_mismatch NMR PDB structure under the hybrid CS-ROSETTA-RNA all-atom energy function:

        ~/Rosetta/rosetta_git/Rosetta/main/source/bin/cs_rosetta_rna.graphics.macosgccrelease \
            -mode score_pdb \
            -pdb rosetta_inputs/GA-AG_mismatch/1MIS_NMR.pdb \
            -score:rna_chemical_shift_exp_data rosetta_inputs/GA-AG_mismatch/1MIS_exp_1H_chem_shifts.str \
            -params_file rosetta_inputs/GA-AG_mismatch/1MIS_params \
            -analytic_etable_evaluation false


5.  Refine (minimize) and score the UAAC tetraloop NMR PDB structure under the CS-ROSETTA-RNA all-atom energy function:

        ~/Rosetta/rosetta_git/Rosetta/main/source/bin/cs_rosetta_rna.graphics.macosgccrelease \
            -mode score_pdb \
            -pdb rosetta_inputs/UAAC_loop/4A4R_NMR.pdb \
            -score:rna_chemical_shift_exp_data rosetta_inputs/UAAC_loop/4A4R_exp_1H_chem_shifts.str \
            -params_file rosetta_inputs/UAAC_loop/4A4R_params \
            -analytic_etable_evaluation false


6. Score (but not refine) the UAAC tetraloop NMR PDB structure under the hybrid CS-ROSETTA-RNA all-atom energy function:

        ~/Rosetta/rosetta_git/Rosetta/main/source/bin/cs_rosetta_rna.graphics.macosgccrelease \
            -mode minimize_pdb \
            -pdb rosetta_inputs/UAAC_loop/4A4R_NMR.pdb \
            -score:rna_chemical_shift_exp_data rosetta_inputs/UAAC_loop/4A4R_exp_1H_chem_shifts.str \
            -params_file rosetta_inputs/UAAC_loop/4A4R_params \
            -analytic_etable_evaluation false

Optional arguments
------------------

*   `-score::rna_chemical_shift_H5_prime_mode MODE`

    Specify how to handle assignment of the diastereotopic H5' and H5'' proton pair.
    Valid modes:

    * LEAST_SQUARE_IGNORE_DUPLICATES (default)

      In this mode, the assignments of H5' and H5'' protons will be based on which values give better agreement between the experimental and back-calculated chemical shifts.
      Uses this mode, if the experimental non-exchangeable 1H chemical shift are not unambiguously assigned.

    * UNIQUE

      In this mode, the assignments H5' and H5'' proton will be used "as is" and the cs_rosetta_rna will not attempt to switch the H5' and H5'' assignments.
      Use this mode only if the the experimental non-exchangeable 1H chemical shift data have unambiguous assignments of the diastereotopic 1H5´ and 2H5´ protons.
      Note that this is uncommon.

Outputs
-------

1.  A breakdown of the hybrid CS-ROSETTA-RNA all-atom energy terms,
    e.g:

        ------------------------------------------------------------
         Scores                       Weight   Raw Score Wghtd.Score
        ------------------------------------------------------------
         fa_atr                       0.230    -125.447     -28.853
         fa_rep                       0.120       8.314       0.998
         fa_intra_rep                 0.003      81.488       0.236
         fa_intra_RNA_base_phos_atr   0.230       0.000       0.000
         fa_intra_RNA_base_phos_rep   0.120       0.000       0.000
         lk_nonpolar                  0.320       2.123       0.679
         lk_nonpolar_intra_RNA        0.320       3.768       1.206
         fa_elec_rna_phos_phos        1.050      -0.074      -0.078
         ch_bond                      0.420     -30.523     -12.820
         rna_torsion                  2.900       2.721       7.892
         rna_sugar_close              0.700       3.171       2.220
         fa_stack                     0.125    -199.844     -24.981
         geom_sol_fast                0.620      56.483      35.020
         geom_sol_fast_intra_RNA      0.620       1.978       1.226
         hbond_sr_bb_sc               0.620       0.000       0.000
         hbond_lr_bb_sc               2.400       0.000       0.000
         hbond_sc                     2.400     -20.116     -48.279
         hbond_intra                  2.400       0.000       0.000
         atom_pair_constraint         1.000       0.000       0.000
         angle_constraint             1.000       0.000       0.000
         rna_bulge                    0.450       0.000       0.000
         rna_chem_shift               4.000       1.232       4.928
         linear_chainbreak            5.000       0.009       0.047
        -----------------------------------------------------------
         Total weighted score:                              -60.558

2.  The total hybrid CS-ROSETTA-RNA all-atom energy, e.g:

        hybrid_CS-ROSETTA-RNA_all-atom energy: -60.5579

3.  The chemical shift RMSD, e.g:

        chem_shift_RMSD: 0.143299

    The chem_shift_RMSD (in ppm unit) is the root-mean-deviation between
    the 'back-calculated' and the experimental 1H chemical shift. A low
    chem_shit_RMSD indicates that the RNA 3D structure agrees well with
    the experimental 1H chemical shift data

4.  The RNA PDB structure after refinement under the hybrid CS-ROSETTA-RNA all-atom energy function (if -mode minimize_pdb).
    The refined PDB is outputted to the run directory under the filename: `<in_pdb_basename>_out`.

Best Practices
--------------

Figure 1: Breakdown of the secondary structure of the tandem GA:AG mismatch
internal loop:

                           1    6                              1
                        5'-CGGACG-3'                        5'-CG
    Entire structure:      ||**||               H1 helix:      ||
                        3'-GCAGGC-5'                        3'-GC
                           12   7                              12

                                                                    6
                             GA                                    CG-3'
    2x2 mismatch:            **                 H2 helix:          ||
                             AG                                    GC-5'
                                                                    7

* How many canonical base-pairs should be included at each helical boundary?

  2 base-pairs should be included at each helical boundary (for rationale,
  see ref. [2]).

  For example, in the case of the tandem GA:AG mismatch internal loop,
  the structure consists of the a 2 base-pairs H1 helix, a 2x2 mismatch,
  and a 2 base-pairs H2 helix.

* Which atoms' chemical shift data should be included?

  The chemical shift data of all non-exchangeable proton should be 
  included in the chemical shift data file.

  The non-exchangeable protons consist of the H1', H2', H3', H4', H5' and
  H5'' ribose protons, and the H2, H5, H6 and H8 base protons.

  Data lines belonging to other atom types will be ignored.

* Which nucleotides' chemical shift data should be included?

  The chemical shift data of all nucleotides EXCEPT those that are right
  at 5' and 3' edges should be included in the chemical shift data file.

  For example, in the case of the tandem GA:AG mismatch internal loop,
  the chemical shift data of all nucleotides except C1, G6, C7 and G12
  should be included.

* How to prepare the parameters file.

  1. Add a "OBLIGATE  PAIR" line for each helical base-pair located
     right at the 5' and 3' edges of the structure. 

     In the case of the tandem GA:AG mismatch internal loop, the OBLIGATE
     PAIRS are "C1-G12" and "G6-C7":

          OBLIGATE   PAIR 1 12 H H A
          OBLIGATE   PAIR 6 7 H H A

     Note that the cs_rosetta_rna app will refine (minimize) ALL nucleotides
     EXCEPT nucleotides that are specified as "OBLIGATE  PAIR", which will
     be be kept static.

  2. Add "ALLOW_INSERT" lines to include all non-canonical loop
     nucleotides position:

     In the case of the tandem GA:AG mismatch internal loop, the
     "ALLOW_INSERT" nucleotide positions are G3, A4, G9 and A10:

          ALLOW_INSERT 3 4
          ALLOW_INSERT 9 10

  3. Add "CUTPOINT_CLOSED" line to include the position intermediately 5'
     of the first non-canonical loop nucleotide position.

     In the case of the tandem GA:AG mismatch internal loop, the
     first non-canonical loop nucleotide position is G3. The
     "CUTPOINT_CLOSED" position is the position intermediately 5' of
     G3, which is G2:

          CUTPOINT_CLOSED 2

     Note that if the "CUTPOINT_CLOSED" line was not included in the 
     parameter line, the cs_rosetta_rna app will still be able run by
     selecting a random loop position as the cutpoint_closed position.
     However, it is recommended that the "CUTPOINT_CLOSED" line be
     explictly included to prevents this random selection.

  4. Add "CUTPOINT_OPEN" line for all position intermediately 5' of 
     chain-breaks.

     In the case of the tandem GA:AG mismatch internal loop, there is
     chain-break between G6 and C7. The "CUTPOINT_OPEN" is the position
     intermediately 5' of the chain-break which is G6:

          CUTPOINT_OPEN 6

  Adding all the above parameter lines together, we get the parameter file
  for the tandem GA:AG mismatch ("rosetta_inputs/GA-AG_mismatch/1MIS_params"):

      OBLIGATE   PAIR 1 12 H H A
      OBLIGATE   PAIR 6 7 H H A
      ALLOW_INSERT 3 4
      ALLOW_INSERT 9 10
      CUTPOINT_CLOSED 2
      CUTPOINT_OPEN 6

  Finally, the cs_rosetta_rna app can also run WITHOUT an input parameter
  file, although this is not recommended. For this case, a simple
  fold-tree with not chain-break/cutpoints will be used and all
  nucleotides will be refined (minimized).

* How to specify the chemical shift data for the diastereotopic H5' and H5'' proton pairs.

  * If two chemical shift data points are measured for the diastereotopic H5' and  H5'' protons pair and unambiguous assignment is possible, then include correct the unambiguous assignment in the data lines, e.g.:

          1  1  1 G H5'  H   4.180 . . 
          2  1  1 G H5'' H   4.540 . .  

    In this case, please explicitly include the command line option:
      -score:rna_chemical_shift_H5_prime_mode UNIQUE

  * If two chemical shift data points are measured for the diastereotopic H5' and  H5'' protons pair BUT unambiguously assignment is not possible, then include either of the two possible assignments in the data lines, e.g.:

          1  1  1 G H5'  H   4.180 . . 
          2  1  1 G H5'' H   4.540 . .  

                    OR

          1  1  1 G H5'  H   4.540 . . 
          2  1  1 G H5'' H   4.180 . .  

    The cs_rosetta_rna app with automatically select the assignments which leads to better agreement between the experimental and back-calculated chemical shift.

  * If only one chemical shift data point is measured for the diastereotopic H5' and  H5'' proton pair AND unambiguous assignment is not possible, then please still include two chemical shift data lines (with same cs-value), one for each proton, e.g.:

          1  1  1 G H5'  H   4.180 . . 
          2  1  1 G H5'' H   4.180 . . 

References
----------

[1] Das, R., Karanicolas, J. & Baker, D. Nat Methods 7, 291-294 (2010).

[2] Sripakdeeving, P. et al. Nature Methods 11, 413–416 (2014).

[3] Sripakdeevong, P., Kladwang, W. & Das, R. Proc Natl Acad Sci U S A 108, 
    20573-20578 (2011).
# RNA mutation

# Author
Rhiju Das, rhiju@stanford.edu

# Brief Description

Quickly mutate an RNA, from command line.

# Abstract

Easy preparation of mutated versions of an RNA. Not optimized currently to remember residue numbers, etc., but that should be easy to fix up.

# Running

```
rna_thread -s rosetta_inputs/1zih_RNA.pdb  -seq gggcgcgagccu -o 1zih_A7G.pdb 
```

Changes the seventh residue to g. (Original sequence was gggcgcaagccu.)

See also demo for rna_thread.

Prediction of peptide binding using FlexPepBind
-----------------------------------------------
Predicting the structure of a peptide-protein interaction (and an interaction in general) does not guarantee that the peptide indeed binds to the receptor.
However, an accurate model of a peptide-protein interaction can be used to evaluate the binding ability of a given peptide sequence to a given receptor, e.g. relative to other peptides. We have found that this works particularly well when the characteristic binding features of an interactions are defined and reinforced during modeling: sequences for which low-energy models can be generated that fullfill these constraints are assumed to be binders.

We have developed FlexPepBind - a protocol that uses the FlexPepDock framework to model the structure of a range of different peptide-receptor complexes and identifies binders / non-binders.

NOTE: FlexPepBind is NOT a general protocol - it needs to be calibrated for each system based on a small set of experimentally characterized substrates/non-substrates. Once calibrated, new substrates can be identified on a large scale.

This demo illustrates how to run the FlexPepBind protocol to predict peptide binding for two systems: (1) substrates that are deacetylated by Histone Deacetylase 8 (HDAC8), and (2) substrates that are farnesylated (a farnesyl moiety is attached to their c-terminus) by Farnesyl Transferase (FTase).

Protocol overview
-----------------
The FlexPepBind protocol predicts peptide binding by modeling different peptide sequences onto a template peptide-receptor complex. Given such a peptide-receptor complex, a set of peptide sequences are threaded onto the template peptide backbone and further optimized using FlexPepDock protocols.

System-specific calibration: Depending on the system, either simple minimization (of the full peptide conformation and rigid body orientation, and the receptor side chains) or extensive optmization using the FlexPepDock refinement protocol might be needed. Also, the optimal measure to rank the different peptides has to be determined (i.e., a score that highlights the interface energy; see below).

In this demo we show how to use the minimization only protocol to predict binders and non-binders in a set of peptide sequences.

Running the FlexPepBind protocol
--------------------------------
The two main stages of FlexPepBind are:

 1. Thread the query peptide sequence onto the template peptide backbone (using Rosetta fixbb design).
 2. Optimize the threaded complex structure. We use FlexPepDock protocols (either FlexPepDock refinement or minimization only) for this stage.

To run on the specific systems provided here, go to the relevant directory, and:

 a. Run the fpbind_run.sh script (located in the scripts directory). This will minimize the peptide-protein complexes after threading the peptide sequence (listed in the input_files/peptide.list file) one by one onto the template.
 b. After the run finished, run fpbind_analysis.sh (located in the scripts directory). It will extract the relevant scores of the minimized structures & save it in the score_analysis/ directory.

Scores currently considered are:
 1. Interface score - the energy contributed at the interface (I_sc)
 2. peptide score - the energy contributed by the internal energy of the peptide + the interface score (pep_sc)
 3. peptide score without reference energy - the same peptide score, but without the amino acid dependent energy terms (Eaa) that were optimized to generate designs with natural amino acid content (pepe_sc_noref; we found that removing this term can significantly improve distinction of substrates from non-substrates in certain systems)
 4. Reweighted score: = sum (total_score + peptide_score + interface score) (reweighted_sc)

You need to change the paths of the Rosetta executables and database directories in the fpbind_run script.

 ROSETTA_BIN="rosetta/main/source/bin"

 ROSETTA_DB="rosetta/main/database/"

Example run:

 $ cd HDAC8/
 $ bash ../scripts/fpbind_run.sh
 $ bash ../scripts/fpbind_analysis.sh

The output is the score with different scoring terms:
For example, in score_analysis/I_sc

 $ paste input_files/peptide.list score_analysis/I_sc

 GYKFGC	-16.757

 GFKWGC	-17.108

 GFKFGC	-16.352

 GMKDGC	-13.870

 GDKDGC	-13.179

 GQKIGC	-13.408

The above lines show the interface scores for 6 peptides (the top 3 are HDAC8 substrates & the bottom 3 are not HDAC8 substrates).


Further information
-------------------
London N., Lamphear C.L., Hougland J.L., Fierke C.A., and Schueler-Furman O. (2011).Identification of a novel class of farnesylation targets by structure-based modeling of binding specificity. PLoS Comput Biol, 7: e1002170.

London N., Gulla S., Keating A.E., and Schueler-Furman O. (2012). In silico and in vitro elucidation of BH3 binding specificity toward Bcl-2. Biochemistry, 51: 5841-50.

Alam N., Zimmerman L., Wolfson N.A., Joseph C.G., Fierke C.A., and Schueler-Furman O. (2016). Structure-Based Identification of HDAC8 Non-histone Substrates. Structure, 24: 458-68.
This file contains information specific for FTase FlexPepBind. For general information about FlexPepBind, see README file in the upper level directory.

Input files needed to run the FlexPepBind protocol on FTase:
------------------------------------------------------------

1tn6.pdb      : The Ftase-peptide complex structure
template.pdb  : The template structure created from 1tn6
substrates    : list of Ftase peptide substrates
nonsubstrates : list of FTase peptide non-substrates
peptide.list  : list of peptide sequence to the the FlexPepBind protocol on. It contains both the substrates and the non-substrates
resfile       : threading instruction (see more on resfile at https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d1/d97/resfiles.html)
fixbb_flags   : list of flags for running peptdide threading
peptide.list  : list of peptide sequence on which to run the FlexPepBind protocol. It contains both the substrates and the non-substrates
resfile       : threading instructions (see more on resfile at https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d1/d97/resfiles.html)
minimization_flags : file containing list of flags for running FlexPepDock minimization
fpp.params    : param file for the small non-peptidic molecule bound to FTase at the active site (see more on how to create param files for non-protein molecules here: https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/df/de9/preparing_ligands.html)
constraints.cst : list of constraints used during the minimization


Generation of the template structure (template.pdb) for FTase FlexPepBind:
--------------------------------------------------------------------------
Starting with the complex structure 1tn6, a template pdb was created following the steps mentioned below:
a. The input pdb was cleaned by removing solvent molecules.
b. A param file for FPP was created.
c. The receptor structure was prepacked to remove internal clashes. (see the demo refinement_of_protein_peptide_complex_using_FlexPepDock/ for more information on how to run prepacking). The FPP molecule was read in using the flag -extra_res_fa FPP.params.
d. The prepacked structure was used for running FlexPepBind.


How to run the FTase FlexPepBind protocol:
------------------------------------------

Run the FTase FlexPepBind protocol as:
$ bash ../scripts/fpbind_run.sh
$ bash ../scripts/fpbind_analysis.sh

The scoring term that worked best for FTase FlexPepBind is pep_sc_noref. The output is in score_analysis/pep_sc_noref

$ paste input_files/peptide.list score_analysis/pep_sc_noref

CSII	-4.9

CLIT	-3.2

CFLS	-2.2

CKKP	 7.6

CTKR	10.7

CSIP	 5.4

The pep_sc_noref scores are listed for 6 peptides (the top 3 are FTase substrates & the bottom 3 are not FTase substrates; it might differ depending on the rosetta version and scoring function used to run the protocol)
This file contains information specific for HDAC8 FlexPepBind. For general information about FlexPepBind, see README file in the upper level directory.

Input files needed to run the FlexPepBind protocol on HDAC8:
------------------------------------------------------------
template.pdb  : The template structure created from 2v5w
substrates    : list of HDAC8 peptide substrates
nonsubstrates : list of HDAC8 peptide non-substrates
peptide.list  : list of peptide sequence to the the FlexPepBind protocol on. It contains both the substrates and the non-substrates
activity      : percent deacetylation of the peptides listed in peptide.list file
resfile       : threading instruction (see more on resfile at https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d1/d97/resfiles.html)
fixbb_flags   : file containing list of flags for peptide threading
minimization_flags : file containing list of flags for running FlexPepDock minimization
constraints.cst : list of constraints used during the minimization

Generation of the template structure (template.pdb) for HDAC8 FlexPepBind:
--------------------------------------------------------------------------
The solved structure of a complex of HDAC8 bound to a substrate (PDB id 2v5w) contains only a short peptide-coumarine ligand. Therefore, we needed to generate a structure to serve as template for hexamer-HDAC8 interactions. It was generated following the steps mentioned below:
a. The input pdb was cleaned by removing solvent molecules.
b. Only a monomer bound to a peptide was copied from the solved dimeric structure.
c. The receptor structure was prepacked to remove internal clashes (see the demo refinement_of_protein_peptide_complex_using_FlexPepDock/ for more information on how to run prepacking).
d. The peptide substrates GYKacFGC was threaded onto the template peptide backbone (missing residues were added in extended conformation), and optimize extensively using the FlexPepDock refinement protocol (see the demo refinement_of_protein_peptide_complex_using_FlexPepDock/ for more information on how to run FlexPepDock refinement), and the top-scoring optmized structure (according to interface score; I_sc) was used as the template.

The underlying assumption was that a structural model of a strong substrate (GYKacFGC) would be a good representative of the preferred binding conformation of HDAC8 substrates. Therefore, we used it as a template to evaluate binding ability of other peptides.

How to run the HDAC8 FlexPepBind protocol:
------------------------------------------
$ bash ../scripts/fpbind_run.sh
$ bash ../scripts/fpbind_analysis.sh

The output contains scores for each peptide, according to different scoring terms. For HDAC8 FlexPepBind, we found that the interface score (I_sc) provides the best ranking.

$ paste input_files/peptide.list score_analysis/I_sc
GYKFGC	-16.757
GFKWGC	-17.108
GFKFGC	-16.352
GMKDGC	-13.870
GDKDGC	-13.179
GQKIGC	-13.408

The Interface scores are given for 6 peptides (the top 3 are HDAC8 substrates & the bottom 3 are not HDAC8 substrates; the score might vary depending on the rosetta version and socirng function being used).
This is an example for Rosetta 3D RNA modeling protocol used by the Das lab.

It should be used with the following tutorial in google doc.

https://docs.google.com/document/d/18MaRXkq1a_J6BQY9CKkF2swc7MZnFy63-421I5CzzFc/edit?usp=sharing

Analyzing Structure Quality
===========================
One of the simplest initial steps to evaluate structural quality is to score the structure using the Rosetta score function and to evaluate the results.

Scoring structures
------------------
To score structures, use the score application. This application can accept any input structure that Rosetta recognizes.

    /path/to/rosetta/bin/score.macosgccrelease -database /path/to/minirosetta_database/ -s ../starting_files/1ubq.pdb.gz -in:file:fullatom -ignore_unrecognized_res -out:output

Description of command-line options
-----------------------------------
* `-database`: Location of the Rosetta database on your system
* `-s`: One or more input strucutre to score, in PDB format.  To use a silent file as input, use `-in:file:silent` instead of `-s` (optionally with `-in:file:tags`).
* `-in:file:fullatom`: Score the structure using the full atom energy function.
* `-out:file:score`: Name of the summary scorefile.
* `-out:output`: Force production of the scored PDBs.

Description of score terms
--------------------------
Rosetta scoring values are reported in REU (Rosetta Energy Units). One study reported a conversion factor of 0.58 from REU to kcal/mol (Kellogg, Leaver-Fay, Baker; Proteins 2011). This conversion may not apply for other applications. Lower (more negative) values indicate a more favorable conformation.

* `fa_atr`:         Lennard-Jones attractive
* `fa_rep`:         Lennard-Jones repulsive
* `fa_sol`:         Lazaridis-Karplus solvation energy
* `fa_intra_rep`:   Lennard-Jones repulsive between atoms in the same residue
* `fa_pair`:        statistics based pair term, favors salt bridges
* `fa_plane`:       pi-pi interaction between aromatic groups, by default = 0
* `fa_dun`:         internal energy of sidechain rotamers as derived from Dunbrack's statistics
* `ref`:            reference energy for each amino acid
* `hbond_lr_bb`:    backbone-backbone hbonds distant in primary sequence
* `hbond_sr_bb`:    backbone-backbone hbonds close in primary sequence
* `hbond_bb_sc`:    sidechain-backbone hydrogen bond energy
* `hbond_sc`:       sidechain-sidechain hydrogen bond energy
* `p_aa_pp`:        Probability of amino acid at phipsi
* `dslf_ss_dst`:    distance score in current disulfide
* `dslf_cs_ang`:    csangles score in current disulfide
* `dslf_ss_dih`:    dihedral score in current disulfide
* `dslf_ca_dih`:    ca dihedral score in current disulfide
* `pro_close`:      proline ring closure energy
* `rama`:           ramachandran preferences
* `omega`:          omega dihedral in the backbone
* `total`:          total energy for structure/residue

Evaluating results
------------------
The score application outputs one summary file (`scorefile.sc`), containing summary information for all input files, and then one scored PDB file per input structure.
In the summary file, there is one line per input structure, listing the total energy from each term over the whole conformation.
At the end of each of the output PDB structures is a table with a per-residue breakdown of score terms.
The table is useful for identifying poorly scoring residues or analyzing the contribution of different score terms to the total for a particular residue or residues.
Design with Phosphoresidues
===========================

Phosphoresidues are very important in biology.  Many signaling pathways depend 
on the fact that phosphorylation of residues like serine, threonine, and 
tyrosine can introduce dramatic allosteric changes.  This demo shows how to 
use phosphoresidues in Rosetta.

Make sure the phosphoresidue is read properly
---------------------------------------------

We create the parameter file for the phosphoresidue as follows:

1.  Extract a phosphotyrosine (PTR) residue from a PDB file

    We use one model of the NMR ensembles of the 2lct.pdb by deleting other 
    models and keep the first one (3lct.pdb refers to first model). From the 
    "3lct.pdb" file, extract the part containing the atoms of the 
    phosphoresidues 342 and 346 to generate "PT1.pdb" and "PT2.pdb". You can 
    use the sample command line

        grep PTR 3lct.pdb | grep 342 | grep HETATM > PT1.pdb
        grep PTR 3lct.pdb | grep 346 | grep HETATM > PT2.pdb

2.  Create a MOLFILE for PTR

    Either use avogadro(free chemical software) or go to the adress: 
    www.molecular-networks.com/online_demos/convert_demo to generate an 
    "PT1.mdl" file from the "PT1.pdb" with the original geometry of the 
    phosphoresidue.

3.  Create a Rosetta PARAMS file PTR

    Copy the "PTR.mdl" file (molfile format) into your working directory.  Then 
    run:

        mini/src/python/apps/public/molfile_to_params.py PT1.mdl -n PT1
        mini/src/python/apps/public/molfile_to_params.py PT2.mdl -n PT2

    This script creates PT1.params and PT1_0001.pdb. You'll need to delete the 
    older PTR residues in 3lct.pdb and replace them with the co-ordinates for 
    generated new residues PT1_0001.pdb and PT2_0001.pdb (4lct.pdb refers to 
    3lct.pdb with new phoshoresidues). Also, change the "HETATM" tags for 
    phoshoresidues in the PDB file to "ATOM".

Design around the phosphoresidue
--------------------------------

1.  Using pyMol, identify the neighbouring residues, within the desired radius 
    (5 A in our case, PRT.resfile)

2.  Write the resfile. Use the option ALLAA next to the sequence position to 
    design the selected residue to all of posibble amino acids.

3.  Run the following command-line:

        ~/mini/bin/fixbb.default.macosgccrelease -s 4lct.pdb -database ~/minirosetta_database -extra_res_fa PT1.params PT2.params -resfile PTR.resfile

    This will generate a single pdb file with designed residues around the 
    phosphoresidues. The sample outputs have been copied to /outputs. If you 
    are using this protocol for a real design application, use the flag 
    `-nstruct 1000` or `-nstruct 10,000` to generate enough designs to give 
    meaningful results.



# StepWise Monte Carlo (examples for RNA)

# Author
Rhiju Das, rhiju@stanford.edu

# Brief Description

Steps to build a model of a complex RNA fold

# Abstract

This code allows build-up of three-dimensional de novo models of RNAs of sizes up to ~300 nts, given secondary structure and experimental constraints. It can be carried out reasonably automatically, but human curation of submodels along the build-up path may improve accuracy. A fully automated pipeline is also in preparation.

# More docs
This documentation (and more) are available in the on-line docs at:

https://www.rosettacommons.org/docs/latest/rna-denovo-setup.html

# Example Rosetta Command Lines

## Make helices

Example in step1_helix/

```
rna_helix.py  -o H2.pdb -seq cc gg -resnum 14-15 39-40
replace_chain_inplace.py  H2.pdb 
```

## Use threading to build sub-pieces

Example in step2_thread/

In the problem above, there is a piece which is a well-recognized motif, 
the UUCG apical loop. Let's model it by threading from an exemplar
of the motif from the crystallographic database. There is one here:

Download 1f7y.pdb from 
http://pdb.org/pdb/explore/explore.do?structureId=1f7y

Slice out the motif of interest:
```
pdbslice.py  1f7y.pdb  -subset B:31-38 uucg_
```

Thread it into our actual sequence:
```
rna_thread -s uucg_1f7y.pdb  -seq ccuucggg -o uucg_1f7y_thread.pdb
```

Let's get the numbering to match our actual test case:
```
renumber_pdb_in_place.py uucg_1f7y_thread.pdb 24-31
```

Done!

## Build models of a sub-piece denovo

In step3_farfar/, we will see how to setup the Rosetta job for motifs between H2 and H4, using our starting H2 and H4 helices as fixed boundary conditions. 

There is currently a wrapper script that sets up the job for the rna_denovo executable, which actually runs fragment assembly of RNA with full atom refinement (FARFAR) is not yet equipped to map numbers from our full modeling problem into the subproblem. We have to create it a little sub-problem and map all the residue numberings into the local problem.

There's a file called README_SETUP which has the wrapper command to set up the job. For completeness, the command there is:

```
rna_denovo_setup.py -fasta RNAPZ11.fasta \
    -secstruct_file RNAPZ11_OPEN.secstruct \
    -working_res 14-25 30-40 \
    -s H2.pdb H4.pdb \
    -fixed_stems \
    -tag H2H3H4_run1b_openH3_SOLUTION1 \
    -native example1.pdb 
```

You don't need to supply a native if you don't have it -- just useful
to compute RMSDs as a reference.

You can run the command by typing:

```
 source README_SETUP
```

Then try this:

```
 source README_FARFAR
```

Example output after a couple of structures is in example_output/.

[You should probably do a full cluster run -- some tools are available for
 condor, qsub, slurm queueing systems, documented here:

https://www.rosettacommons.org/docs/latest/RNA-tools.html

]

Extract 10 lowest energy models:

```
extract_lowenergy_decoys.py H2H3H4_run1b_openH3_SOLUTION1.out 10
```

Inspect in pymol.
(For an automated workflow, you can also cluster these runs and just carry forward the top 5 clusters.)

## Graft together models for the full-length RNA

Example in step4_graft/:

These were threading and FARFAR solutions that we liked for each submotif -- now we can graft:

```
rna_graft -s H2H3H4_run1b_openH3_SOLUTION1.pdb  uucg_1f7y_thread.pdb  H1H2_run2_SOLUTION1.pdb -o full_graft.pdb
```

Done! 

# StepWise Monte Carlo (example for protein loop)

## Author
Rhiju Das, rhiju@stanford.edu

## Brief Description

Solve structure of a protein loop

## Abstract

Ab initio and comparative modeling of biopolymers (RNA, protein, protein/RNA) often involves solving well-defined small puzzles (4 to 20 residues), like RNA aptamers, RNA tertiary contacts, and RNA/protein interactions. If these problems have torsional combinations that have not been seen previously or are not captured by coarse-grained potentials, most Rosetta approaches will fail to recover their structures.  This app implements a stepwise ansatz, originally developed as a 'stepwise assembly' enumeration that was not reliant on fragments or coarse-grained modeling stages, but was computationally expensive. The new mode is a stepwise monte carlo, a stochastic version of stepwise assembly. 


## Running

### Example Rosetta Command Line

```
stepwise -s rosetta_inputs/noloop_mini_1alc_H.pdb -fasta rosetta_inputs/mini_1alc.fasta -native rosetta_inputs/mini_1alc.pdb -score:weights stepwise/protein/protein_res_level_energy.wts -silent swm_rebuild.out -from_scratch_frequency 0.0 -allow_split_off false -cycles 200 -nstruct 20
```

Most of the simulation may be spent flickering bits of secondary structure -- in the future, we will probably setup some precomputation of these bits so that computation can be focused on build up of the complete mini-protein structure.



Shortening loops in Rosetta
===========================

To shorten loops in Rosetta, you will edit the PDB of the original loop to 
delete the undesired residues, build a loop and fragment file for the newly 
short protein, and run the loop through loop modeling with build_initial 
active.

Manually removing the loop from the PDB
----------------------------------------

Take your input PDB and delete the residue lines related to the residue you are 
deleting.  Here we are shortening 1FNA.pdb by deleting the ALA with the residue 
number 83. For the deletion, open the PDB file with your favorite text editor 
and remove the ATOM entry lines of residue 83. In addition remove also all non 
ATOM entry lines, i.e. lines not starting with the term ATOM.

Preparing fragments and loop files
----------------------------------

1FNA_del.pdb is the primary Rosetta input. 

We also need a fragments file and a loop file.  The fragments file is best 
created via the [[Robetta server|http://robetta.bakerlab.org/fragmentqueue.jsp]] 
or using the fragment insertion tutorial.  You'll need a FASTA file of your 
protein to generate fragments; don't forget to delete the deleted-residue from 
your FASTA file and also note that their is shift in the amino acid number 
between the fasta file and the PDB file. You can find the FASTA file of the PDB 
protein in starting_files.

The loop file format is:

    LOOP START STOP CUT SAMPLE EXTENDED

In our case, the loop start and stop are 71 (T76) and 81 (87P), respectively. 
Note that this file is in Rosetta numbering, not in PDB numbering.  The 
cutpoint is 77, the position before the removed residue (A83 was removed; 77 is 
P82).

Rebuilding the loop with Rosetta
--------------------------------

The only executable we'll need here is loop modeling. Briefly, we will run 
loop modeling with the build_initial mode active. This option triggers Rosetta 
to rebuild the loop in a closed state. To do the bulk of the loop remodeling, 
we will choose KIC remodeling.  CCD would work as well.

    path/to/loopmodel.[platform][compiler][mode] @rosetta_inputs/options -database path/to/rosetta_database

An options file has been provided (rosetta_inputs/options), annotated with a 
description of what each flag is doing.

* `-database`  
  Specify the path to the Rosetta database, required for any Rosetta 
  simulation.

* `@rosetta_inputs/options`  
  File holding all rosetta commandline flags. See section "Option file" below.

* `-in:file:fullatom`  
  Necessary for pretty much all loop modeling runs to read in PDBs properly.

* `-loops:input_pdb rosetta_inputs/1FNA_del.pdb`  
  Path to input pdb.

* `-loops:loop_file rosetta_inputs/loop_file`  
  Path to loops file.

* `-loops:frag_sizes 9 3 1`  
  what sizes are the fragments?  9 and 3 are traditional.  The flag seems to 
  require a third argument, but you can pass no fragments in that size.

* `-loops:frag_files rosetta_inputs/aa1FNA_09_05.200_v1_3 rosetta_inputs/aa1FNA_03_05.200_v1_3 none`  
  Paths to the fragments in the same vein as previous - none for 1mer 
  fragments.

* `-loops::build_initial`  
  This flag triggers build initial mode, which fixes the broken loop before 
  re-solving it.

* `-loops:remodel perturb_kic -loops:refine refine_kic`  
  These flags specify KIC loop modeling to remodel the loops.

* `-loops:remodel perturb_ccd -loops:refine refine_ccd`  
  These flags could be used instead of those above to specify CCD remodeling.

* `-out:path sample_output`  
  Output directory.

* `-out:prefix 1FNA_del_`  
  Prefix for output.  Would not be necessary if someone would rewrite loop 
  modeling to use jd2.

* `-nstruct 1`  
  This option controls how many output structures you get; larger is better!  1 
  is used here because the tutorial can only take so long; in production you'd 
  use 10000 or more.

Interpreting the results
------------------------

You will get PDB files with energy scores at the bottom of the file as your 
results. To find the best structures, sort the PDB by their total score (using 
script sort_by_score), and manually examine the top 5% (or top 100, or whatever 
you have time for) in PyMOL or another viewer.  Using a combination of total 
score and your protein intuition (this is an art, not a science), pick which 
you think is best.

Scientifically, this tutorial is underpowered - at best you could use these 
results to determine IF a loop can be closed after residues have been deleted. 
To do true loop modeling, you would first close the loop using a quick protocol 
like this, then run a more rigorous loop modeling protocol (covered in the main 
loop modeling documentation)
Refinement of peptide-protein complexes using FlexPepDock refinement
--------------------------------------------------------------------
This demo illustrates how to run FlexPepDock refinement of a peptide-protein complex. The FlexPepDock Refinement protocol is designed to create high-resolution models of complexes between a flexible peptide and a globular protein, with side chains of binding motifs modeled at nearly atomic accuracy. It is intended for cases where an approximate, coarse model of the interaction is available.

Protocol Overview
-----------------
The input to the refinement protocol is a coarse model of the peptide-protein complex in PDB format. The protocol iteratively optimizes the peptide backbone and its rigid-body orientation relative to the receptor protein, including periodic on-the-fly side-chain optimization. The protocol is able to account for a considerable diversity of peptide conformations within a given binding site. However it is important to note that the refinement protocol can only refine models which are close to the correct solution both in terms of Cartesian and dihedral (phi/psi) distance. For the cases where no information regarding the peptide backbone is availavle, we recommend to use the FlexPepDock ab-initio protocol (See the demo 'abinitio_fold_and_dock_of_peptides_using_FlexPepDock').


Running the FlexPepDock refinement protocol:
--------------------------------------------
1. Create an initial complex structure: A coarse model of the peptide-protein interaction can be obtained from a low resolution peptide-protein docking protocol or can be built using a homologous structure. Here in this demo we have provided 1AWR.ex.pdb in which the peptide is in extended conformation. This will serve as the coarse model and our goal will be to refine it to obtain a near-native model. The native structure is 1AWR.pdb. Both 1AWR.ex.pdb and 1AWR.pdb are located in the input directory.

2. Prepack the input model: This step involves the packing of the side-chains in each monomer to remove internal clashes that are not related to inter-molecular interactions. The prepacking guarantees a uniform conformational background in non-interface regions, prior to refinement. The prepack_flags file contains the flags for running the prepacking job. The run_prepack script will run prepacking of the input structure 1AWR.ex.pdb located in the input directory.

You need to change the paths of the Rosetta executables and database directories in the run_prepack script (also for run_refine; see below).

 ROSETTA_BIN="rosetta/main/source/bin"
 ROSETTA_DB="rosetta/main/database/"

After changing the paths run the run_prepack script as:
 $./run_prepack

The output will be a prepacked structure, 1AWR.ex.ppk.pdb located in the input directory; a scorefile named ppk.score.sc and a log file named prepack.log file located in the output directory. This prepacked structure will be used as the input for the refinement step.

3. Refine the prepacked model: This is the main part of the protocol. In this step, the peptide backbone and its rigid-body orientation are optimized relative to the receptor protein using the Monte-Carlo with Minimization approach, including periodic on-the-fly side-chain optimization. An optional low-resolution (centroid) pre-optimization will increase the sampling range and may improve performance further. The file refine_flags contains flags for running the refinement job. The run_refine script will run refinement of the prepacked structure generated in the prepacking step located in the input directory.

After changing the Rosetta related paths run the run_refine script as:
 $./run_refine

The output will be a refined structure (1AWR.ex.ppk_0001.pdb) located in the output directory; a scorefile named refine.score.sc and a log file named refine.log file located in the output directory. This script has to be modified to run on a cluster during a production run.

Further information
-------------------
A detailed documentation on FlexPepDock is available at https://www.rosettacommons.org/docs/latest/application_documentation/docking/flex-pep-dock
Raveh B, London N, Schueler-Furman O. (2010).Sub-angstrom modeling of complexes between flexible peptides and globular proteins. Proteins 78:2029-40.

AbInitio Structure Prediction Using Chemical-Shift Generated Fragments, NOE Distance Restraints and RDC Restraints
==================================================================================================================

Written by Lei Shi.
Nikolas Sgourakis drafted the previous version

---

We will use the chemical shifts to improve the fragments from which Rosetta builds up structures, and the NOEs+RDCs to guide the Rosetta calculations towards the native structure. 

Please see references at:
* rosetta abinitio: Bradley, P et al Science 2005
* chemical shift fragments: Shen Y et al. PNAS 2008;105:4685-4690
* chemical shift+NOE+RDC: Raman S, et al Science 2010

These Rosetta calculation steps are also described separately:
* Sgourakis NG et al JACS,2011,133(16):6288-98:

In this demo, we will use PDB 2JY7, which is a small protein (for demo purpose) and has experimental data deposited. Several scripts are provided in the scripts folder for formattting purposes:

    bmrb2talos.com
    cst_map_toCB.py
    upl2mini.csh
    scores.score.cfg

If you are from David Baker lab, there are scripts available to make setup easier without going through public servers. The following instructions should work just fine without having direct access to any Baker lab cluster.

Running the demo
----------------
1. Create following folders:
    ```
    mkdir starting_inputs
    mkdir rosetta_inputs
    mkdir rosetta_inputs/talos_output
    mkdir rosetta_inputs/pick_cs_fragments
    ```

2. Download protein fasta and experimental data.  
Download fasta from http://www.pdb.org/pdb/explore/explore.do?structureId=2JY7
    ```
    wget http://www.pdb.org/pdb/files/fasta.txt?structureIdList=2JY7 -O starting_inputs/t000_.fasta
    ```
Download chemical shift data from http://www.bmrb.wisc.edu/data_library/summary/index.php?bmrbId=15591
    ```
    wget http://rest.bmrb.wisc.edu/bmrb/NMR-STAR2/15591 -O starting_inputs/raw.cs.bmrb
    ```
Download NOE data from http://restraintsgrid.bmrb.wisc.edu/NRG/MRGridServlet?pdb_id=2JY7&show_blocks=true&min_items=0:
    ```
    #http://restraintsgrid.bmrb.wisc.edu/NRG/MRGridServlet?db_username=wattos1&format=ambi&mrblock_id=434910&pdb_id=2jy7&program=DYANA%2FDIANA&request_type=block&subtype=general+distance&type=distance
    echo "save file as starting_inputs/NOE_data.upl"
    ```
Download RDC data (used for a different demo) from http://restraintsgrid.bmrb.wisc.edu/NRG/MRGridServlet?pdb_id=2JY7&show_blocks=true&min_items=0:
    ```
    #http://restraintsgrid.bmrb.wisc.edu/NRG/MRGridServlet?db_username=wattos1&format=n%2Fa&mrblock_id=23753&pdb_id=2jy7&request_type=archive&subtype=n%2Fa&type=dipolar+coupling
    echo "save file as starting_inputs/nh_xplor.rdc"
    ```

3. Format data for Rosetta use:  
Formatting NOE: Note only residues separated by more than 3 are kept in constraint.
This script `scripts/upl2mini.csh` only works with cyana format NOE.
    ```
    scripts/upl2mini.csh starting_inputs/NOE_data.upl > rosetta_inputs/NOE.cst
    scripts/cst_map_toCB.py rosetta_inputs/NOE.cst > rosetta_inputs/NOE.centroid.cst
    ```
Formmatting chemical shift data for TALOS:
    ```
    scripts/bmrb2talos.com starting_inputs/raw.cs.bmrb > rosetta_inputs/cs.talos
    ```
Formmatting RDC data for Rosetta:
Note that rename HN to H in the rdc file:
    ```
    awk '{print $4,$7,$11,$14,$16}' starting_inputs/nh_xplor.rdc > rosetta_inputs/nh.rdc
    ```

4. Generating talos predictions using http://spin.niddk.nih.gov/bax/nmrserver/talosn/ using rosetta_inputs/cs.talos
Save/copy pred.tab and predSS.tab to rosetta_inputs/talos_output

5. Generate fragment/profile from RobettaServer http://www.robetta.org/fragmentqueue.jsp using starting_inputs/t000_.fasta
Save/copy t000_.checkpoint to rosetta_inputs/

6. Pick fragments using secondary structure profile and chemical shift data
    ```
    Rosetta/main/source/bin/fragment_picker -database Rosetta/main/database/ -in::file::vall Rosetta//tools/fragment_tools/vall.apr24.2008.extended.gz -frags::n_frags 200 -frags::frag_sizes 3 9 -frags::sigmoid_cs_A 2 -frags::sigmoid_cs_B 4 -out::file::frag_prefix rosetta_inputs/pick_cs_fragments/frags.score -frags::describe_fragments rosetta_inputs/pick_cs_fragments/frags.fsc.score -frags::scoring::config scripts/scores.score.cfg -in:file:fasta starting_inputs/t000_.fasta -in:file:checkpoint rosetta_inputs/t000_.checkpoint -in:file:talos_cs rosetta_inputs/cs.talos -frags::ss_pred rosetta_inputs/talos_output/predSS.tab talos -in::file::talos_phi_psi rosetta_inputs/talos_output/pred.tab
    ```

7. Run Rosetta with the fragments made above and use NOEs to guide search in both centroid sampling and full-atom optimization
```
Rosetta/main/source/bin/minirosetta -database Rosetta/main/database/ -in:file:rdc rosetta_inputs/nh.rdc -cst_fa_file rosetta_inputs/NOE.cst -cst_file rosetta_inputs/NOE.centroid.cst -abinitio:stage1_patch scripts/patch_atom_pair_constraint_rdc -abinitio:stage2_patch scripts/patch_atom_pair_constraint_rdc -abinitio:stage3a_patch scripts/patch_atom_pair_constraint_rdc -abinitio:stage3b_patch scripts/patch_atom_pair_constraint_rdc -abinitio:stage4_patch scripts/patch_atom_pair_constraint_rdc -score:patch scripts/patch_atom_pair_constraint_rdc -in:file:fasta starting_inputs/t000_.fasta -file:frag3 rosetta_inputs/pick_cs_fragments/frags.score.200.3mers -file:frag9 rosetta_inputs/pick_cs_fragments/frags.score.200.9mers -nstruct 1 -out:file:silent csrosetta_noe_rdc.out -run:protocol abrelax -abinitio::relax -overwrite
```
If you want to use non-linear-square fitting with Da or R fixed, use:
```
-rdc:fit_method nls -rdc:fixDa 12.6 -rdc:fixR 0.42
```
You can/should adjust the weights of NOE and rdc constraints in scripts/patch_atom_pair_constraint_rdc.
You can provide either `-rdc:fixDa`, `-rdc:fixR` or both of them
Change nstruct to generate desired number of models. Larger is better depending on your available computer time, etc.
Note that the demo in abinitio_w_chemicalshift_only, you can add flags such as:
```
    -abinitio::rg_reweight 0.5 
    -abinitio::rsd_wt_helix 0.5 
    -abinitio::rsd_wt_loop 0.5 
    -disable_co_filter true 
    -abinitio::increase_cycles 10
```
to help sampling in centroid stage.
They are not used probably NOE constraint helps guided the search.

Processing the output
---------------------
1. Extract the low energy models:
    ```
    grep SCORE csrosetta_noe_rdc.out | sort –nk2 | head
    ```
The second column contains the energies of the lowest energy 10 models.
Select as the cutoff the energy on the last line.
You should also use NOE constraint energy and RDC energy as a criteria to select structures.
Example is only provided for total score.

2. This script:
    ```
    cull_silent.pl csrosetta_noe_rdc.out “score < cutoff”
    ```
will produce csrosetta.select.silent which contains the lowest total energy 10 models.

3. extract pdbs from silent files for a given tag in the silent file
    ```
    Rosetta/main/source/bin/extract_pdbs -database Rosetta/main/database/ -in::file::silent csrosetta_noe_rdc.out -in::file:tags S_00000001  
    ```

4. Check convergence by superimposing the ten low energy models in pymol or your favorite molecular graphics

5. Check convergence by clustering the lowest energy models (see clustering demo for instructions)

6. To see how NOE constraints are satisfied by a model:
    ```
    Rosetta/main/source/bin/r_cst_tool.linuxgccrelease -database Rosetta/main/database/ -in:file:s lowscore_1.pdb -cst_file rosetta_inputs/NOE.cst
    r_cst_tool is a pilot program by Oliver Lange in Rosetta/main/source/src/apps/pilot/olli/
    ```

7. To see how NOE constraints are satisfied by a model:
Rosetta/main/source/bin/r_cst_tool.linuxgccrelease -database Rosetta/main/database/ -in:file:s lowscore_1.pdb -cst_file rosetta_inputs/NOE.cst
    ```
    r_cst_tool is a pilot program by Oliver Lange in Rosetta/main/source/src/apps/pilot/olli/
    ```

8. To see how rdc constraints are satisfied by a model
    ```
    Rosetta/main/source/bin/score_jd2.default.linuxgccrelease -database Rosetta/main/database/ -in:file:s lowscore_1.pdb -in::file::rdc rosetta_inputs/nh.rdc -score::patch scripts/patch_atom_pair_constraint_rdc -rdc:fit_method nls 1 -out:file:scorefile rdc_score_for_lowscore_1.sc -out:level 999 -mute core.chemical core.conformation -rdc:print_rdc_values lowscore_1.rdc_corr
    ```
AbInitio Structure Prediction Using Chemical-Shift Generated Fragments
======================================================================

Written by Lei Shi.
Ray Wang drafted the previous version.

---

We will use the chemical shifts to improve the fragments from which Rosetta builds up structures, and the NOEs to guide the Rosetta calculations towards the native structure.

Please see references at:
* Rosetta abinitio: Bradley, P et al Science 2005
* Chemical shift fragments: Shen Y et al. PNAS 2008;105:4685-4690
* Chemical shift+NOE+RDC: Raman S, et al Science 2010

In this demo, we will use PDB 2JY7, which is a small protein (for demo purpose) and has experimental data deposited. Several scripts are provided in the scripts folder for formatting purposes:

    bmrb2talos.com
    cst_map_toCB.py
    upl2mini.csh
    scores.score.cfg

If you are from David Baker lab, there are scripts available to make setup easier without going through public servers. The following instructions should work just fine without having direct access to any Baker lab cluster.

Running the Demo
----------------
1. Create following folders:
    ```
    mkdir starting_inputs
    mkdir rosetta_inputs
    mkdir rosetta_inputs/talos_output
    mkdir rosetta_inputs/pick_cs_fragments
    ```

2. Download protein fasta and experimental data:  
Download fasta from http://www.pdb.org/pdb/explore/explore.do?structureId=2JY7
    ```
    wget http://www.pdb.org/pdb/files/fasta.txt?structureIdList=2JY7 -O starting_inputs/t000_.fasta
    ```
Download chemical shift data from http://www.bmrb.wisc.edu/data_library/summary/index.php?bmrbId=15591
    ```
    wget http://rest.bmrb.wisc.edu/bmrb/NMR-STAR2/15591 -O starting_inputs/raw.cs.bmrb
    ```

3. Format data for Rosetta use:  
Formmatting chemical shift data for TALOS:
    ```
    scripts/bmrb2talos.com starting_inputs/raw.cs.bmrb > rosetta_inputs/cs.talos
    ```

4. Generating talos predictions using http://spin.niddk.nih.gov/bax/nmrserver/talosn/ using `rosetta_inputs/cs.talos`.
Save/copy `pred.tab` and `predSS.tab` to `rosetta_inputs/talos_output`

5. Generate fragment/profile from RobettaServer http://www.robetta.org/fragmentqueue.jsp using `starting_inputs/t000_.fasta`.
Save/copy `t000_.checkpoint` to `rosetta_inputs/`

6. Pick fragments using secondary structure profile and chemical shift data:
    ```
    Rosetta/main/source/bin/fragment_picker -database Rosetta/main/database/ -in::file::vall Rosetta//tools/fragment_tools/vall.apr24.2008.extended.gz -frags::n_frags 200 -frags::frag_sizes 3 9 -frags::sigmoid_cs_A 2 -frags::sigmoid_cs_B 4 -out::file::frag_prefix rosetta_inputs/pick_cs_fragments/frags.score -frags::describe_fragments rosetta_inputs/pick_cs_fragments/frags.fsc.score -frags::scoring::config scripts/scores.score.cfg -in:file:fasta starting_inputs/t000_.fasta -in:file:checkpoint rosetta_inputs/t000_.checkpoint -in:file:talos_cs rosetta_inputs/cs.talos -frags::ss_pred rosetta_inputs/talos_output/predSS.tab talos -in::file::talos_phi_psi rosetta_inputs/talos_output/pred.tab
    ```

7. Run Rosetta with the fragments chemical shift fragments:
    ```
    Rosetta/main/source/bin/AbinitioRelax -database Rosetta/main/database/ -in:file:fasta starting_inputs/t000_.fasta -file:frag3 rosetta_inputs/pick_cs_fragments/frags.score.200.3mers -file:frag9 rosetta_inputs/pick_cs_fragments/frags.score.200.9mers -nstruct 1 -abinitio::increase_cycles 10 -abinitio::relax -score::weights score13_env_hb -abinitio::rg_reweight 0.5 -abinitio::rsd_wt_helix 0.5 -abinitio::rsd_wt_loop 0.5 -disable_co_filter true -out:file:silent csrosetta.out
    ```
Change nstruct to generate desired number of models. Larger is better depending on your available computer time, etc.

Processing the output
---------------------
1. Extract the low energy models:
    ```
    grep SCORE csrosetta.out | sort –nk2 | head
    ```
The second column contains the energies of the lowest energy 10 models.
Select as the cutoff the energy on the last line.

2. This script:
    ```
    cull_silent.pl csrosetta.out “score < cutoff”
    ```
will produce `csrosetta.select.silent` which contains the lowest models below cutoff.
  
3. Extract pdbs from selected silent file
    ```
    Rosetta/main/source/bin/extract_pdbs -database Rosetta/main/database/ -in::file::silent csrosetta.select.silent
    ```

4. Check convergence by superimposing the ten low energy models in pymol or your favorite molecular graphics

5. Check convergence by clustering the lowest energy models (see clustering demo for instructions)
Design with flexible loops
==========================

Method 1: RosettaRemodel
------------------------

Documentation: https://wiki.rosettacommons.org/index.php/Remodel

Requirement:

    svn co https://svn.rosettacommons.org/source/trunk/rosetta_scripts/remodel

1.  Preparing the PDB

    Take the PDB 3k2m.pdb. Select chain B and C. Get rid of HETATM. Do fast 
    relax.

        ~/mini/bin/relax.macosgccrelease -database /Users/rjha/minirosetta_database -s 3k2m_bc.pdb -ignore_unrecognized_res -use_input_sc -constrain_relax_to_start_coords -relax:fast -out:file:renumber_pdb

    Let's say the PDB is `3k2m_bc_relax.pdb`. Delete the Chain ID column.

2.  Generating BluePrint

        ./getBluePrintFromCoords.pl -pdbfile 3k2m_bc_relax.pdb > test.blueprint

    The blueprint file has the information to direct the protocol on which 
    residue to design and remodel.

3.  Editing BluePrint

    The interface on the monobody is Residue 182-187. We decided to rebuild with 
    1 residue insertion starting at 183. The blueprint has to be modified in 
    the following way

        ...
        178 G .
        179 E .
        180 D .
        181 S .
        182 A L
        183 G L
        0 x L
        184 Y E
        185 M E
        186 F E
        187 M E
        188 Y .
        189 S .
        190 P .
        191 I .
        ...

    In the above example, "0 x L" will mean eXtension and the secondary 
    structure assined for the inserted region in Loop.

    * Column 1 is the residue postion
    * Column 2 is the residue identity
    * Column 3 is the backbone behavior 

4.  Running the remodel application

        /mini/bin/remodel.macosgccrelease -s 3k2m_bc_relax.pdb -remodel:blueprint test.blueprint -extrachi_cutoff 1 -ex1 -ex2 -use_input_sc -num_trajectory 3 -save_top 1 -use_clusters false -database ~/minirosetta_database -find_neighbors


Method 2: Loopmodel and Fixbb
-----------------------------
We can use `loopmodel` with fragment files and `fixbb` with a resfile to design 
the monobody part of interface.

1.  Preparing the starting PDB

    Same as Step 1 in METHOD 1. We do not need to get rid of the Chai ID column 
    in this case. The relaxed PDB `3k2m_bc_relax.pdb` will be the input PDB.

2.  Creating fragment libraries

    Take the fasta file of 3k2m_bc_relax.pdb. Create  fragment libraries of 
    sizes 9 and 3 locally or through Robetta Server.

3.  Other input files

    Loop file: `3k2m.loop_file`

    The format is as follows:

        #LOOP  start end cutpoint skip-rate extend
        LOOP 85 89 0 0.0 0
        LOOP 179 185 0 0.0 0

    where

        column1  "LOOP":     The loop file identify tag
        column2  "integer":  Loop start residue number
        column3  "integer":  Loop end residue number
        column4  "integer":  Cut point residue number, >=startRes, <=endRes. default - let LoopRebuild choose cutpoint
        column5  "float":    Skip rate. default - never skip
        column6  "boolean":  Extend loop. Default false.

4.  Running the loopmodel application:

        ~/mini/bin/loopmodel.macosgccrelease  @flags

    where the flags file consist of following options (edit path to database!):

        -database ~/minirosetta_database
        -in:file:fullatom
        -loops:input_pdb 3k2m_bc_relax.pdb
        -loops:loop_file 3k2m.loop_file
        -loops:frag_sizes 9 3 1
        -loops:frag_files aat000_09_05.200_v1_3 aat000_03_05.200_v1_3 none
        -loops:remodel quick_ccd
        -loops:ccd_closure
        -loops:random_loop
        -out:prefix 3k2m_
        -mut core.io.database
        -nstruct 5

    The actual experiment should have 1000-10000 nstruct.

4.  Running the fixbb application with a resfile:

    Create a list of the output PDBs from loop modeling.

        ~/mini/bin/fixbb.macosgccrelease -database ~/minirosetta_database  -l list -resfile resfile -extrachi_cutoff 1 -ex1 -ex2 -nstruct 5

    Resfile format:

        NATAA
        start
        179 B ALLAAxc
        180 B ALLAAxc
        181 B ALLAAxc
        182 B ALLAAxc
        183 B ALLAAxc
        184 B ALLAAxc
        185 B ALLAAxc

5.  Select PDBs based on total_score

6.  Optional run for binding energy.

Interface Design Demo
=====================

Author: Kevin Houlihan (khouli at unc dot edu)

This demo contains scripts to redesign residue sequences at the interface of
existing heterodimers. Protocols for running with PackRotamers, MinMover,
and MinMover with a scorefunction produced by optE are included.

Setup
-----

This demo contains scripts with absolute paths that need to be updated. A
script is included to do this for you. (If you don't want to use it a list of
files with absolute paths that need to be set is at the bottom of this
section.)

Run the configure script. Pass your rosetta directory as an argument if it
isn't in your home directory. After setting your rosetta directory, the script
tries to add a weight file.

Example:

    [khouli@killdevil-login2 demo]$ ./configure 
    Use /nas02/home/k/h/khouli/rosetta as Rosetta directory for scripts? [Y/n]y
    /nas02/home/k/h/khouli/rosetta
    Add optE_inf_premin.wts to /nas02/home/k/h/khouli/rosetta/rosetta_database/scoring/weights/? [Y/n]y

Alternatively:

    [khouli@killdevil-login2 demo]$ ./configure other_dir/rosetta
    other_dir/rosetta
    Add optE_inf_premin.wts to other_dir/rosetta/rosetta_database/scoring/weights/? [Y/n]

This fixes absolute filepaths within the interface_design_demo/ directory.
Adding the optE_inf_premin.wts to your weights folder is necessary for the
protocol in scripts/design/minpac_optE_premin/ to run.

The configure script makes a few assumptions. You may want to change these.
What those assumptions are and where they can be found:

- The relative paths downwards from interface_design_demo/ have not been 
  changed e.g. scripts within the demo can access other scripts within the demo 
  at their expected paths

- Your rosetta path contains:

  `rosetta_database/`: referenced by design scripts, e.g. pacrot_s12p_nomin.sh 
  (run by run.sh)

  `rosetta_source/bin/sequence_recovery.linuxgccrelease`: referenced by 
  scripts/analyze/seqRec.sh (run by infpro.sh)

  `rosetta_source/bin/rosetta_scripts.mpi.linuxgccrelease`: referenced by 
  design scripts referenced by scripts/analyze/seqRec.sh (run by infpro.sh)

- Your mpi system is the same as mine and you want to use the same options 
  referenced in design scripts

This can also be done manually. Files with absolute paths:

- design scripts in /scripts/design/ e.g. 
  `scripts/design/pacrot_s12p_nomin/pacrot_s12p_nomin.sh`

- `scripts/analyze/seqRec.sh`

- `inputs/selected_chains/selected_chains.list`

- `inputs/min_nats/nats.list`

Running Design Protocols
------------------------

From within one of the protocol directories in scripts/design/, run the run.sh
script.

What this script is doing:

This script runs the actual design script (the one with the same name as the
directory) and passes it 3 arguments. The first argument names the job after
the directory (and hence creates an identically named sub-directory with output
pdbs), the second passes selected_chains.list as the list of proteins to use,
and the third passes the name of the .xml file containing the Rosetta Script to
be used.

Example:

    [khouli@killdevil-login2 pacrot_s12p_nomin]$ pwd
    /nas02/home/k/h/khouli/interface_design_demo/scripts/design/pacrot_s12p_nomin
    [khouli@killdevil-login2 pacrot_s12p_nomin]$ ./run.sh 
    Group   (-G) : bkuhlman_pi
    Project (-P) : bkuhlman_pi
    Memory Limit (-M) : 4 GB
    Job <287057> is submitted to queue <day>.

Change Rosetta run parameters within the <name-of-protocol>.sh script file.
Change the protocol described by the Rosetta script in the xml file.
Interface residues are determined by the RestrictToInterfaceVector task
operator in the Rosetta script.

Analysis
--------

Once a design run completes, analyze it by running post_run.sh from the same
directory the run.sh script was executed it. Inside the output directory, this
will generate the files top_interface_score, sequence_recovery.txt, and
submatrix.txt as well as some others of much less interest.

What this script is doing:

This script is a convenient way to execute the script
interface_design_demo/scripts/analyze/infpro.sh with the proper working
directory. `../../../analyze/infpro.sh` from within the output directory is
awkward.

In turn, infpro.sh is runs a script to score interfaces of output pdbs,
identify the best trajectories of each protein based on those scores, and
then runs a script to analyze sequence recovery and to aggregate substitutions
in residues between the inputs and designs.

This script assumes your rosetta path contains:

    rosetta_source/bin/sequence_recovery.linuxgccrelease

This can be changed in scripts/analyze/seqRec.sh.

To run infpro.sh, set the directory that conntains the output pdbs as your
working directory and then run interface_design_demo/scripts/analyze/infpro.sh.

Example:

    [khouli@killdevil-login2 pacrot_s12p_nomin]$ pwd
    /nas02/home/k/h/khouli/interface_design_demo/scripts/design/pacrot_s12p_nomin/pacrot_s12p_nomin
    [khouli@killdevil-login2 pacrot_s12p_nomin]$ ls sequencerecovery.txt
    ls: sequencerecovery.txt: No such file or directory
    [khouli@killdevil-login2 pacrot_s12p_nomin]$ ~/interface_design_demo/scripts/analyze/infpro.sh 
    [khouli@killdevil-login2 pacrot_s12p_nomin]$ ls sequencerecovery.txt 
    sequencerecovery.txt

The outputs worth looking at are sequence_recovery.txt, submatrix.txt,
and top_interface_scores.


Example ouputs
--------------

The directory example_outputs/ contains outputs generated by running the
configure script and then running each protocol as described above. Each
protocol produces one output folder as well as a benchmark.txt file. The
benchmark.txt files have been prefixed with their protocol but otherwise all
outputs are exactly as produced within the scripts/design/name_of_protocol/
directories.

Where to change parameters
--------------------------

* Adding/removing input heterodimer structures:

  In the run.sh file for the design, change the variable $input_pdb_list to the 
  filename of a text file containing a list of pdbs to use as inputs. 
  Alternatively, add pdbs to inputs/selected_chains/ and re-run configure. The 
  new pdbs will be added to inputs/selected_chains/selected_chains.list Either 
  way, if you want sequence recovery to function you'll need to add a 
  corresponding structure to the set used for sequence recovery.

* Adding/removing native structures for sequence recovery:

  In the post_run.sh file for the design, change the assignment of 
  $native_pdb_list to the filename of a text file containing a list of pdbs to 
  use as natives. Alternatively, add pdbs to inputs/min_nats/ and re-run 
  configure. The new pdbs will be added to inputs/min_nats/nats.list. If you 
  want sequence recovery to function you'll need to have added a corresponding 
  structure to the set used for input. (They have to align by simple ls, i.e. 
  cap-sensitive non-numerical alphabetical order). 

* Number of trajectories used for sequence recovery:

  In the run.sh file, change the n_traj assignment at the top.

* Options passed to Rosetta and MPI:

  In a protocol directory, edit name_of_protocol.sh

Troubleshooting
---------------

PMGR_COLLECTIVE ERROR following post_run.sh: Fix the path to the 
sequence_recovery binary in `interface_design_demo/scripts/analyze/seqRec.sh`

Epilogue
--------

Happy Rosetta-ing!
To run this script as written, the optE_inf_premin.wts weights file must be
added to "rosetta/rosetta_database/scoring/weights/". The configure script in
the top-most directory of this demo attempts to add it automatically.
These are the scripts/design/*/ directories cp'd after a successsful design
run. The <name-of-protocol>.sh scripts were lightly modified before the runs
to meet the whims of the one particular LSF system and so they may not run for
you. Since the directories have been moved, the selected_chains.list
sym link files are broken since they use relative paths. They are included
as part of a full demonstration of the final state of the directories.
To run this script as written, the optE_inf_premin.wts weights file must be
added to "rosetta/rosetta_database/scoring/weights/". The configure script in
the top-most directory of this demo attempts to add it automatically.
# optE

This demo contains the input files and command lines
necessary to refit the reference energies for the
Talaris2013 score function using optE's sequence-
profile-recovery procedure.

The demo should be run from within the
  run/
directory, using the launch_optE.scr script.  The
demo is set up to run on 27 processors, but can be
run on fewer processors.  If the demo is running
too slowly, comment out the -ex1 and -ex2 flags
from the optE_seqprof.flags file. There  are several
things about the launch_optE.scr script that will have
to be edited.

1) the path to the executable, in launch_optE.scr
2) the mpi launch command, in launch_optE.scr, and
3) the path to the rosetta database, in score.flags

After optE has finished running, use the python script,
find_best_weight_set.py, in the
  scripts/
directory, to determine which of the iterations was
optimal.
These are the 54 PDBs in HiQ54 dataset, curated by Jane Richardson. You can see a full description of the curation process here:

http://kinemage.biochem.duke.edu/databases/HiQ54.phpThis is the README.dox file for the run_rosettaholes_protocol_on_ubiquitin.
# Supercharge

If you want to run supercharge now, the application is called 'supercharge' in `src/apps/public/supercharge.cc`.

Here are four examples:

```
./supercharge.default.macosgccrelease @rosetta_inputs/options1 -database <path>           // Rosetta-mode, positive-charge, fixed surface cutoff and input ref energies
./supercharge.default.macosgccrelease @rosetta_inputs/options2 -database <path>		  // Rosetta-mode, negative-charge, fixed surface cutoff and target net charge
./supercharge.default.macosgccrelease @rosetta_inputs/options3 -database <path>		  // AvNAPSA-mode, negative-charge, target net charge
./supercharge.default.macosgccrelease @rosetta_inputs/options4 -database <path>		  // AvNAPSA-mode, positive-charge, fixed surface cutoff
```


Rosetta-mode and AvNAPSA-mode are explained below...

# Why supercharge protein surfaces?

Reengineering protein surfaces to have high net charge, called supercharging, can improve reversibility of unfolding by preventing aggregation of partially unfolded states. Aggregation is a common obstacle for use of proteins in biotechnology and medicine.  Additionally, highly cationic proteins and peptides are capable of nonviral cell entry, and highly anionic proteins are filtered by kidneys more slowly than neutral or cationic proteins.  

Optimal positions for incorporation of charged side chains should be determined, as numerous mutations and accumulation of like-charges can also destabilize the native state.  A previously demonstrated approach deterministically mutates flexible polar residues (amino acids DERKNQ) with the fewest average neighboring atoms per side chain atom (AvNAPSA: Lawrence MS, Phillips KJ, Liu DR, 2007, Supercharging proteins can impart unusual resilience, JACS).  Our approach uses Rosetta-based energy calculations to choose the surface mutations.  Both automated approaches for supercharging are implemented in this online server.

# Two Approaches
There are two automated approaches, **Rosetta supercharge (Rsc)** and **AvNAPSA supercharge (Asc)**

**AvNAPSA supercharge philosophy (Asc):** mutate the most exposed polar residues to minimize structural change or destabilization.  Only DE-RK-NQ residues can be mutated.

**Rosetta supercharge philosophy (Rsc):** mutate residue positions that preserve and/or add favorable surface interactions.  Hydrophobic and small polar surface residues can also be mutated.

**AvNAPSA drawbacks:** mutating surface polar residues can eliminate hydrogen bonds.  Helix capping, edge-strand interaction, and loop stabilization all result from surface hydrogen bonds.  Furthermore, this automated protocol mutates N to D and Q to E, but N and Q sometimes act simultaneously as a donor and acceptor for hydrogen bonds.

**Rosetta drawbacks:** mutating less-exposed positions can lead to better computed energies, but mistakes at these positions can be destabilizing.  AvNAPSA favors charge swaps, so Rosetta requires more mutations to accomplish the same net charge.

The AvNAPSA approach varies net charge by adjusting the surface cutoff.  The Rosetta approach varies net charge by adjusting reference energies of the positive or negatively charged residues.

The supercharge server can run in four different modes:
-AvNAPSA with a target net charge
-AvNAPSA with a surface cutoff
-Rosetta with a surface cutoff and target net charge
-Rosetta with a surface cutoff and input reference energies for charged residue types


**What does AvNAPSA stand for:** average number of neighboring atoms per sidechain atom.  This is a value that measures the extent of burial/accessibility.  It's similar to the residue neighbors by distance that Rosetta typically uses to define the surface, but it's on the atom-level rather than residue-level.  AvNAPSA-mode calculates an AvNAPSA value for every residue.  'surface_atom_cutoff' indicates the cutoff AvNAPSA value that defines surface residues.  AvNAPSA values of 50-150 are typical for surface residues.  AvNAPSA values >150 are typical for core residues.  A surface_atom_cutoff of 100 will lead to moderate supercharging.  A surface_atom_cutoff of 150 will lead to heavier supercharging.


## Workflow of Each Mode

### AvNAPSA-mode, target charge
1. Define surface.  sort NQ and RK/DE residues by AvNAPSA value (low to high)
2. Next residue in sorted list: Positive: mutate DENQ-->K, Negative: mutate RKQ-->E and N-->D
3. If net charge = target net charge, output pdb

### AvNAPSA-mode, surface cutoff
1. Define surface by AvNAPSA value (<100 default)
2. For each NQ and DE/RK residue in the surface: Positive: mutate DENQ-->K, Negative: mutate RKQ-->E and N-->D
3. Output pdb

### Rosetta-mode, surface cutoff and target charge
1. Define surface.  Neighbor by distance calculator (CB dist.), <16 neighbors default or Define surface by AvNAPSA value (<100 default)
2. Set design task
    - read user resfile, if provided
    - dont_mutate gly, pro, cys
    - dont_mutate h-bonded sidechains
    - dont_mutate correct charge residues
3. Set reference energies for RK/DE, starting at user input values
4. pack rotamers mover
5. check net charge, increment/decrement reference energies (back to step 3.)
6. Once a pack rotamers run results in the correct net charge, output pdb

### Rosetta-mode, surface cutoff and input reference energies for charged residue types
1. Define surface.  Neighbor by distance calculator (CB dist.), <16 neighbors default or Define surface by AvNAPSA value (<100 default)
2. Set design task
    - read user resfile, if provided
    - dont_mutate gly, pro, cys
    - dont_mutate h-bonded sidechains
    - dont_mutate correct charge residues
3. Set reference energies for RK/DE, using the user input values
4. pack rotamers mover
5. Output pdb


# Options

## AvNAPSA Mode

```
AvNAPSA_positive  BOOL def(false);				//run positive-charge AvNAPSA
AvNAPSA_negative  BOOL def(false); 				//run negative-charge AvNAPSA
target_net_charge  SIGNED_INT def(0);  				//residue positions will be mutated one at a time from most exposed to least exposed until target net charge is achieved
surface_atom_cutoff  UNSIGNED_INT def(100); 			// if you have no target net charge in mind, AvNAPSA will mutate all surface DE-RK-NQ residues on the surface, with this surface cutoff
```

## Rosetta Mode

```
surface_residue_cutoff  UNSIGNED_INT def(16);  //residues with <16 neighboring residues within 10 Å are considered part of the surface

include_arg  BOOL def(false);  //use arginine in Rosetta supercharge
include_lys  BOOL def(false);  //use lysine in Rosetta supercharge
include_asp  BOOL def(false);  //use aspartate in Rosetta supercharge
include_glu  BOOL def(false);  //use glutamate in Rosetta supercharge

//the reference energies of the charged residue types will govern the net charge of Rosetta designs.  Rosetta can choose between the allowed charged residue types and the native residue.  More negative reference energies will result in more charge mutations.
refweight_arg  FLOAT def(-0.98);
refweight_lys  FLOAT def(-0.65);
refweight_asp  FLOAT def(-0.67);
refweight_glu  FLOAT def(-0.81);

dont_mutate_glyprocys  BOOL def(true);		 //glycine, proline, and cysteine often serve special structural roles in proteins
dont_mutate_correct_charge  BOOL def(true);      //i.e., don’t mutate arginine to lysine
dont_mutate_hbonded_sidechains  BOOL def(true);  //don’t mutate residues with sidechains forming a hydrogen bond
pre_packminpack  BOOL def(false);                //Packrotamers is always done as the first step.  This option will go one step further and run packrotamers, sidechain+backbone minimization, packrotamers on the input structure before performing the supercharge design step.

nstruct  UNSIGNED_INT def(1);  			 //Monte Carlo sequence design of a protein surface is often convergent but it is still stochastic, multiple design runs can be performed if desired.
target_net_charge  UNSIGNED_INT def(0);  	 //a target net charge can be achieved if desired, this is done in an automated way by incrementing/decrementing charged residue reference energies until the desired net charge results from the Monte Carlo design step.
```

## AvNAPSA and Rosetta Mode

```
surface_atom_cutoff  UNSIGNED_INT def(100); // this is how AvNAPSA defines surface, can be used in either approach
compare_energies  BOOL def(false);  	 	      		//prints a full residue-by-residue energy analysis in the log file
only_compare_mutated_residues  BOOL def(false);  		//only includes mutated residues in the energy analysis
resfile  FILE;  	       	    				//this is how you can specify which residues to not mutate.  Default setting must be ALLAA, and residue-by-residue settings should be NATAA, as shown below:

ALLAA
start
  20  A  NATAA
  24  A  NATAA
  26  A  NATAA


Note: an input resfile is optional.  However, every supercharge run generates an output resfile that governs the design run.  The default of this output resfile is NATAA, which prevents core residues from mutating (see below).  The input resfile is read first, the output resfile (see below) is read second, and this is why ALLAA must be the default for the input resfile.  If the default were NATRO, for example, no design would occur!
```

# Output

As output, a log file, the residue file that governed the design run, and the output PDB are provided.  First, the log file contains the exact Rosetta command line, the residue positions identified as located on the surface, a list of charged residues in the final sequence, the net charge, a list of mutations, text for a PyMOL selection to easily view the mutations in PyMOL, and optionally, a full energetic comparison of repacked native versus supercharged structures.  Secondly, the Rosetta residue file indicates which residue positions could possibly mutate, and to what residue types.  The third output file is the atomic coordinate file of the supercharged protein, in PDB format, and the naming of the output PDB is intended to facilitate self-documentation of the inputs for a given design run.  For Rosetta designs, the name includes the final reference energies that were used and the final net charge, and for AvNAPSA designs, the name includes the net charge and the largest AvNAPSA value of the mutated residues. 


This is what an output resfile looks like for AvNAPSA-positive supercharging, which always chooses lysine:

```
NATAA
start
   6 A  PIKAA  K
   19 A  PIKAA  K
   21 A  PIKAA  K
   32 A  PIKAA  K
   34 A  PIKAA  K
   39 A  PIKAA  K
```

This is what an output resfile looks like for Rosetta positive supercharging, which allows choice between native and RK, and preserves h-bonds:

```
NATAA
start
   6 A  PIKAA ERK
   9 A  PIKAA TRK
   11 A  PIKAA VRK
   21 A  PIKAA DRK
   25 A  PIKAA HRK
   26 A  NATAA  #same charge
   30 A  NATAA  #same charge
   32 A  NATRO  #has sc hbond energy=-1.15844
   38 A  PIKAA TRK
   39 A  NATRO  #has sc hbond energy=-1.33149
   43 A  PIKAA TRK
   50 A  NATRO  #has sc hbond energy=-0.536622
   52 A  NATAA  #same charge
   76 A  PIKAA DRK
   77 A  PIKAA HRK
```

# Symmetric Docking of Insulin Trimer of Dimers

Input files:

1. one pdb file
2. one symmetry definition file

(1) The pdb file represents the monomer entity; regardless if this file has multiple physical chains or its own internal symmetry this, is the monomer. However, for safe running of this script the pdbfile must be rewritten to have a single chain label.

(2) The symmetry definition file informs Rosetta about the symmetry type, number of monomers to dock, how to calculate the score and some more things.

# Example 1: Docking of an Insulin Dimer (chains A, B, C, D in 1ZEH) into a Trimer of Dimers (C3 symmetry)

## Clean Up the PDB File

IF the pdb file is multi-chain we must normalize it to a pseudo single chain with the following convenience script "addChain.pl". You can find that script in the scripts directory of this tutorial. This step is only required when the pdb is multi-chain!

```
	> ./scripts/addChain.pl 1ZEH.pdb A | grep '^ATOM'  >  rosetta_inputs/1ZEH_monomer_c3.pdb
```

where 1ZEH.pdb is the pdb file input and A is the the dummy new chain ID.
This simple script only normalizes the ATOM records; it removes the TER chain separators and deletes unnecessary SEQRES, REMARK, and other lines.

In the example here, 1ZEH.pdb, contains 4 chains, and the 1ZEH_monomer.pdb has a single chain joining all the ATOM records of these into one chain called A. Now we have the pdb input consisting of only our chosen monomer entity.

## Generating the Symmetry File

To generate a de novo symmetry file use script: `make_symmdef_file_denovo.py` that comes with the Rosetta code.
The minimal inputs of this script are (1) the symmetry: cn (circular symmetry) or dn (dihedral symmetry) and (2) the number of "monomers" (i.e. subunits).

```
	> src/apps/public/symmetry/make_symmdef_file_denovo.py -symm_type cn -nsub 3  > rosetta_inputs/1ZEH.c3.symm
```

where "cn" is cyclic symmetry and "3" is the number of "monomers" (that gives a 3 fold single axis rotational symmetry). This produces the symmetry definition file for Rosetta.

**Note:** Pay no attention to the fact that 1ZEH has it's own internal symmetries. It is just considered a monomer subunit.

The output file 1ZEH.c3.symm contains the following:

```
	symmetry_name c3
	subunits 3
	recenter
	number_of_interfaces  1
	E = 3*VRT0001 + 3*(VRT0001:VRT0002)
	anchor_residue COM
	virtual_transforms_start
	start -1,0,0 0,1,0 0,0,0
	rot Rz 3
	virtual_transforms_stop
	connect_virtual JUMP1 VRT0001 VRT0002
	connect_virtual JUMP2 VRT0002 VRT0003
	set_dof BASEJUMP x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)
```

## Run Rosetta

```
	> source/bin/SymDock.linuxgccrelease @flags.c3
```

where the provided flags file (flags.c3) looks like this:

```
	-database ~/Rosetta/main/database   # change path to your directory!!!
 	-in:file:s 1ZEH_monomer.pdb
 	-symmetry:symmetry_definition 1ZEH.c3.symm
 	-packing:ex1
 	-packing:ex2aro
 	-out:nstruct 1   # typically this should be 10,000 or more.  reduced to 1 for demo.
 	-out:file:fullatom
 	-symmetry:initialize_rigid_body_dofs
 	#-symmetry:symmetric_rmsd     # this line is only used when a reference pdb at the full C3 symmetry is available for testing purposes.
 	-ignore_unrecognized_res
```

## Results

The Rosetta run generated 1ZEH_monomer_c3_0001.pdb and score.fasc. The pdb file has the trimer decoy of the monomer docked in C3 symmetry.  The 3 monomers are entered as chains A,B and C in the trimer. The score.fasc is a standard score file as one would get for any normal docking run. The score is the score of the complete trimer. (See the PLOS one paper "Modeling Symmetric Macromolecular Structures in Rosetta3" or the symmetrical docking documentation to understand how this is computed. Conceptually you may think of it as the internal score of the monomer plus the score of it's interface; all multiplied by 3.)

Now for a nominally trickier example.

# Example 2: Docking of Six Monomers (chains A, B in 1ZEH) into a Trime of Dimers (D3 symmetry)

# Clean Up the PDB File

If you had taken a peek at the file 1ZEH.pdb beforehand you might have noticed that it already had an internal dimer symmetry.  The four chains were arranged so that chain A and B were in C2 symmetry to C and D.  So the above example just computed a trimer of dimers.  The program did not know anything about the c2 symmetry because it just considered the file 1ZEH_monomer.pdb to be the subunit for the trimer.

Now we are going to run this again but first removing chains C and D from the 1ZEH.pdb file.

```
	> perl -wane 'print if m/^ATOM/ and ($F[4] eq "A" or $F[4] eq "B")'  1ZEH.pdb  > 1ZEH_monomer_d3.pdb (can re removed after the next step)
```

And re-chain it to a single pseudo chain.

```
	> ./scripts/addChain.pl /tmp/1ZEH_monomer_d3.pdb A > rosetta_inputs/1ZEH_monomer_d3.pdb
```

## Generating the Symmetry File

Then we will use d3 symmetry for docking.  As you know, d3 symmetry is c3 symmetry with a mirror plane. In other words this has a total of 6 subunits.

```
	> src/apps/public/symmetry/make_symmdef_file_denovo.py -symm_type dn -nsub 6  > rosetta_inputs/1ZEH.d3.symm
```

**Note:** the number of subunits in d3 symmetry is 6.

Here is the symmetry definition file 1ZEH.d3.symm

```
symmetry_name d3
subunits 6
recenter
number_of_interfaces  4
E = 6*VRT0001 + 6*(VRT0001:VRT0002) + 3*(VRT0001:VRT0004) + 3*(VRT0001:VRT0005) + 3*(VRT0001:VRT0006)
anchor_residue COM
virtual_transforms_start
start -1,0,0 0,1,0 0,0,0
rot Rz_angle 120.0
rot Rz_angle 120.0
rot Rx_angle 180.0
rot Rz_angle 120.0
rot Rz_angle 120.0
virtual_transforms_stop
connect_virtual JUMP1 VRT0001 VRT0002
connect_virtual JUMP2 VRT0002 VRT0003
connect_virtual JUMP3 VRT0003 VRT0004
connect_virtual JUMP4 VRT0004 VRT0005
connect_virtual JUMP5 VRT0005 VRT0006
set_dof BASEJUMP x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)
set_dof JUMP3 z(50) angle_z(0:60.0)
```

## Run Rosetta

```
	> source/bin/SymDock.linuxgccrelease @flags.c3
```

The result is again a pdbfile 1ZEH_haptomer_00001.pdb which has 6 monomers in d3 symmetry.   Note, again, the monomer in this second example had half as many atoms as the first example but since we used 6 of these instead of 3 the final structure has the same number of total atoms.  The score is appended into score.fasc

# Remarks

Not shown in this example is the optional but recommended practice of "pre-packing".  The docking protocol packs the rotamers of all residues in the interfaces.  However it does not pack any residues that are not in the interface.  Therefore a recommended practice is the pack the monomer ahead of time.  This is simply running the packer on the monomer file with any desired packing options.  Then follow the above protocol as written to the symmetric docking.  A possible point of confusion here is that there is a flag for docking named "prepack" however this does not pack the interior as desired here.  Instead use the Docking_prepack application (see below).


There are other flags associated with both "docking" and "symmetry" available beyond the ones show in this tutorial.  These flags are documented in their respective formal documentation section.

# Links

- https://www.rosettacommons.org/docs/latest/symmetry.html
- https://www.rosettacommons.org/docs/latest/sym-dock.html
- https://www.rosettacommons.org/docs/latest/docking-protocol.html
- https://www.rosettacommons.org/docs/latest/docking-prepack-protocol.html

# Appendix

## Summary of the Commands

```
./scripts/addChain.pl 1ZEH.pdb A | grep '^ATOM'  >  1ZEH_monomer.pdb
src/apps/public/symmetry/make_symmdef_file_denovo.py -symm_type cn -nsub 3  > 1ZEH.c3.symm
bin/SymDock.macosgccrelease @flags.c3  >& 1ZEH.c3.log
perl -wane 'print if m/^ATOM/ and ($F[4] eq "A" or $F[4] eq "B")'  1ZEH.pdb  > /tmp/1ZEH_haptomer.pdb
./scripts/addChain.pl /tmp/1ZEH_haptomer.pdb A > 1ZEH_haptomer.pdb
src/apps/public/symmetry/make_symmdef_file_denovo.py -symm_type dn -nsub 6  > 1ZEH.d3.symm
bin/SymDock.macosgccrelease @flags.d3  >& 1ZEH.d3.log
```

## Authors
- Sebastian Rämisch
- Charlie E. M. Strauss
- Bruno Correia

**Rosetta Revision: 43611**

**Date: 08/05/2011**
# RNA De Novo

```
rna_denovo.linuxgccrelease @flags -database <PATH TO ROSETTA DATABASE>
```
Ligand Dock Demo
================

To run:

    [path]/rosetta/rosetta_source/bin/ligand_dock.linuxgccrelease -database ~/rosetta/rosetta_database/ @flags

Real run has nstruct 500 e.g. and run 10 separate runs for 5000 total runs.

Gives output:

    silent.out

To extract the scores, run this:

    [path]/rosetta/rosetta_source/src/apps/public/ligand_docking/get_scores.py
    <silent.out > scores.tab

Then you can extract a structure:

    ~/rosetta/rosetta_source/bin/extract_atomtree_diffs.macosgccrelease
    @flagsextract -database ~/rosetta/rosetta_database

Here there's only one structure and one tag for that structure. In a
real case you'd have many more:

    [path]/rosetta/rosetta_source/src/apps/public/ligand_docking/best_ifaceE.py
    -n 1 silent.out

For more information, see the ligand dock entry in the manual:  
http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/app_ligand_docking.html

Author: Barak Raveh

Overview:
=========
This demo demonstrates ab-initio folding and docking of the peptide ssrA (chain D) to the receptor protein sspB (chain B) at a specified binding sites. The user supplies Rosetta with a model of the sspB receptor protein (bound model, for the purpose of this demo), the sequence of the ssrA peptide, and a list of atoms on the surface of sspB that are known to bind the peptide.

Note on identifying the binding site:
================================================
At this stage, Rosetta FlexPepDock allows ab-initio folding and docking of a peptide only at a specifie\
d binding site (completely blind global docking is under development in the Furman lab). That said, for\
tunately peptides have a nice property - many of them tend to bind to large grooves and pockets over th\
e receptor protein surface (London et al., Structure 2010 - The Structural Basis of Peptide-Protein Bin\
ding Strategies). While this does not hold true in 100% of the cases, this is certainly the case for th\
e provided example. Therefore, this demo provides a good example for ab-initio folding and docking of a\
 peptide at a binding site that is specified as a list of atoms or residues over the surface of the rec\
eptor protein. In this case, the binding site residues (specified in <input_files/site_constraints.cst)\
 is based on the Pocket-Finder server (http://www.modelling.leeds.ac.uk/pocketfinder/), which uses the \
Ligsite pocket detection program (Hendlich et al., 1997). Note that other pocket detection prorgrams mi\
ght also work for the same purpose. This would not work in all cases, but here the binding pocket is well defined.

General Sceheme:
================
In order to identify the peptide binding site at the largest pocket, and fold-and-dock the peptide at the binding site, we worked as follows:
(1) We used the Pocket-Finder server (http://www.modelling.leeds.ac.uk/cgi-bin/pocketfinder/pfmage.cgi) on the receptor chain (chain B), and identified the largest pocket. In principle, any other pocket detector can be used for this task. The resulting residues list copied from the server is found in <input_files/pocket_finder_site.txt>.
(2) Based on the site identified in step 1, we created a constraint file that biases FlexPepDock ab-initio to sample peptide conformations at the vicinity of the identified pocket, using soft constraints that can be violated quite generously, to prevent over-bias in case of noise in the input. The constraints for side-chain atoms were repklaced with constraints on the centroid atom, since these constraints are used in the centroid low-resolution step of FlexPepDock. The resulting constrains file is found in <input_files/site_constraints.cst>.
(3) We ran FlexPepDock ab-initio + FlexPepDock Refinement to create the high-resolution model of the interaction. In the starting model, we positioned the peptide in an arbitrary initial orientation relative to the receptor protein, and set its backbone conformation to an ideal extended backbone conformation. The constraint file was used in the initial low-resolution stage of FlexPepDock ab-initio to force the peptide to contact the receptor at the vicinity of the specified binding pocket. 

Details:
========need to fix broken link:
global_dock_ssrA_peptide_against_sspB/input_files/frags/vall.dat.2006-05-05
The file is removed for now.Membrane Homology Modeling
==========================

This demo shows how to generate a homology model of a membrane protein.  It 
only works for helical membrane proteins.  In the example, a dopamine receptor 
structure is used as the template for a homologous dopamine receptor.  This is 
a relatively easy homology modeling problem.

Preparation of input files
--------------------------

There is a non-trivial amount of preparatory work that is required before 
running a homology modeling task in Rosetta.

1.  Find a structural template for homology modeling.

    The first step is to find a structural alignment to a template protein.  In 
    this demo, the fasta sequence file for the target (`DXDR_.fasta`) was 
    pasted into the input panel of the HHpred web server.  The server gives you 
    a list of structural homologs ranked by the quality of alignment.  For 
    each, you get a secondary structure prediction and a template structure. 
    The next task is to manipulate this information to create a combined fasta 
    file that shows the alignment between the target and the template file.  In 
    this demo this combined fasta file is named `DXDR_D3DR.HHp`.

2.  Process the template structure to accept the target sequenced.

    Next you need to acquire the pdb file for the template, and you need to 
    renumber the residues in the file to match the target protein numbering. 
    Scripts for performing this are provided.  First, the perl script 
    fastAln2zones.pl is invoked on the combined fasta file:

        fastAln2zones.pl DXDR_D3DR.HHp <zones_file>

    This creates a "zones" file that shows the correspondence between residues 
    in the target and residues in the template.

    Using several of these pieces of information, the script 
    `createTemplate.pl` creates the renumbered template structure:

        createTemplate.pl -zonesfile <zones_file> -fastafile DXDR_.fasta -parentpdb 3PBL_A_renum.pdb -outpdb DXDR_D3DR.pdb

    3PBL_A_renum.pdb was previously created by taking 3PBL.pdb and cutting out 
    the T4 lysozyme insertion from a loop in the membrane domain, and 
    renumbering the remaining residues with the awk script `renum.awk`. 
    Similar manipulations may be required for each specific application.

    The end product of step 2 is the file DXDR_D3DR.pdb

3.  Create a fragment library for the target protein.

    You also need to generate a fragments file as for any other loop/homology 
    modeling application.  See the appropriate demo for fragment generation.

4.  Identify the regions of the protein that must be rebuilt by Rosetta.

    Now we need to tell Rosetta the regions that must be rebuild by Rosetta. 
    These are the regions that are poorly defined.  This is defined in 
    DXDR.loopfile.  The format is:  first line, a comment, other lines are 
    `LOOP start-res end-res cutpoint skip-rate extend` for each flexible region 
    to rebuild.  The cutpoint column is used to tell Rosetta where to cut the 
    loop.  A zero tells Rosetta to use its default.  Skip rate allows you to 
    spend more time on one loop versus the others.  Here we don't ask for this. 
    Extend "X" indicates not to use any information from the starting template 
    as part of the flexible regions.

5.  Create membrane specific input files.

    There are two membrane protein-specific input files that are required.  The 
    first is a weight file with specialized scoring function weights.  It also 
    triggers the initialization of a number of membrane-specific scoring 
    function capabilities.  The second membrane specific file tells Rosetta 
    where the membrane is.  This file contains more information than is 
    necessary.  It is called DXDR.span.  The format is:

    * Line 1 is a comment.
    * Line 2 is the number of transmembrane helices and total number of 
      residues.
    * Line 3 tells about the topology, only important for ab initio folding, 
      but part of the format.
    * Line 4+ tell the bounds of the transmembrane helices in residue pairs. 
      Only the first pair is used for homology modeling.  The number of lines 
      must match the number of helices given above.

To sum up, we now have:

* `DXDR_D3DR.pdb` – the template pdb file generated above
* `DXDR.span` – the membrane helix information
* `Frags/` a directory with fragment data for the target protein
* `DXDR.loopfile` – the file that tells which regions are to be rebuilt

Running the demo
----------------

The following command line options are used as arguments to the minirosetta 
executable (minirosetta.(os)(options)(mode), for example:

    minirosetta.linuxgccrelease
    minirosetta.macosclangrelease

Provide the path to your own rosetta database:

    -database /Users/patrickbarth/RosettaCon2011/tutorial/trunk_r43621/minirosetta_database

Membrane proteins use the same protocol as regular proteins:

    -run:protocol looprelax

Information for fragments files:

    -loops:frag_sizes 9 3 1
    -loops:frag_files ./frags/aaDXDR_09_05.200_v1_3 ./frags/aaDXDR_03_05.200_v1_3 none

Some input information.  For benchmark purposes, replace the following PDB by 
the experimentally-determined native structure:

    -in:file:native ./rosetta_inputs/DXDR_D3DR.pdb
    -in:file:fullatom
    -s DXDR.pdb

The following PDB corresponds to the starting template:

    -loops:input_pdb /Users/patrickbarth/RosettaCon2011/tutorial/dopamine/input/DXDR_D3DR.pdb

Loop file generated from the zone file:

    -loops:loop_file /Users/patrickbarth/RosettaCon2011/tutorial/dopamine/input/DXDR.loopfile

Membrane specific input information:

    -score:weights ./rosetta_inputs/membrane_highres_t1.wts
    -in:file:spanfile ./rosetta_inputs/DXDR.span

Output format:

    -out:file:fullatom
    -out:file:silent_struct_type binary
    -out:file:silent DXDR_lprlxmb.out

Number of structures generated:

    -nstruct 1

Loop remodeling step: for additional details, see the regular looprelax 
protocol options:

    -loops:remodel  quick_ccd
    #-loops:refine  refine_ccd
    -loops:random_order
    -loops:idealize_before_loop_close

If you want FULL STRUCTURE relax then do this, otherwise don't:

    -loops::relax    fastrelax

Fail on bad H-bonds:

    -fail_on_bad_hbond false

The expected output is a silent file with the model coordinates.  In this demo, 
the name of the file is:

    DXDR_lprlxmb.out.

Membrane ab initio modeling
===========================

This document was created on August 5, 2011 by:

* Vladimir Yarov-Yarovoy (yarovoy@ucdavis.edu)
* Ingemar Andre (ingemar.andre@biochemistry.lu.se)
* Joe Jardine (jardinejg@gmail.com)
* Daniel Keedy (daniel.keedy@gmail.com)

This protocol was developed to predict helical membrane protein structures.

Algorithm
---------

This protocol will only generate low-resolution centroid models.  The protocol 
is using transmembrane region predictions from OCTOPUS server 
(http://octopus.cbr.su.se/) to set initial membrane normal and membrane center 
vectors and define membrane-specific environment (hydrophobic core, interface, 
polar, and water layers). For multispan helical membrane proteins the protocol 
starts from embedding only two randomly selected adjacent transmembrane helical 
regions and then continues folding by inserting one of adjacent helices until 
all transmembrane helices will be embedded into the membrane.

Input Files
-----------

1.  Generate structure fragments:

        make_fragments.pl -verbose -id BRD4 BRD4_.fasta -nojufo -nopsipred

    Documentation at 
    http://www.rosettacommons.org/manuals/archive/rosetta3.3_user_guide/file_fragments.html. 
    Note: use only SAM secondary structure prediction file (\*.rdb). jufo and 
    psipred predict transmembrane helical regions poorly.

2.  Genarate transmembrane regions (OCTOPUS) file:

    Input OCTOPUS topology file is generated at http://octopus.cbr.su.se/ using 
    protein sequence as input.

    Sample OCTOPUS topology file:

        ##############################################################################
        OCTOPUS result file
        Generated from http://octopus.cbr.su.se/ at 2008-09-18 21:06:32
        Total request time: 6.69 seconds.
        ##############################################################################

        Sequence name: BRD4
        Sequence length: 123 aa.
        Sequence:
        PIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFV
        WWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGA
        GIV

        OCTOPUS predicted topology:
        oooooMMMMMMMMMMMMMMMMMMMMMiiiiMMMMMMMMMMMMMMMMMMMMMooooooMMM
        MMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMM
        ooo

3.  Convert OCTOPUS file to .span file format:

    `BRD4.span`: transmembrane topology prediction file generated using 
    octopus2span.pl script as follows:

        octopus2span.pl <OCTOPUS topology file>

    Example:

        <path to mini>/mini/src/apps/public/membrane_abinitio/octopus2span.pl BRD4.octopus

    Sample .span file:

        TM region prediction for BRD4 predicted using OCTOPUS
        4 123
        antiparallel
        n2c
           6    26     6    26
          31    51    31    51
          58    78    58    78
          97   117    97   117

    Format for .span files:

    * 1st line is comment line.

    * 2nd line shows number of predicted transmembrane helices (4 in the 
      command lines example below) and total number of residues (123 in the 
      example below).

    * 3rd line shows predicted topology of transmembrane helices in the 
      membrane (currently only antiparallel topology is implemented).

    * 4th line and all lines below show start and end residue numbers of each 
      of the predicted transmembrane helices (current format repeats these 
      numbers twice).

4.  Generate .lips4 file.

    `BRD4.lips4`: lipophilicity prediction file created using run_lips.pl 
    script as follows (note that blastpgp and nr database are necessary to run 
    run_lips.pl script:

        run_lips.pl <fasta file> <span file> <path to blastpgp> <path to nr database> <path to alignblast.pl script>

    Example:
    
        <path to mini>/mini/src/apps/public/membrane_abinitio/run_lips.pl BRD4.fasta BRD4.span /work/bjornw/Apps/blast/bin/blastpgp /scratch/shared/genomes/nr ~bjornw/mini/src/apps/public/membrane_abinitio/alignblast.pl

    Sample lips4 file:

        Lipid exposed data: resnum mean-lipo lipophil entropy
              6  -1.000   3.004   1.211
              9  -1.000   2.268   2.137
             10  -1.000   4.862   1.095
             13  -1.000   1.304   1.552
             16  -1.000   3.328   2.025
        ...

5.  Run membrane *ab initio* application with the following flags (see an 
    example in scripts/membrane-abinitio.cmd):

        ./bin/membrane_abinitio2.linuxgccrelease
        -in:file:native BRD4.pdb                  Native structure (optional)
        (or -in:file:fasta BRD4_.fasta)           Protein sequence in fasta format (required if native structure is not provided)
        -in:file:spanfile BRD4.span               Octopus transmembrane prediction (see above)
        -in:file:lipofile BRD4.lips4              Lipophilicity prediction (see above)
        -in:file:frag3 aaBRD4_03_05.200_v1_3      3-residue fragments
        -in:file:frag9 aaBRD4_09_05.200_v1_3      9-residue fragments
        -in:path:database ~minirosetta_database   Path to rosetta database
        -abinitio:membrane                        Membrane ab initio application
        -score:find_neighbors_3dgrid              Use a 3D lookup table for residue neighbors calculations
        -membrane:no_interpolate_Mpair            Switch off the interpolation between the two layers for the Mpair term
        -membrane:Menv_penalties                  Switch on the following penalties:
                                                  1. no non-helical secondary structure fragments in predicted transmembrane helical regions in the hydrophobic layer of the membrane.
                                                  2. no N- and C- termini residues of predicted transmembrane helical regions in the hydrophobic layer of the membrane.
                                                  3. no transmembrane helices with orientation >45 degrees relative to the membrane normal vector.
        -nstruct 1                                Number of output structures

    You have a choice to use either Monte Carlo (default) or discrete search of membrane normal and membrane center.

    Optional settings used for Monte Carlo based membrane normal and center search protocol:

        -membrane:normal_cycles (default=100)     Total number of membrane normal cycles
        -membrane:normal_mag (default=5)          Magnitude of membrane normal angle search step size (degrees)
        -membrane:center_mag (default=1)          Magnitude of membrane center search step size (Angstroms)

    Tip: to speedup Monte Carlo based membrane normal and center search use the following settings:

        -membrane:normal_cycles 40
        -membrane:normal_mag 15
        -membrane:center_mag 2

    Optional settings for alternative discrete search of membrane normal and membrane center:

        -membrane:center_search (default= false) - perform membrane center search within "center_max_delta" deviation (see below).
        -membrane:normal_search (default= false) - perform membrane normal search with normal_start_angle, normal_delta_angle, and normal_max_angle values (see below).
        -membrane:center_max_delta (default= 5 A) - magnitude of maximum membrane width deviation during membrane center search (Angstroms).
        -membrane:normal_start_angle (default= 10 degrees) - magnitude of starting angle during membrane normal search (degrees).
        -membrane:normal_delta_angle (default= 10 degrees) - magnitude of angle deviation during membrane normal search (degrees).
        -membrane:normal_max_angle (default= 40 degrees) - magnitude of maximum angle deviation during membrane normal search (degrees).

6.  Expected Outputs

    Convert output silent file into pdb file using score application as follows 
    (see an example in scripts/membrane-centroid-score.cmd):

        /Users/vyarovyarovoy/mini/bin/score.macosgccrelease \
        -in:file:native rosetta_inputs/BRD4.pdb \
        -in:file:centroid_rosetta_inputs/ \
        -in:file:silent BRD4_silent.out \
        -in:file:silent_struct_type binary \
        -in:file:spanfile rosetta_inputs/BRD4.span \
        -in:file:lipofile rosetta_inputs/BRD4.lips4 \
        -membrane:no_interpolate_Mpair \
        -membrane:Menv_penalties \
        -score:find_neighbors_3dgrid \
        -score:weights score_membrane.wts \
        -out:nstruct 1 \
        -database /Users/vyarovyarovoy/minirosetta_database \
        -out:output \
        -nstruct 1

    Membrane ab initio application specific score outputs in the output score file are:

        Mpair                         membrane pairwise residue interaction energy
        Menv                          membrane residue environment energy
        Mcbeta                        membrane residue centroid density energy
        Mlipo                         membrane residue lipophilicity energy
        Menv_hbond                    membrane non-helical secondary structure in the hydrophobic layer penalty
        Menv_termini                  membrane N- and C-temini residue in the hydrophobic layer penalty
        Menv_tm_proj                  transmembrane helix projection penalty

Post Processing
---------------

Generate at least 10,000 models and then use rosetta Cluster application to 
identify most frequently sampled conformations.  In general case, at least one 
of top 5-10 clusters will have models with the lowest rmsd to the native 
structure.

Good luck :)!
DARC Demo
=========

Docking Approach using Ray Casting (DARC) is structure-based computational method for carrying out virtual screening by docking small-molecules into protein surface pockets.
In this demo, DARC is used to dock a small molecule in a pocket centered around residue 61 of the protein, E3 ubiquitin-protein ligase Mdm2 (PDB: 4ERF).

Generating input files
----------------------

* Protein files:

  4ERF.pdb was downloaded from The Protein Databank. Remove any water or ligand 
  or other HETATM present in the protein and remove redundant chains. Although 
  Rosetta will build any missing atoms including hydrogens in input protein PDB 
  files on the fly defore it uses them, we need to prepare the protein to 
  generate the electrostatic potential grid. This is one way to dump the 
  protein after building missing atom and added hydrogens:

      $ Rosetta/main/source/bin/score.linuxgccrelease -in:file:s 4ERF.pdb -out:output -no_optH false

  which gives the output protein 4ERF_0001.pdb

* Ligand files:

  Input files are included in the folder inputs, but the following explains how 
  they were obtained or generated. 4ERF.pdb was downloaded from The Protein 
  Databank. The small molecule used for docking was found in ZINC, the database 
  of commercially available compounds at http://zinc.docking.org/. We download 
  ZINC13989607 in mol2 format. The file name is `zinc_13989607.mol2`.

  Since the charges assigned to the atom is same for all conformers, it is fast 
  to assign charges before generating conformers. The charges were added to the 
  compound using the molcharge program from OpenEye (http://www.eyesopen.com/) 
  using the following command:

      $ OpenEye/bin/molcharge -in zinc_13989607.mol2 -out zinc_13989607_charged.mol2 -method am1bccsym

  If you have the compound in 2D SMILES format, then first convert to 3D format and then add charges.
  You can use openbabel or OpenEye program for converting from SMILES format to mol2 format

      $ OpenEye/bin/omega2 -in zinc_13989607_charged.mol2 -out zinc_13989607_conformers.mol2 -maxconfs 100

  then we need to generate parameter files for the compound to use in Rosetta.
  save the names of the compounds in a text file for generating parameters. Next, parameter files for the compounds are generated using batch_molfile_to_params.py which is a python app in rosetta. This command is used:For Ex:

      $ echo zinc_13989607_conformers.mol2 > molfile_list.txt

      $ Rosetta/main/source/src/python/apps/public/batch_molfile_to_params.py -d Rosetta/main/database --script_path=Rosetta/main/source/src/python/apps/public/molfile_to_params.py molfile_list.txt

  The expected output files are:

      params/zinc_13989607_conformers/000_conformers.pdb
      params/zinc_13989607_conformers/000.params

* Other input files:

  To run DARC we need to generate a RAY file for the input protein. To generate 
  this ray-file we need to input the protein in PDB format and specify a target 
  residue at the interface. We can specify more than one residue at the 
  interface. The command to run DARC is as follows:

      $ Rosetta/main/source/bin/make_ray_files.linuxgccrelease -database Rosetta/main/database/ -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only 

  Expected output from this command will be a ray-file named 
  `ray_4ERF_0001_61.txt`.

  To use a center of mass of any residue as origin points for casting rays:

      $ Rosetta/main/source/bin/make_ray_files.linuxgccrelease -database Rosetta/main/database/ -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only -set_origin 5 -origin_res_num 85:A

  Expected output: `ray_4ERF_0001_61.txt`

  To use multiple origin points for casting rays:

      $ Rosetta/main/source/bin/make_ray_files.linuxgccrelease -database Rosetta/main/database/ -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only -set_origin 5 -origin_res_num 85:A -multiple_origin

  Expected output: `ray_4ERF_0001_61,54.txt`

  To use a bound ligand to center the grid:

      $ Rosetta/main/source/bin/make_ray_files.linuxgccrelease -database Rosetta/main/database/ -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only -set_origin 5 -origin_res_num 85:A -bound_ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -lig_grid

  Expected output: `ray_4ERF_0001_61.txt`

  To include electrostatics calculations:

  For electrostatics calculations, first we need to resize the electrostatic 
  potential grid (generated from openeye '4ERF.agd') to match the size of the 
  interface pocket grid. This step can be carried out while generating the ray 
  file.

  To generate the electrostatic potential grid, we can use the examples 
  Listings provided in the zap toolkit manual.

      $ OpenEye/bin/Listing_2 -in 4ERF_0001.pdb -out 4ERF.agd -buffer 2 -grid_spacing 0.5 -epsout 80 -epsin 1

  Expected output: `4ERF.agd`

      $ Rosetta/main/source/bin/make_ray_files.linuxgccrelease -database Rosetta/main/database/ -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61,54 -set_origin 5 -origin_res_num 85:A -multiple_origin -bound_ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -lig_grid -espGrid_file 4ERF.agd

  The output from this command will be a ray-file named 
  `ray_4ERF_0001_61,54.txt` and an electrostatic potential grid file named 
  `DARC_4ERF.agd` which we use as input for running docking using DARC. 

Running DARC
------------

In the second step, we are actually running the docking calculations using the 
pre-generated ray-file. Call DARC for running the docking calculations using 
the pre-generated ray-file and the matched electrostatic potential grid. use 
the flag -num_particles to increase or decrease the sampling. high number of 
particles will take a longer time to finish. For best results, We suggest 
alleast 100 particles and 100 runs for docking single conformer or iterative 
docking (one-by-one) of ligand conformers. For on-the-fly sampling of ligand 
conformers we suggest at least 500 particles and 500 runs. Here we give the 
input ligands for screening against the ray-file, as follows: 

Docking single ligand conformer:

    $  Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd 

Expected output: `DARC_4ERF_0001_LG1.pdb`

Docking multiple ligand conformers:

    $  Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd 

Expected output: `DARC_4ERF_0001_000.pdb`

To run DARC with shape only (without including electrostatics score):

    $  Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file -darc_shape_only

Expected output: `DARC_4ERF_0001_000.pdb`

To search conformers on-the-fly:

    $  Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -search_conformers true

Expected output: `DARC_4ERF_0001_000.pdb`

The output of the DARC run is a docked model of the protein-ligand complex 
named `DARC_4ERF_0001_000.pdb`. To minimize the docked model: We can do the 
fullatom minimization of the DARC models separately in Rosetta or as an 
additional option while running DARC itself. To minimize the DARC models 
immediately after docking we add the flag `-minimize_output_complex`:

    $  Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -minimize_output_complex

Expected output files are:

* DARC docked model complex: `DARC_4ERF_0001_000.pdb`
* DARC docked and minimized model complex: `mini_4ERF_0001_000.pdb`

Other options for DARC:

* `-use_connolly_surface`:  use this flag to make ray files with connolly surface insted of pocketshell generated by pocketgrid
* `-num_particles`:         set the number of particles used durin particle swarm optimization
* `-num_runs`:              set the number of particles used durin particle swarm optimization
* `-missing_point_weight`:  weight for rays misses the pocket
* `-steric_weight`:         weight for rays missesthe ligand
* `-extra_point_weight`:    weight for ligand moves away from the pocket
* `-esp_weight`:            electrostatics weight
* `-pocket_static_grid`:    no autoexpansion of the pocket grid (set ON for DARC)

Building Rosetta
----------------
Build Rosetta from the source directory

    $ cd /Rosetta/main/source/
    $ ./scons.py mode=release bin

To build with GPU enabled:

    $ ./scons.py mode=release extras=opencl bin

Note: Building with opencl results in different binary extension (eg: 
DARC.opencl.linuxgccrelease) run on gpu systems, we need to use the 
'opencl.linuxgccrelease' binary extension.

Example: Commands for sample DARC run for protein MDM2 [PDB:4ERF]
-----------------------------------------------------------------

    $ cd run_dir

Copy 4ERF.pdb, 000_conformers.pdb, and 000.params, 4ERF.agd, DARC_4ERF.agd, 4ERF_0001.pdb, 4ERF_XTL_0001.pdb, 4ERF_XTL.params to run_dir. 
Run DARC with the following command:

* DARC shape only:

        $ Rosetta/main/source/bin/score.linuxgccrelease -in:file:s 4ERF.pdb -out:output -no_optH false
        $ Rosetta/main/source/bin/make_ray_files.linuxgccrelease -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only
        $ Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -darc_shape_only -num_particles 50 -num_runs 50 

  Expected output: `DARC_4ERF_0001_LG1.pdb`

* DARC shape and electrostatics:

        $ Rosetta/main/source/bin/make_ray_files.linuxgccrelease -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -espGrid_file 4ERF.agd
        $ Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 50 -num_runs 50

  Expected output: `DARC_4ERF_0001_LG1.pdb`

* DARC sampling conformers on the fly

        $ Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 500 -num_runs 500 -search conformers true

  Expected output: `DARC_4ERF_0001_000.pdb`

* DARC sampling conformers iteratively one-by-one

        $ Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 500 -num_runs 500 -search conformers false

  Expected output: `DARC_4ERF_0001_000.pdb`

Running DARC in GPU systems
---------------------------

DARC is adapted to run on GPU systems. On average, Running DARC on GPU system 
is 190 fold faster than running on CPU. The results from DARC that runs on CPU 
and GPU systems should be same (with same input and constant seed).

First make sure you built rosetta with the flag `extras=opencl` to run on GPU 
systems

    $ cd Rosetta/main/source
    $ ./scons.py mode=release extras=opencl bin

Then use the opencl executables instead of default executables with the flag '-gpu 1', i.e call DARC.opencl.linuxgccrelease 
instead of DARC.linuxgccrelease. All other options are same as CPU.
Here is an example for GPU command:

    $ Rosetta/main/source/bin/score.opencl.linuxgccrelease -in:file:s 4ERF.pdb -out:output -no_optH false -gpu 1
    $ Rosetta/main/source/bin/make_ray_files.opencl.linuxgccrelease -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only -gpu 1
    $ Rosetta/main/source/bin/DARC.opencl.linuxgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -darc_shape_only -gpu 1 -num_particles 500 -num_runs 500 -search_conformers true

Expected output: `DARC_4ERF_0001_000.pdb`

Running DARC will result in these output files:

* DARC_4ERF_0001_LG1.pdb is the final docked pose for single conformer.
* DARC_4ERF_0001_000.pdb is the final docked pose for multiple conformers. 
  Note:The flag '-minimize_output_complex' can be used to minimized the output 
  model
* mini_4ERF_0001_LG1.pdb is final docked and minimized pose for single 
  conformers.
* mini_4ERF_0001_000.pdb is final docked and minimized pose for multiple 
  conformers.

Analyzing the results
---------------------

The output should look similar to:

    Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 10 -num_runs 10
    core.init: Rosetta version 1e0e3a1462814b2fa93580e9e1b1f7ab81736efc 2014-11-19 15:22:07 -0600 from git@github.com:RosettaCommons/main.git
    core.init: command: /Users/ragul/Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 10 -num_runs 10
    core.init: 'RNG device' seed mode, using '/dev/urandom', seed=590325598 seed_offset=0 real_seed=590325598
    core.init.random: RandomGenerator:init: Normal mode, seed=590325598 RG_type=mt19937
    core.init: Resolved executable path: /Users/ragul/Rosetta/main/source/bin/DARC.linuxgccrelease
    core.init: Looking for database based on location of executable: /Users/ragul/Rosetta/main/database/
    core.chemical.ResidueTypeSet: Finished initializing fa_standard residue type set.  Created 738 residue types
    core.pack.task: Packer task: initialize from command line() 
    origin_space : 7.25
    Reading ligand 4ERF_XTL_0001.pdb
    core.pack.task: Packer task: initialize from command line() 
    Starting PSO
    Time for DARC runs: 0.589119
    SCORES : 4ERF_0001_LG1	Conformer 1 has DARC Score : 4.29993

In this case the DARC score for the final docked pose is 4.29993. This score 
can be use to compare with other ligand that are docked using the same ray 
file. The lower the score the better the match between the protein pocket and 
ligand. `DARC_4ERF_0001_LG1.pdb` is the out file that contains the final 
docked pose for single conformer.

If we use the flag `-minimize_output_complex` the model will be minimized and 
we get the file `mini_4ERF_0001_LG1.pdb` as output.

Design with non-canonical amino acids (NCAA)
============================================

In the first part of this tutorial we show how to make a point mutation in a 
PDB file (1l2y) including two non-canonical amino acids (NVL and HLU). 
In the second part of this tutorial we show how to incorporate novel 
non-canonical amino acid side-chains into rosetta for design.

Designing with NCAA
===================

Below are two command-lines for making point-mutations including the possibility 
of incorporating non-canonical amino acids. In this procedure the backbone is 
kept fixed and the Trp6 residue is allowed to be altered.

* Trp6 → Canonical amino acids plus NVL and HLU

        ~/svn/mini/bin/fixbb.linuxiccrelease -s 1l2y.pdb -use_input_sc -nstruct 1 -database ~/minirosetta_database/ -ex1 -ex2 -overwrite -minimize_sidechains -resfile ncaa_resfile_pluscanon

* Trp6 → NVL or HLU

        ~/svn/mini/bin/fixbb.linuxiccrelease -s 1l2y.pdb -use_input_sc -nstruct 1 -database ~/minirosetta_database/ -ex1 -ex2 -overwrite -minimize_sidechains -resfile ncaa_resfile

Adding new sidechains to Rosetta
================================

Generating an initial sidechain model
-------------------------------------
I usually do this in PyMOL but any molecular editor will work. PyMOL has 
rudimentary but useful editing functionality. The goal is to get the 
structure close to the final structure. Don't stress too much because the QM 
minimization will clean up most things. Make sure you double check chirality. 
To make sure the geometry for the atoms near the ends of the molecule is 
correct you will need to put capping groups on the molecule. I would recommend 
an acetyl (ACE) on the nitro-terminus of the betapeptide(BP) and a n-methyl 
(NME) on the carboxy terminus. This structure (ACE-X-NME) is usally called a 
dipeptide despite the fact that it only has one complete residue.

Pro-tips:
- Within PyMOL go to the "Mouse" drop down and select "3 button editing".
- Within PyMOL type "help edit_keys" in the command box will bring up the instructions for the editor. 
- Within PyMOL type "set valence, 1" in the command box to show single, double, triple bonds.
- Be consistent with atom order. It make future steps easier. I like to put all the atoms for the capping groups first, followed by backbone heavy atoms, side chain heavy atoms, and then then hydrogens.

See the example pdb files in the folder `stage_01_initial_structures` and shown 
below:

    <<<<< ornithine.pdb start >>>>>
    ATOM      1  C   ACE     1       0.984   0.045  -0.578  1.00  0.00           C
    ATOM      2  O   ACE     1       1.815  -0.721  -1.083  1.00  0.00           O
    ATOM      3  CH3 ACE     1      -0.445  -0.396  -0.349  1.00  0.00           C
    ATOM      4 1HH3 ACE     1      -0.468  -1.251   0.313  1.00  0.00           H
    ATOM      5 2HH3 ACE     1      -1.013   0.407   0.098  1.00  0.00           H
    ATOM      6 3HH3 ACE     1      -0.904  -0.668  -1.288  1.00  0.00           H
    ATOM     26  N   NME     3       4.937   2.266   0.684  1.00  0.00           N
    ATOM     27  CH3 NME     3       5.341   3.324   1.592  1.00  0.00           C
    ATOM     28  H   NME     3       5.534   1.511   0.316  1.00  0.00           H
    ATOM     29 1HH3 NME     3       4.689   4.177   1.478  1.00  0.00           H
    ATOM     30 2HH3 NME     3       5.286   2.977   2.614  1.00  0.00           H
    ATOM     31 3HH3 NME     3       6.355   3.627   1.379  1.00  0.00           H
    ATOM      7  N   ALA     2       1.392   1.405  -0.195  1.00  0.00           N
    ATOM      8  CA  ALA     2       2.806   1.518  -0.539  1.00  0.00           C
    ATOM      9  C   ALA     2       3.512   2.478   0.390  1.00  0.00           C
    ATOM     10  O   ALA     2       2.929   3.435   0.912  1.00  0.00           O
    ATOM     11  CB  ALA     2       2.895   1.945  -2.014  1.00  0.00           C
    ATOM     12  C01 ALA     2       2.214   0.926  -2.946  1.00  0.00           C
    ATOM     13  N01 ALA     2       1.680   0.434  -5.295  1.00  0.00           N
    ATOM     14  C02 ALA     2       2.334   1.404  -4.405  1.00  0.00           C
    ATOM     15  H   ALA     2       0.770   2.167   0.251  1.00  0.00           H
    ATOM     16  HA  ALA     2       3.283   0.529  -0.413  1.00  0.00           H
    ATOM     17 2HB  ALA     2       2.409   2.923  -2.195  1.00  0.00           H
    ATOM     18 3HB  ALA     2       3.943   2.035  -2.355  1.00  0.00           H
    ATOM     19  H01 ALA     2       2.698  -0.045  -2.840  1.00  0.00           H
    ATOM     20  H02 ALA     2       1.461  -0.405  -4.777  1.00  0.00           H
    ATOM     21  H03 ALA     2       1.161   0.837  -2.679  1.00  0.00           H
    ATOM     22  H04 ALA     2       2.300   0.206  -6.059  1.00  0.00           H
    ATOM     23  H05 ALA     2       3.387   1.492  -4.674  1.00  0.00           H
    ATOM     24  H06 ALA     2       1.850   2.375  -4.509  1.00  0.00           H
    ATOM     25  H07 ALA     2       0.828   0.835  -5.661  1.00  0.00           H
    END
    <<<<< ornithine.pdb end >>>>>

Making GAUSSIAN input
---------------------
We need to minimize the initial structure we made in PyMOL to get a good set of 
ideal bond lengths and angles. To do this we will take the coordinates from the 
pdb file and put them into a gaussian input file. Like the one shown below and 
in stage_02_gaussian_input. A discussion of the complete gaussian input 
structure is beyond the scope of this document. Documentation for gaussian can 
be found at http://www.gaussian.com/.

In short the input is as follows:

- line 1 sets the path for the checkpoint file
- line 2 describes the level of theory, options, and convergence criteria. You may need to change the basis set to something smaller if you are using stuff bellow the 4th line of the periodic table. 
- lines 3-5 are comments
- line 6 is the charge and multiplicity (this is usually 0 1 but ornithine is charged)
- line 7-37 are the elemental type and xyz coordinates of the atoms from the pdb file in the same order as the pdb file
- line 38 is blank
- line 39-40 is the modredundant input that says that we want to keep the torsion formed by atoms 1 13 14 and 15 fixed at 150.00 degrees
- line 41 is blank.

Running Gaussian is simple but the minimizations take a long time (a few hours 
per structure). The command bellow will run gaussian on the input file in the 
`stage_02` folder and put the output in the `stage_03` folder.

    $ g03 stage_02_gaussian_input/ornithine.com stage_03_gaussian_output/ornithine.log

There is a new version of Gaussian called Gaussian09. If you are using 
this version, the commands and options might be slightly different. The input 
below will probably work with Gaussian09 but was written for and tested with 
Gaussian03.

The Gaussian output is very verbose. It's not included in this document but it is 
in the demo folder.

    <<<<< ornithine.com start >>>>>
    %Chk=stage02_ace_nme_res_ordered_pdbs/ornithine.chk
    # HF/6-31G(d) Opt=ModRedundant SCF=Tight Test 

    scan rotamers

    1  1
    C        0.984   0.045  -0.578
    O        1.815  -0.721  -1.083
    C       -0.445  -0.396  -0.349
    H       -0.468  -1.251   0.313
    H       -1.013   0.407   0.098
    H       -0.904  -0.668  -1.288
    N        4.937   2.266   0.684
    C        5.341   3.324   1.592
    H        5.534   1.511   0.316
    H        4.689   4.177   1.478
    H        5.286   2.977   2.614
    H        6.355   3.627   1.379
    N        1.392   1.405  -0.195
    C        2.806   1.518  -0.539
    C        3.512   2.478   0.390
    O        2.929   3.435   0.912
    C        2.895   1.945  -2.014
    C        2.214   0.926  -2.946
    N        1.680   0.434  -5.295
    C        2.334   1.404  -4.405
    H        0.770   2.167   0.251
    H        3.283   0.529  -0.413
    H        2.409   2.923  -2.195
    H        3.943   2.035  -2.355
    H        2.698  -0.045  -2.840
    H        1.461  -0.405  -4.777
    H        1.161   0.837  -2.679
    H        2.300   0.206  -6.059
    H        3.387   1.492  -4.674
    H        1.850   2.375  -4.509
    H        0.828   0.835  -5.661

    1 13 14 15  -150.00 F
    13 14 15 7  150.00 F
    <<<<< ornithine.com end >>>>>

Converting GAUSSIAN output to molfile
-------------------------------------
The program bable we built in the first step can convert the gaussian output to 
a molfile.

Bable says it can handle output form Gaussian09. I am not sure when it was 
added, the current version is 2.3.0 and the modified version I have included is 
based on version 2.2.0. The changes to 2.2.0 were not extensive and could 
probably be ported to 2.3.0.

The command bellow will convert the gaussian output to molfile format for the 
peptide and peptoid examples:

    ```
    $ openbabel/install/bin/babel -i g03 stage_03_gaussian_output/ornithine.log -o mol stage_04_molfile/ornithine.mol
    $ openbabel/install/bin/babel -i g03 stage_03_gaussian_output/amino2.log -o mol stage_04_molfile/amino2.mol
    ```

    <<<<< ornithine.mol start >>>>>
     OpenBabel01101117083D

     31 30  0  0  0  0  0  0  0  0999 V2000
        0.4866    2.2296   -0.2354 C   0  0  0  0  0
        1.2081    1.9315   -1.1549 O   0  0  0  0  0
        0.4369    3.6268    0.3352 C   0  0  0  0  0
       -0.3187    4.1933   -0.1996 H   0  0  0  0  0
        0.1852    3.6362    1.3888 H   0  0  0  0  0
        1.3920    4.1075    0.1799 H   0  0  0  0  0
       -2.6920   -1.1572   -0.6650 N   0  0  0  0  0
       -4.0550   -1.5847   -0.3922 C   0  0  0  0  0
       -2.3500   -1.2554   -1.5934 H   0  0  0  0  0
       -4.1238   -1.9440    0.6231 H   0  0  0  0  0
       -4.7601   -0.7727   -0.5220 H   0  0  0  0  0
       -4.3066   -2.3874   -1.0709 H   0  0  0  0  0
       -0.3415    1.3293    0.3431 N   0  0  0  0  0
       -0.6023    0.0286   -0.2270 C   0  0  0  0  0
       -2.0235   -0.3660    0.1877 C   0  0  0  0  0
       -2.4547   -0.0133    1.2526 O   0  0  0  0  0
        0.3558   -1.0672    0.2829 C   0  0  0  0  0
        1.8184   -0.8036   -0.0848 C   0  0  0  0  0
        4.1559   -1.6146    0.0004 N   0  0  0  0  0
        2.7190   -1.9434    0.3678 C   0  0  0  0  0
       -0.9580    1.6065    1.0765 H   0  0  0  0  0
       -0.5197    0.1005   -1.3044 H   0  0  0  0  0
        0.2486   -1.1369    1.3608 H   0  0  0  0  0
        0.0400   -2.0197   -0.1352 H   0  0  0  0  0
        1.9103   -0.6656   -1.1566 H   0  0  0  0  0
        4.2542   -1.4797   -0.9970 H   0  0  0  0  0
        2.1487    0.1188    0.3755 H   0  0  0  0  0
        4.7999   -2.3436    0.2758 H   0  0  0  0  0
        2.4931   -2.8811   -0.1191 H   0  0  0  0  0
        2.7087   -2.0883    1.4385 H   0  0  0  0  0
        4.4568   -0.7571    0.4437 H   0  0  0  0  0
      2  1  2  0  0  0
      1  3  1  0  0  0
      1 13  1  0  0  0
      4  3  1  0  0  0
      6  3  1  0  0  0
      3  5  1  0  0  0
      9  7  1  0  0  0
      7  8  1  0  0  0
      7 15  1  0  0  0
     12  8  1  0  0  0
     11  8  1  0  0  0
      8 10  1  0  0  0
     14 13  1  0  0  0
     13 21  1  0  0  0
     22 14  1  0  0  0
     14 15  1  0  0  0
     14 17  1  0  0  0
     15 16  2  0  0  0
     24 17  1  0  0  0
     18 17  1  0  0  0
     17 23  1  0  0  0
     25 18  1  0  0  0
     18 20  1  0  0  0
     18 27  1  0  0  0
     26 19  1  0  0  0
     19 28  1  0  0  0
     19 20  1  0  0  0
     19 31  1  0  0  0
     29 20  1  0  0  0
     20 30  1  0  0  0
    M  END
    <<<<< ornithine.mol end >>>>>

Modifying the molfiles
----------------------
The `molfile2params_polymer.py` script requires some additional data to be added 
to the end of the molfile. This data is specified at the end of the file after 
the bond information. It is a list of variable names and then a list of values. 
The variable are described bellow. 

* ROOT: Single numerical value. Atom number (acording to the order above). 
  Where the atom tree is rooted for this residue type. Should probably be the 
  nitrogen of the central residue.

* POLY_N_BB, POLY_CA_BB, POLY_C_BB, POLY_CO_BB: Single numerical value. The 
  backbone nitrogen, alpha-carbon, carbonyl-carbon and carbonyl-oxygen. These 
  get special rosetta atom types and so are listed here special. Note: You will 
  need to add some additional code to work with your additional backbone atom.

* POLY_IGNORE: List of numerical values.  Atom number (acording to the order in 
  the molfile). These are the atoms for the capping groups with the exception of 
  the upper and lower connect atoms. They will not be listed in the atoms and 
  bonds in the params file but are used in determining the atom types.

* POLY_UPPER, POLY_LOWER: Single numerical value. Atom number (according to the 
  order in the molfile). These are the atoms in the capping groups that connect 
  to the residue. They will not be listed in the atoms and bonds in the params 
  file but are used in determining the atom types and they are listed in the 
  internal coordinate section.

* POLY_CHG: Single numerical value. Overall charge on the residue. 

* POLY_PROPERTIES: List of alpha-numerical values. These get used by rosetta at 
  various places in the program. You can say something like:
  ```
  if (pose.residue(10).type().is_protein() ) { /* do someting */ }
  ```
  You will want to create a new property called "BETAPEPTIDE" and corresponding 
  functions to the residue type class. 

* END: The end of the file.

Note that there are 2 spaces between the "M" and the variable name. If you have to 
make multiple residue types (which it sounds like you will), keeping the atoms 
in the same order makes assigning all these numbers easier. 

Additional info for peptide ornithine:

    M  ROOT 13
    M  POLY_N_BB 13
    M  POLY_CA_BB 14
    M  POLY_C_BB 15
    M  POLY_O_BB 16
    M  POLY_IGNORE 2 3 4 5 6 8 9 10 11 12
    M  POLY_UPPER 7
    M  POLY_LOWER 1
    M  POLY_CHG 1
    M  POLY_PROPERTIES PROTEIN POLAR CHARGED
    M  END

Additional info for amino2:

    M  ROOT 1
    M  POLY_N_BB 1
    M  POLY_CA_BB 2
    M  POLY_C_BB 3
    M  POLY_O_BB 4
    M  POLY_IGNORE 9 10 22 23 24 27 28 29 30 31 32 33 34
    M  POLY_UPPER 26
    M  POLY_LOWER 7
    M  POLY_CHG 1
    M  POLY_PROPERTIES PEPTOID POLAR CHARGED
    M  END

Running molfile2params.py
-------------------------
Finally, the `molfile2params_polymer.py` script will convert the modified molfile 
to a params file. The commands bellow work for the example files:

    $ python molfile_to_params_polymer.py --clobber --polymer --no-pdb --name C40 -k ornithine.kin stage_05_modified_molfile/ornithine.mol
    $ python molfile_to_params_polymer.py --clobber --polymer --peptoid --no-pdb --name P01 -k amino2.kin stage_05_modified_molfile/amino2.mol

There is additional tweaking that needs to happen to the params files to make 
them work. Compare these to the ones in the database for reference. 

The `molfile2params.py` script can produce a kinemage file. This is super handy 
as it lets you check how Rosetta will build the atom tree, and the atom type 
assignments, and other stuff. You can open kinemage files using the KiNG 
program from the Richardson lab at Duke 
(http://kinemage.biochem.duke.edu/software/king.php).

Using your shiny new params file
================================
Rosetta params files live in the database at database/chemical/residue_type_sets/(fa_standard)/..., but you can just pass your new parameters in on command line via `-extra_res_fa` or `-extra_res` (the latter for centroid).

Rotamer library design
======================

For questions on backbone-dependent rotamer library design please contact 
either Doug Renfrew (dougrenfrew at gmail dot com) or Tim Craven (twc254 at nyu dot edu).

Docking with constraints
========================

In this demo, a peptide (chain B) is docked onto a protein (chain A) with 
user-defined constraints. The purpose of this demo is to illustrate the use of 
constraint files and as many different types of constraints as possible. The 
input files are based on PDB ID (2o2m) which is a "Crystal Structure Of Human 
Gαi1 Bound To The Goloco Motif Of Rgs14". The protein is chain A. The peptide 
is chain B. 

Protocol
--------

A general overview of the protocol is available at:  
http://www.rosettacommons.org/manuals/archive/rosetta3.2.1_user_guide/app_docking.html

A write-up of the file format for Constraints is available at:  
http://www.rosettacommons.org/manuals/archive/rosetta3.2.1_user_guide/constraint_file.html

High level overview:

1. Docking Prepack
    - Useful to optimize all of the sidechains on both partners
    - http://graylab.jhu.edu/Rosetta.Developer.Documentation/all_else/de/d69/docking_prepack_protocol.html

2. Docking with constraints. In general, constraints are typically defined for CA atoms of the partner sidechains.
    - Docking will honor AtomPairContraint, AmbiguousConstraint, SiteContraint. 
      For additional information on contraints, see "Constraints" section 
      below. 

3. Analysis
    - The expected output of the docking run is (1) a decoy with both docked 
      partners and (2) a score file with score breakdown for each decoy 
    - Typical analysis would be to sort the generated score file with a command 
      line i.e. 'sort -nk 28` where the integer identifies the column number 
      for the quality to sort on.

Constraints
-----------

A good overview of constraints is available in the documentation:  
http://www.rosettacommons.org/manuals/archive/rosetta3.2.1_user_guide/constraint_file.html

The three constraints which are honored by the Docking Protocol are: 

* AtomPairConstraint
* AmbiguousConstraint
* SiteConstraint.

AtomPair and AmbiguousConstraint are described the above documentation for 
Rosetta 3.2.1.  A SiteConstraint allows you to specify that a particular 
residue should be in contact with a particular chain. An example of a 
SiteConstraint is:

    SiteConstraint CA 4A D FLAT_HARMONIC 0 1 5

This will add a FLAT_HARMONIC potential with the parameters 0 1 5 (recommended) 
around the distance between the CA of residue 4 (PDB numbering) on chain A and 
the closest CA on chain D to the ScoreFunction. 

# Oop Design

## Author
written by Kevin Drew (Bonneau Lab), kdrew@nyu.edu

## General Description
This demo shows how to run the oop\_design application.  An oligooxopiperazine (oop) is a helical mimetic scaffold used for inhibiting protein interactions. The demo shows the design of an oop inhibitor for the MDM2-P53 protein interaction.

## Algorithm
1. Pertubation phase: rigid body movement of oop wrt target, oop small moves to oop ring conformation, oop puck move to change oop ring pucker, small moves to ring linkers
2. Design phase: design user specified residues on oop scaffold and minimize
3. Repeat 10x

## Command
```
oop_design<.exe> -database <path to your database> @input/flags
```

## Input files
```
./input/flags = user specified options
./input/mdm2_oopAAAA.pdb = input structure where target is chain 1 and oop is chain 2
```

## Options
```
-oop_design_positions = residues on oop to design (numbering is relative to oop, for example 3 is the third residue on oop), default repacks with no design
-pert_num = number of pertubations during pertubation phase, default 10, production 100
-design_loop_num = number of pertubation + design cycles, default 10, production 10
-nstruct = production 1000
```

## Post Processing
Similar to other multi chain design protocols, the ddG is computed and is a good indicator of a good design.  First sort by total score, take top 5 percent and then sort by REPACK_ENERGY_DIFF.

## Limitations: 
This app is inflexible to adjusting monte carlo temperatures, score functions, degree of rigid body pertubations, designing noncanonical amino acids, etc. The app also requires the oop is close to a plausible binding mode with respect to the target. 

# Prepare PDBs for Rosetta

## Introduction
For the most part Rosetta does a pretty good job at reading
structures, however with unusual chemical modifications, sometimes
Rosetta runs into problems. Here are some of the problems with pdb
files that we have encountered that currently cause problems with
Rosetta:


## 0: Check the PDB
If you run into problems with processing a structure, it would
be a good idea to go in and look the text of the pdb file. 

## 1: Necessary Unrecognized Residue
If there is an unrecognized residue,
  1. Go look into the structure and see what it is and what's going
  on. Some structures in the PDB may be C-alpha only like, 2rdo.pdb or
  3htc.pdb.

  2. If it is a HET group that you do not need, such as a non-covalent
  small molecule, like a NAD, HEME group or glycerol that is included
  with the structure, then you can include -ignore_unrecognized_res
  option to the command line.

  3. If there is a covalently linked ligand or modified amino acid or
  base, either there is a parameter file in the
  database/chemical/residue_type_sets and it is looked up
  correctly. If it is important, (for instance if it is part of the
  covalent chain like the CSP residue in 1qu9.pdb) but not in the
  database, it is possible to convert a MDL, MOLfile, or SDF file to a
  Rosetta residue topology file.
    ```
     src/python/apps/public/molfile_to_params.py --help 
    ```
  Please see documenation on preparing ligands for using
  molfile_to_params correctly.

## 2: Unnecessary Unrecognized Residue
If there are unrecognized residues but you don't actually care about them,

  1. Manually remove the HETATM atom records and checking in PyMOL
  that the structure is still reasonable (no broken chains etc.)

  2. Use a script to clean the PDB appropriately. For example
  tools/protein_util/scripts/clean_pdb.py can handle some basic
  cleanup proceedures. This script is a little hacky so use with
  caution.

  3. Run with `-ignore_unrecognized_res` or just `-ignore_water`. Simply
  ignoring unrecognized residues can coding assumptions that are not
  rigerously checked (i.e. expose bugs in Rosetta). So be prepared
  Rosetta to fail in unusual ways: 
     1. if the first residue of a chain is ignored (e.g. with 1mhm.pdb)
    then Rosetta may fail with this error:

        ```
        ERROR: ( anchor_rsd.is_polymer() && !anchor_rsd.is_upper_terminus() ) && ( new_rsd.is_polymer() && !new_rsd.is_lower_terminus() )
        ERROR:: Exit from: src/core/conformation/Conformation.cc line: 428
        ```

    2. if the residue type three letter code mis identified, such as
    SUC (sucrose in 3KB8.pdb) for the SUCK residue type, used in a
    score term to reduce buried cavities.

        ```
        core.io.pdb.file_data: [ WARNING ] discarding 23 atoms at position 207 in file /home/momeara/scr/data/pdb/kb/pdb3kb8.ent.gz. Best match rsd_type:  SUCK
        core.io.pdb.file_data: [ WARNING ] discarding 23 atoms at position 750 in file /home/momeara/scr/data/pdb/kb/pdb3kb8.ent.gz. Best match rsd_type:  SUCK
        core.io.pdb.file_data: [ WARNING ] discarding 23 atoms at position 980 in file /home/momeara/scr/data/pdb/kb/pdb3kb8.ent.gz. Best match rsd_type:  SUCK
        core.conformation.Conformation: [ WARNING ] missing heavyatom:  OXT on residue SER_p:CtermProteinFull 201
        core.conformation.Conformation: [ WARNING ] missing heavyatom: ORIG on residue SUCK 202
        ...
        ERROR: too many tries in fill_missing_atoms!
        ERROR:: Exit from: src/core/conformation/Conformation.cc line: 2590
        ```

## 3: If You Care About Hydrogens
If you care about the coordinates of hydrogens in the PDB file, you
should be careful that they are named correctly. Rosetta does not use
the most uptodate hydrogen atom naming convention. If the hydrogen
atom names are not recognized, then they will be stripped off and
replaced automatically by the -optH protocol when the pdb file is
processed.

To make sure H-atoms are named correctly you can use the
tools/convert_hatom_names.py script.

```
   python tools/convert_hatm_names.py --data_dir input_files --output_dir prepared_input_files
```

# Relax Around Chemically Bound Ligand

# Metadata
This document was last edited on January 4, 2012 by Ron Jacak.  The document was originally written by Andrew Leaver-Fay, and completed, expanded and modified for Doxygen by Ron Jacak.

# Quick Guide
1. Create parameter files for the ligand, and modified residue (if one doesn't already exist)

2. Replace ligand HETATM lines with Rosetta ATOM lines, and change modified residue lines to proper type

3. Make constraints file

4. Run relax protocol

# Introduction

The purpose of this demo is to relax a protein with a chemically attached small molecule. Currently, one example on the LOV2 domain is provided. In time, this README will be expanded to include a second example on a glycosylated protein.

Rosetta requires a "params" file for each of the residue types that are used in a trajectory. These are already available for all of the amino acids, nucleic acids, and most of the metals.  These parameter files live in minirosetta_database/chemical/residue_type_sets/fa_standard/residue_types/. Before relax can be run on a liganded structure, a parameter file for that ligand needs to be created.

The "LOV2" domain binds a flavin and, when exposed to blue light, forms a thiol bond at CYS 450 to this flavin.  One question that would be nice to answer is what effect this chemical bond has on the stability of the LOV2 domain, since it is known that after this bond forms, the C-terminal
"J-alpha" helix unfolds.  To address this with Rosetta, we would like to run relax on the starting structure while preserving this chemical bond.  The crystal structure 2V0W was made by crystalizing the LOV2 domain in the dark, and then exposing the crystals to light.  The thiol bond between CYS 450 and the flavin is visible in the structure. 

# Creating a Ligand Parameter File

We downloaded the .sdf file for the flavin attached to LOV2 from the PDB (http://www.pdb.org/pdb/explore/explore.do?structureId=2V0W), Ligands_noHydrogens_withMissing_1_Instances.sdf (an awful name for this file). Also, the file doesn't contain hydrogens, which we do want. So, first, we need to add hydrogens to this file.  The program REDUCE (which can be obtained from the Richardson lab website) can be used to add hydrogens.  To run REDUCE, you have to put a copy (or symlink) of the reduce_het_dict.txt file in the /usr/local directory.

```
$ cd /usr/local
$ sudo ln -s /Users/andrew/reduce/reduce_het_dict.txt
$ ~/reduce/reduce.3.03.061020.macosx.i386 2V0W.pdb > blah.pdb

$ grep FMN blah.pdb | tail -n 9
HETATM    0 3HM8 FMN A1547      15.201   6.333   6.801  1.00 11.04           H   new
HETATM    0 3HM7 FMN A1547      13.095   4.009   5.709  1.00 11.55           H   new
HETATM    0 2HM8 FMN A1547      14.999   5.775   8.496  1.00 11.04           H   new
HETATM    0 2HM7 FMN A1547      12.825   3.407   7.380  1.00 11.55           H   new
HETATM    0 1HM8 FMN A1547      16.624   6.243   7.893  1.00 11.04           H   new
HETATM    0 1HM7 FMN A1547      13.426   5.078   7.114  1.00 11.55           H   new
HETATM    0  HN3 FMN A1547      19.774  -2.545   4.778  1.00 12.52           H   new
HETATM    0  H9  FMN A1547      18.030   4.577   7.483  1.00 12.40           H   new
HETATM    0  H6  FMN A1547      14.354   1.560   6.108  1.00 12.59           H   new
```

Now open up the PDB file of the ligand in pymol and have PyMOL write it out as a .mol file.
```
$ grep FMN blah.pdb > FMN_w_h.pdb
$ pymol FMN_w_h.pdb
```
--> save molecule

--> save as one file

--> file format: MOL

produces FMN_w_h.mol

We can convert the MOL file to a .params file with a script named molfile_to_params.py (located in ~/rosetta_source/src/python/apps/public/). This script can be run with the "-h" flag to list all the flags that are applicable.
```
$ ~/rosetta_source/src/python/apps/public/molfile_to_params.py --name=FMN FMN_w_h.mol
```

The --name flag specifies what 3-letter name will used to identify this ligand in PDB files.  We chose "FMN" which stands for flavin mononucleotide.  The python script has the possibility of mis-assigning atom types, so we opened up the FNM.params file to look and see that we got the chemistry right.

Example output from a molfile_to_params.py run:
```
Centering ligands at (  19.791,    2.871,    7.017)
Atom names contain duplications -- renaming all atoms.
WARNING:  atom  P1  has valence > 4
WARNING:  structure contains double bonds but no aromatic bonds
  Aromatic bonds must be identified explicitly --
  alternating single/double bonds (Kekule structure) won't cut it.
  This warning does not apply to you if your molecule really isn't aromatic.
Total naive charge -7.250, desired charge 0.000, offsetting all atoms by 0.234
WARNING: fragment 1 has 31 total atoms including H; protein residues have 7 - 24 (DNA: 33)
WARNING: fragment 1 has 31 non-H atoms; protein residues have 4 - 14 (DNA: 22)
WARNING: fragment 1 has 7 rotatable bonds; protein residues have 0 - 4
Average 31.0 atoms (31.0 non-H atoms) per fragment
(Proteins average 15.5 atoms (7.8 non-H atoms) per residue)
WARNING:  no root atom specified, using auto-selected NBR atom instead.
Wrote params file FNM.params
Wrote PDB file FNM_0001.pdb
```

The .params file looks good, and the PDB file FMN_0001.pdb looks good. The other thing we have to do is to tell Rosetta that there is a chemical bond between this flavin and some other residue.  A chemical connection must be placed between the flavin atom C8 and SG atom of the modified CYS
residue. To do this, add a "CONNECT C8" line to the FMN params file (anywhere after the ATOM records is fine).

Finally, we need to add an "ICOOR" line for this connection at the bottom of the file (and to get this, we'll measure the distance to the CYS SG in the 2V0W pdb). The format for this line is:
```
ICOOR_INTERNAL CONN1 <torsion> <180 - N3:C8:SG angle> <C8:SG distance> <parent atom name> <grandparent name> <great-grandparent name>
```

The dihedral isn't particularly useful, but we do need to know the C8-SG distance, C8-SG-CB angle, and the last angle from the "parent atom" of C8 to C8 to SG. we can find the parent atom in the "ICOOR" line for C8 in the .params file.
```
ICOOR_INTERNAL C8 27.068913 64.302080 1.415450 N3 C9 C16
```

This says "N3" is the parent, C9 is the grandparent, and C16 is the great-grandparent.  So we need the N3-C8-SG angle.
```
C8-SG: 1.9 A
C8-SG-CB: 111.2 degrees
N3-C8-SG: 110.7 degrees
```

We'll need some of these parameters for the ICOOR, and some of them later for the constraints we'll use to enforce good geometry for this chemical bond (Rosetta does not do that for you!).

Thus, the line to be added should be the following:<br>
```
ICOOR_INTERNAL  CONN1  180.000000   69.300000    1.900000   C8    N3    C9
```

The name "CONN1" says that this is the location for the atom that is covalently connected to this residue. The number is taken as the order in which the connection was listed in the .params file, and in this case, there is only one CONNECT connection. This modified .params file is FMN_modded.params.

# Using the Newly Created Ligand and Modified Residue Parameter Files

The next step is to swap the newly generated FMN with the existing FMN residue: the point is to use the newly generated names for the atoms.  Atom names have to agree (perfectly!) between the input PDB file and the .params file.  Careful: the .params file has a strange format. The
atom names given on the ATOM lines are column formatted, so don't add extra whitespace. The rest of the line (not the names) is whitespace delimited.  Paste the contents of the FMN_0001.pdb file to the bottom of the original 2V0W file, and remove the HETATM lines using the original
flavin molecule atom names. To add the newly created ligand residue to the database of residue types Rosetta uses, use the flag -extra_res_fa on the command line:
```
-extra_res_fa FMN_modded.params
```

It is also necessary to change the CYS which is forming the chemical bond with the flavin in the 2V0W pdb file to a modified CYS conformation. The version of CYS in Rosetta which forms a chemical bond to something else is named CYX. So rename residue 450 to CYX from CYS. This residue
type is defined in rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/sidechain_conjugation/CYX.params. It is not automatically read in when Rosetta loads: its entry in the rosetta_database/chemical/residue_type_sets/fa_standard/residue_types.txt file is commented
out, so we also have to uncomment this line. ALTERNATIVELY: we can explicitly add this file to the list of files that Rosetta loads by including it on the command line as we did above. The flag for this would be:
```
-extra_res_fa /path/to/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/sidechain_conjugation/CYX.params
```

# Making a constraints file

Now we need to generate a set of constraints to make Rosetta preserve the bond geometry between the CYS and C8.  Constraint files rely on Rosetta numbering which starts at 1 and does not use the PDB numbering. The first residue in this PDB is 401, so CYS 450 will be Rosetta-residue 50. The last residue in the protein is 546, so the ligand (i.e. last) will be Rosetta-residue 147. The distance constraint is given by this line:
```
AtomPair SG 50 C8 147 HARMONIC 1.9 0.01
```

where HARMONIC is the functional form, 1.9 is the x0 (the ideal value, in Angstroms), and 0.01 is the standard deviation, i.e. the inverse of the spring constant. The two angle constraints are given by these lines:
```
Angle SG 50 C8 147 N3 147 HARMONIC 1.93207948 0.034906585
Angle CB 50 SG 50 C8 147 HARMONIC 1.94080613 0.034906585
```

where for the first constraint, 1.93207948 represents the ideal angle (x0) of 110.7 degrees, given in radians, and 0.034906585 represents a standard deviation of 2 degrees, given in radians.

** From google: 110.7 * pi / 180 = 1.93207948<br>
** From google:   2.0 * pi / 180 = 0.034906585<br>

# Running the relax protocol

To run the relax protocol, we need to pass in a PDB with the correct FMN and CYX lines, the parameter files for the modified CYS and the FMN, and the constraint file.  We also need to activate constraints during scorefunction evaluation, which can be done using the score:weights flag
on the command line.  A complete command line for the protocol would be as follows:
```
~/rosettabin/relax.linuxgccrelease -database ~/rosettadb/ -s 2V0W.pdb -in:file:fullatom -extra_res_fa FMN_modded.params -extra_res_fa ~/rosettadb/chemical/residue_type_sets/fa_standard/residue_types/sidechain_conjugation/CYX.params -overwrite -cst_fa_file chemical_bond.cst -score:weights score12+constraints.wts -mute basic core.init core.scoring
```

Truncated output from the command:
```
protocols.relax.FastRelax: ================== Using default script ==================
protocols.jd2.PDBJobInputter: Instantiate PDBJobInputter
protocols.jd2.PDBJobInputter: PDBJobInputter::fill_jobs
protocols.jd2.PDBJobInputter: pushing 2V0W.pdb nstruct index 1
protocols.jd2.PDBJobInputter: PDBJobInputter::pose_from_job
core.chemical.ResidueTypeSet: Finished initializing fa_standard residue type set.  Created 6480 residue types
core.conformation.Conformation: Connecting residues: 50 ( CYX ) and 147 ( FMN ) at atoms  SG  and  C8 
core.conformation.Conformation:  with mututal distances: 2.63043 and 2.09338
core.import_pose.import_pose: Can't find a chemical connection for residue 147 FMN
core.pack.task: Packer task: initialize from command line() 
protocols.jd2.PDBJobInputter: filling pose from PDB 2V0W.pdb
core.io.constraints: read constraints from chemical_bond.cst
HARMONIC 1.93208 0.0349066

HARMONIC 1.94081 0.0349066

core.pack.task: Packer task: initialize from command line() 
core.pack.dunbrack: Dunbrack library took 0.028403 seconds to load from binary
protocols.relax.FastRelax: CMD: repeat  -80.6309  0  0  0.44
core.pack.interaction_graph.interaction_graph_factory: Instantiating DensePDInteractionGraph
core.pack.pack_rotamers: built 1514 rotamers at 147 positions.
core.pack.pack_rotamers: IG: 1208136 bytes
core.pack.pack_rotamers: pack_rotamers_run(): simulated annealing took 0.080246 seconds
```

The relax protocol should output a structure named 2V0W_0001.pdb which has lower energy than the starting structure.

Idealize Demo
=============

Compatible with the Rosetta SVN code repository after the Revision 51602

    ~/rosetta/rosetta_source/bin/rosetta_scripts.linuxgccrelease -database ~/rosetta/rosetta_database -s ./input.pdb -parser:protocol ./input.xml -nstruct 100 -ex1 -ex2aro

Principles for designing ideal protein structures Nobuyau Koga, Rie 
Tatumi-Koga, Gaohua Liu, Rong Xiao, Thoma B. Acton, Gaetano T. Montelione, and 
David Baker Nature 2012, Supplementary Data 

MOHCA-Seq Modeling Demo
=======================

Author: Clarence Cheng (cyucheng at stanford dot edu)

---

The topic of this demo is Fragment Assembly of RNA with Full-Atom Refinement 
(FARFAR), guided by pseudoenergy constraints from Multiplexed •OH Cleavage 
Analysis by paired-end Sequencing (MOHCA-seq) and secondary structure 
information from one-dimensional chemical mapping and mutate-and-map (M2).

This demo was taken from the appendix of the Methods in Enzymology 
chapter: "Modeling complex RNA tertiary folds with Rosetta".  The abstract for 
that chapter is included below:

Reliable modeling of RNA tertiary structures is key to both understanding these 
structures’ roles in complex biological machines and to eventually facilitating 
their design for molecular computing and robotics. In recent years, a concerted 
effort to improve computational prediction of RNA structure through the 
RNA-Puzzles blind prediction trials has accelerated advances in the field. 
Among other approaches, the versatile and expanding Rosetta molecular modeling 
software now permits modeling of RNAs in the 100 to 300 nucleotide size range 
at consistent sub-helical (~1 nanometer) resolution. Our lab’s current 
state-of-the-art methods for RNAs in this size range involve Fragment Assembly 
of RNA with Full-Atom Refinement (FARFAR), which optimizes RNA conformations in 
the context of a physically realistic energy function, as well as hybrid 
techniques that leverage experimental data to inform computational modeling. In 
this chapter, we give a practical guide to our current workflow for modeling 
RNA three-dimensional structures using FARFAR, including strategies for using 
data from multidimensional chemical mapping experiments to focus sampling and 
select accurate conformations.

Running the demo
----------------

Documentation for setting up Rosetta and RNA tools:
* https://www.rosettacommons.org/docs/latest/Build-Documentation.html
* https://www.rosettacommons.org/docs/latest/RNA-tools.html 

Step-by-step instructions for running the demo:

1.  Example FASTA file:

        >3P49_RNA.pdb
        ggauaugaggagagauuucauuuuaaugaaacaccgaagaaguaaaucuuucagguaaaaaggacucauauuggacgaaccucuggagagcuuaucuaagagauaacaccgaaggagcaaagcuaauuuuagccuaaacucucagguaaaaggacggag

    The RNA sequence must be lowercase.

2.  Example secondary structure file:

        .((((((((......((((((....)))))).(((...((((.....))))..)))........))))))))........(((((......((((((...)))))).(((...((((....((((....)))).....))))..))).......)))))

3.  Generate command lines for helix pre-assembly:

        helix_preassemble_setup.py -secstruct [secondary structure file] -fasta [FASTA file]

4.  Example command line for helix pre-assembly:

        rna_denovo -nstruct 100 -params_file helix0.params -fasta helix0.fasta  -out:file:silent helix0.out -include_neighbor_base_stacks  -minimize_rna true -rna::corrected_geo -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -score:weights stepwise/rna/rna_helix -cycles 1000 -output_res_num 2-9 65-72

5.  Run command lines for helix pre-assembly (local):

        source CMDLINES

6.  Prepare native/reference structure for Rosetta, if available:

        make_rna_rosetta_ready.py 3P49.pdb

    Outputs reformatted native model as `3p49_RNA.pdb`, to be input to 
    README_SETUP. In the glycine riboswitch example presented here, the 3P49 
    crystal structure includes a protein-binding loop that is not part of the 
    construct used for de novo modeling. Command lines [5] through [14] show 
    how to replace the extraneous residues with a tetraloop matching the 
    experimentally probed construct using a short FARFAR modeling run.

7.  Cut out a segment of a model:

        pdbslice.py [3p49_RNA.pdb] -subset [1-21 36-169] [slice_]

    The first input is the model from which you want to excise regions of 
    interest. The second input is the range of nucleotides that you want to 
    keep in your model. The third input is the prefix that will be added to the 
    beginning of the input model’s filename. Here, the protein-binding loop is 
    excised by specifying the range of residues given in the command line.

8.  Renumber a PDB:

        renumber_pdb_in_place.py [slice_3p49_RNA.pdb] [1-21 26-159]

    The first input is the PDB file to be renumbered and the second input is 
    the desired final ranges of sequence positions. Gaps may be intentionally 
    left in the final sequence range to allow for remodeling in the middle of 
    the RNA. Here, a UUUA tetraloop will be built in place of the excised 
    protein-binding loop.

9.  Example README_SETUP for de novo remodeling with a sliced input PDB:

        rna_denovo_setup.py \
            -fasta fasta \
            -secstruct_file secstruct \
            -tag native \
            -working_res 1-159 \
            -s slice_3p49_RNA.pdb \
            -cycles 20000 \
            -ignore_zero_occupancy false \

    Options:

        -fasta [fasta]               input FASTA file
        -secstruct_file [secstruct]  input secondary structure file
        -tag                         name for output files
        -working_res                 specify range of residues to model
        -s slice_3p49_RNA.pdb        see below

    The `-s` flag allows users to input a list of PDB files to use in the 
    modeling; the residues that are part of the input PDB files will not be 
    moved relative to each other, though if multiple PDB files are input, the 
    orientations of the residues in the separate files may change. In this 
    example, the full-atom refinement algorithm will be applied in the same run 
    as fragment assembly.

10. Generate command line for FARFAR modeling:

        source README_SETUP

11. Example README_FARFAR:

        rna_denovo -nstruct 500 -params_file native.params -fasta native.fasta  -out:file:silent native.out -include_neighbor_base_stacks  -minimize_rna true -s slice_3p49_RNA.pdb -input_res  1-21 26-159 -cycles 20000 -ignore_zero_occupancy false -output_res_num  1-159

12. Test command line for FARFAR modeling:

        source README_FARFAR

    This command runs a single local job on your computer. Wait until sampling 
    begins successfully (command line output similar to "Picked Fragment 
    Library for sequence u and sec. struct H ... found 2308 potential 
    fragments"), then cancel the run and submit the job to the cluster.

13. Submit jobs to the cluster:

    rosetta_submit.py README_FARFAR out [16] [1]

    The first number states how many processors to use for the run, while the 
    second number states the maximum time each job will be allowed to run 
    (walltime, in hours). Note that certain supercomputers only allow requests 
    specific multiples of processors (e.g. the Stampede cluster requires a 
    multiple of 16).

14. Concatenate all models from the out folder:

        easy_cat.py out

    Also outputs the number of models in the final silent file to the screen.

15. Extract lowest-energy models to .pdb files for viewing in PyMOL:

    extract_lowscore_decoys.py native.out [1]

    Input the number of lowest-scoring models to extract from the silent file. 
    Here, extract the single lowest-scoring model to use as the native model 
    input for comparison to the de novo models.

16. Rename lowest-score model for use as reference model

    mv native.out.1.pdb 3p49_native_RNA.pdb

17. Example pseudo-energy constraint file:

        [ atompairs ]
        O2' 2 C4' 38 FADE   0 30 15 -4.00  4.00
        O2' 2 C4' 38 FADE -99 60 30 -36.00 36.00
        O2' 1 C4' 44 FADE   0 30 15 -4.00  4.00
        O2' 1 C4' 44 FADE -99 60 30 -36.00 36.00
        O2' 5 C4' 60 FADE   0 30 15 -4.00  4.00
        O2' 5 C4' 60 FADE -99 60 30 -36.00 36.00
        O2' 2 C4' 64 FADE   0 30 15 -4.00  4.00
        O2' 2 C4' 64 FADE -99 60 30 -36.00 36.00
        O2' 25 C4' 54 FADE   0 30 15 -4.00  4.00
        O2' 25 C4' 54 FADE -99 60 30 -36.00 36.00
        O2' 45 C4' 64 FADE   0 30 15 -4.00  4.00
        O2' 45 C4' 64 FADE -99 60 30 -36.00 36.00
        O2' 45 C4' 75 FADE   0 30 15 -4.00  4.00
        O2' 45 C4' 75 FADE -99 60 30 -36.00 36.00
        O2' 32 C4' 88 FADE   0 30 15 -4.00  4.00
        O2' 32 C4' 88 FADE -99 60 30 -36.00 36.00
        O2' 42 C4' 84 FADE   0 30 15 -4.00  4.00
        O2' 42 C4' 84 FADE -99 60 30 -36.00 36.00
        O2' 48 C4' 84 FADE   0 30 15 -4.00  4.00
        O2' 48 C4' 84 FADE -99 60 30 -36.00 36.00
        O2' 55 C4' 88 FADE   0 30 15 -4.00  4.00
        O2' 55 C4' 88 FADE -99 60 30 -36.00 36.00
        O2' 55 C4' 108 FADE   0 30 15 -4.00  4.00
        O2' 55 C4' 108 FADE -99 60 30 -36.00 36.00
        O2' 58 C4' 118 FADE   0 30 15 -4.00  4.00
        O2' 58 C4' 118 FADE -99 60 30 -36.00 36.00
        O2' 67 C4' 119 FADE   0 30 15 -4.00  4.00
        O2' 67 C4' 119 FADE -99 60 30 -36.00 36.00
        O2' 67 C4' 121 FADE   0 30 15 -4.00  4.00
        O2' 67 C4' 121 FADE -99 60 30 -36.00 36.00
        O2' 78 C4' 113 FADE   0 30 15 -4.00  4.00
        O2' 78 C4' 113 FADE -99 60 30 -36.00 36.00
        O2' 78 C4' 135 FADE   0 30 15 -4.00  4.00
        O2' 78 C4' 135 FADE -99 60 30 -36.00 36.00
        O2' 42 C4' 157 FADE   0 30 15 -4.00  4.00
        O2' 42 C4' 157 FADE -99 60 30 -36.00 36.00
        O2' 74 C4' 156 FADE   0 30 15 -4.00  4.00
        O2' 74 C4' 156 FADE -99 60 30 -36.00 36.00
        O2' 100 C4' 148 FADE   0 30 15 -4.00  4.00
        O2' 100 C4' 148 FADE -99 60 30 -36.00 36.00
        O2' 100 C4' 145 FADE   0 30 15 -4.00  4.00
        O2' 100 C4' 145 FADE -99 60 30 -36.00 36.00
        O2' 113 C4' 153 FADE   0 30 15 -4.00  4.00
        O2' 113 C4' 153 FADE -99 60 30 -36.00 36.00
        O2' 135 C4' 154 FADE   0 30 15 -4.00  4.00
        O2' 135 C4' 154 FADE -99 60 30 -36.00 36.00
        O2' 5 C4' 119 FADE   0 30 15 -4.00  4.00
        O2' 5 C4' 119 FADE -99 60 30 -36.00 36.00
        O2' 25 C4' 88 FADE   0 30 15 -0.80  0.80
        O2' 25 C4' 88 FADE -99 60 30 -7.20  7.20
        O2' 37 C4' 62 FADE   0 30 15 -0.80  0.80
        O2' 37 C4' 62 FADE -99 60 30 -7.20  7.20
        O2' 79 C4' 103 FADE   0 30 15 -0.80  0.80
        O2' 79 C4' 103 FADE -99 60 30 -7.20  7.20
        O2' 15 C4' 88 FADE   0 30 15 -0.80  0.80
        O2' 15 C4' 88 FADE -99 60 30 -7.20  7.20
        O2' 32 C4' 108 FADE   0 30 15 -0.80  0.80
        O2' 32 C4' 108 FADE -99 60 30 -7.20  7.20
        O2' 9 C4' 138 FADE   0 30 15 -0.80  0.80
        O2' 9 C4' 138 FADE -99 60 30 -7.20  7.20
        O2' 25 C4' 118 FADE   0 30 15 -0.80  0.80
        O2' 25 C4' 118 FADE -99 60 30 -7.20  7.20

18. Example README_SETUP:

        rna_denovo_setup.py \
            -fasta fasta \
            -secstruct_file secstruct \
            -fixed_stems \
            -no_minimize \
            -tag glycine_riboswitch \
            -working_res 1-159 \
            -native 3p49_native_RNA.pdb \
            -cst_file constraints \
            -staged_constraints \
            -cycles 20000 \
            -ignore_zero_occupancy false \
            -silent helix0.out helix1.out helix2.out helix3.out helix4.out helix5.out helix6.out helix7.out \
            -input_silent_res 2-9 65-72 16-21 26-31 33-35 54-56 39-42 48-51 81-85 155-159 92-97 101-106 108-110 145-147 114-117 139-142 \

    Options:

        -fasta [fasta]	input FASTA file
        -secstruct_file [secstruct]	input secondary structure file
        -fixed_stems	specify whether helices should be fixed
        -no_minimize	specify not to perform full-atom refinement; minimization will be performed in the next stage of modeling
        -tag	name for output files
        -working_res	specify range of residues to model
        -native [native.pdb]	input reference or native model; used for benchmarking cases and will return rms calculations for all models (see command line [5])
        -cst_file [constraints]	input file with pseudoenergy constraints
        -staged_constraints	apply constraints
        -ignore_zero_occupancy false	
        -silent [helix0.out helix1.out …]	input silent files with pre-assembled helices
        -input_silent_res [2-9 65-72 16-21 26-31 …]		specify position ranges of helices in silent files

19. Generate command line for FARFAR modeling:

        source README_SETUP

20. Example README_FARFAR:

        rna_denovo -nstruct 500 -params_file glycine_riboswitch.params -fasta glycine_riboswitch.fasta  -out:file:silent glycine_riboswitch.out -include_neighbor_base_stacks  -minimize_rna false -native glycine_riboswitch_3p49_native_RNA.pdb  -in:file:silent helix0.out helix1.out helix2.out helix3.out helix4.out helix5.out helix6.out helix7.out -input_res  2-9 65-72 16-21 26-31 33-35 54-56 39-42 48-51 81-85 155-159 92-97 101-106 108-110 145-147 114-117 139-142 -cst_file glycine_riboswitch_constraints -staged_constraints -cycles 20000 -ignore_zero_occupancy false -output_res_num  1-159

21. Test command line for FARFAR modeling:

        source README_FARFAR

22. Submit jobs to the cluster:

        rosetta_submit.py README_FARFAR out [96] [16]

23. Concatenate all models from the out folder:

        easy_cat.py out

24. Extract lowest-energy models to .pdb files for viewing in PyMOL:

        extract_lowscore_decoys.py glycine_riboswitch.out [15]

25. Example MINIMIZE:

        parallel_min_setup.py -silent glycine_riboswitch.out -tag glycine_riboswitch_min -proc [96] -nstruct [2000] -out_folder min_out -out_script min_cmdline "-native glycine_riboswitch_3p49_native_RNA.pdb -cst_fa_file glycine_riboswitch_constraints -params_file glycine_riboswitch.params -ignore_zero_occupancy false -skip_coord_constraints"

    The first number states how many processors to use for the run, while the 
    second number is 1/6 the total number of previously generated FARNA models. 
    If you are running on a supercomputer that only allows specific multiples 
    of processors, use an appropriate number for the first input.

26. Generate command lines for full-atom refinement:

        source MINIMIZE

27. Example command line from min_cmdline to run as test:

        rna_minimize -native glycine_riboswitch_3p49_native_RNA.pdb -cst_fa_file glycine_riboswitch_constraints -params_file glycine_riboswitch.params -ignore_zero_occupancy false -skip_coord_constraints -in:file:silent min_out/0/0.silent -out:file:silent min_out/0/glycine_riboswitch_min.out

28. Submit jobs to the cluster:

        rosetta_submit.py min_cmdline min_out [1] [16]

    The first number states how many processors to use for each line in 
    min_cmdline. Here, enter 1 for the first input so that the total number of 
    processors used will be equal to the number of processors entered with the 
    `-proc` flag in command line [12], above. The second number states the 
    maximum time each job will be allowed to run (walltime).

29. Concatenate all models from the min_out folder:

        easy_cat.py min_out

30. Sort models by Rosetta energy and select a subset for clustering:

        silent_file_sort_and_select.py [glycine_riboswitch_min.out] -select [1-60] -o [glycine_riboswitch_min_sort.out]

    The range of models under the -select tag includes 0.5% of the total number 
    of FARNA models generated previously. Outputs a new silent file containing 
    selected number of lowest-energy models.

31. Cluster models:

        cluster -in:file:silent glycine_riboswitch_min_sort.out -in:file:fullatom -out:file:silent_struct_type binary -export_only_low false -out:file:silent cluster.out -cluster:radius [radius]

    Select a radius so that 1/6 of the models in the input sorted silent file 
    are in the largest cluster (cluster0) of models.

32. Copy clustered .out file to a new file to isolate cluster0:

        cp cluster.out cluster0.out

33. Extract lowest-energy models to .pdb files for viewing in PyMOL:

        extract_lowscore_decoys.py cluster0.out [15] –no_replace_names

    Input the number of models in cluster0. The -no_replace_names tag preserves 
    the filenames of the cluster members to reflect their order in the cluster, 
    rather than renaming them in order of Rosetta energy score.

34. Cut out a segment of a model:

        pdbslice.py [3p49_native_RNA.pdb] -subset [2-72 81-159] [slice_kinkturn_]

    Here, the 3P49 crystal structure includes an additional G at position 0, 
    which must be excised to allow the leader sequence to be added to the 5´ 
    end, and the internal linker that forms the kink-turn motif with the leader 
    sequence is also excised to allow remodeling.

35. Renumber a PDB:

        renumber_pdb_in_place.py [slice_kinkturn_3P49_native_RNA.pdb] [10-80 89-167]

    Here, the PDB is renumbered to allow the leader sequence to be added at the 
    5´ end.

36. Example revised FASTA file:

        >3P49_RNA_kinkturn.pdb
        ucggaugaagauaugaggagagauuucauuuuaaugaaacaccgaagaaguaaaucuuucagguaaaaaggacucauauuggacgaaccucuggagagcuuaucuaagagauaacaccgaaggagcaaagcuaauuuuagccuaaacucucagguaaaaggacggag

37. Example revised secondary structure file:

        (((......((((((((......((((((....)))))).(((...((((.....))))..)))........))))))))...)))..(((((......((((((...)))))).(((...((((....((((....)))).....))))..))).......)))))

38. Example README_SETUP for de novo remodeling with a sliced input PDB:

        rna_denovo_setup.py -fasta fasta2 -secstruct_file secstruct2 \
         -fixed_stems \
         -tag glycine_rbsw_kinkturn \
         -working_res 1-167 \
         -s slice_kinkturn_3P49_native_RNA.pdb \
         -cycles 20000 \
         -ignore_zero_occupancy false \

39. Thread an RNA sequence into a template structure:

        rna_thread –in:file:fasta [fasta] -in:file:s [template PDB] –o [output PDB]

    The first input is a FASTA file containing two RNA sequences: 1. the 
    sequence of interest, onto which the structure of the template sequence 
    will be threaded, and 2. the template sequence. The template sequence 
    should be truncated to the regions into which the sequence of interest will 
    be threaded; use hyphens (‘-’) to align the template sequence with the 
    target sequence in the FASTA file. The second input, the template structure 
    in PDB format, should be similarly truncated, using pdbslice.py if 
    necessary. If the template PDB is not correctly formatted for Rosetta 
    modeling, use make_rna_rosetta_ready.py to reformat it. The last input is 
    the name of the output PDB.

Further documentation for RNA threading in Rosetta can be found on the 
RosettaCommons website:  
https://www.rosettacommons.org/docs/latest/rna-thread.html


#Relax with all-heavy-atom constraints

##Introduction

We looked for a way to simultaneously minimize rosetta energy and keep all heavy atoms in a crystal structure as close as possible to their starting positions. As hard-won experience has shown, simply running relax on a structure will often move the backbone a few Angstroms. The best way we have found to perform the simultaneous optimization is to run relax with constraints always turned on (typically constraints ramp down in the late cycles of a relax run) and to constrain not just backbone but also sidechain atoms. This protocol has been tested on a benchmark set of 51 proteins and found to increase sequence recovery in enzyme design by 5% as compared with design in raw pdb structures. It accomplishes this with only .077 Angstrom RMSD over the set of proteins (C-alpha RMSD) from raw pdb to relaxed-with-csts pdb. A more complete description of the data leading to this protocol is below.


## Protocol

### Prepare structures for relax

The required files are in: `rosetta/rosetta_source/src/apps/public/relax_w_allatom_cst`
Many pdbs have features, such as non-canonical amino acids, which will cause rosetta to fail. This script will be able to process most input pdbs so they are ready for a relax run. The script will "clean" structures to replace non-canonical amino acids with their closest counterparts: 
```
rosetta/rosetta_source/src/apps/public/relax_w_allatom_cst/clean_pdb_keep_ligand.py

python rosetta/rosetta_source/src/apps/public/relax_w_allatom_cst/clean_pdb_keep_ligand.py your_structure_original.pdb -ignorechain
```

### Relax with all-heavy-atom constraints: Short protocol (recommended)

Relax with all-heavy-atom constraints is built into the relax application itself. If this is a new structure you may want to first clean it up using the above script. Relax proceeds as follows:
(Note that these flags will not preserve any ligands in the structure unless the ligand params file is added.)
```
mkdir test; cp starting_inputs/1A99_1A99.pdb test/; cd test/

../../../../rosetta_source/bin/relax.default.linuxgccrelease -database ../../../../rosetta_database/ -s 1A99_1A99.pdb @../starting_inputs/flags2 > log2.txt
```
These flags are required: 
```
-relax:constrain_relax_to_start_coords
-relax:coord_constrain_sidechains
-relax:ramp_constraints false
```
The flags2 file includes a set of recommended flags:
```
-ignore_unrecognized_res
-relax:constrain_relax_to_start_coords
-relax:coord_constrain_sidechains
-relax:ramp_constraints false
-ex1
-ex2
-use_input_sc
-correct
-no_his_his_pairE
-no_optH false
-flip_HNQ
```

### Relax with all-heavy-atom constraints: Longer Protocol (not recommended)

In general the short protocol is preferred for most applications, since this version is more complicated and the two give nearly identical results. In this protocol a separate script first generates sidechain atom constraints from an input pdb, then the relax protocol is run with this pre-generated constraitn file. The shorter protocol does this all in one step, and this longer version is largely deprecated. Certain users might prefer this protocol because it allows you to see a list of all constraints, and perhaps to modify constraints using other scripts/data, prior to relax. 

(a) Generate sidechain coordinate constraints on your pdb using sidechain_cst_3.py:
```
mkdir test2; cp starting_inputs/1A99_1A99.pdb test2/; cd test2/
python ../../../../rosetta_source/src/apps/public/relax_w_allatom_cst/sidechain_cst_3.py 1A99_1A99.pdb 0.1 0.5
[output:1A99_1A99_sc.cst]
```

(b) Run relax using these constraints and with a custom relax script to force constraints to stay on during the entire run. 
Note that these flags will not preserve any ligands in the structure unless the ligand params file is added. If you want to keep the ligand, simply copy it over from the input pdb. 
```
cd test2/
../../../../rosetta_source/bin/relax.default.linuxgccrelease -database ../../../../rosetta_database/ -s 1A99_1A99.pdb -constraints:cst_fa_file 1A99_1A99_sc.cst @../starting_inputs/flags > log.txt
```
flags file:
```
-constrain_relax_to_start_coords
-relax:script ../../../../rosetta_source/src/apps/public/relax_w_allatom_cst/always_constrained_relax_script
-ignore_unrecognized_res
-preserve_header
-ex1
-ex2
-use_input_sc
-correct
-no_his_his_pairE
-score::hbond_params correct_params
-lj_hbond_hdis 1.75
-lj_hbond_OH_donor_dis 2.6
-linmem_ig 10
-nblist_autoupdate true
-dun08 false
-no_optH false
-flip_HNQ
-chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation SpecialRotamer
 VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals
pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated
thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated
lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated
cys_acetylated tyr_diiodinated N_acetylated C_methylamidated
MethylatedProteinCterm
```
note: Including extra rotamers is important if your goal is to keep all sidechain atoms tightly constrained; if that is not important for your applications, by all means exclude the ex1 and ex2.

### Relax with all-heavy-atom constraints: Data
To test protocols for relaxation of input pdbs we used the 51 scaffold test set for enzdes. We run design over the input structures and calculate the percent of residues which come back with the native identity -- sequence recovery, only over the designed residues. This is of course an imperfect metric (the original sequence might not be fully optimal), but it allows us to ask how many residues rosetta will correctly choose, assuming that the input structure is already at a minimum in sequence space for the ligand in question. It had already been found that relax alone (with no constraints) will distort most structures, and that those structures will give a much higher sequence recovery in design, but this is a result of the distortion to the input structure.

We decided to test a few protocols to find the one that would best minimize RMSD from the original pdb while maximizing sequence recovery. All calculations are averages over 50 runs of the 51 scaffold set. All RMSDs are to native C-alpha. It is also possible to read in electron density for a pdb and use that as a constraint during relax, but electron density is not uniformly available, and we were looking for a protocol that would work in every case even if the electron density was lacking. We chose the first protocol as the best because it minimized rmsd (and sidechain motions) while maximizing sequence recovery.

Running relax with sidechain coordinate constraints and bb coord constraints: (ex flags and use native; ex flags on in enzdes)

    0.447 sequence recovery (0.077 RMSD) [-557 totalscore]

Running relax with sc-sc distance constraints at 3.5 distance cutoff. Note that this gets similar rmsd minimization but doesn't maximize sequence recovery or minimize score as well as coordinate constraints. We tested this protocol with a variety of sc-sc distance constraint cutoff values and found it to slightly but systematically underperform coordinate constraints:

    0.436 sequence recovery (0.0706 RMSD) (-534 totalscore)

Running elax with backbone constraints only. Note that this does worse in terms of rmsd and that many sidechains are a few angstroms off:

    0.488 sequences recovery (0.098 RMSD) (-633 totalscore)

No relax benchmark and rosetta scoring of native input structures:

    0.40 (0 RMSD by definition) (-194.7 totalscore)

At this point, the astute reader might ask, what score terms became 438 Rosetta Energy Units (REU) better? We ranked the difference in scores over all structures, comparing the all-atom coordinate constraint protocol and the non-relaxed input structure. The biggest difference is fa_dun (-192.4), followed by fa_rep (-108.1), pro_close (-26.4), hbond_sc (-10.8) and omega (-8.5). Many input rotamers are close to but not in a "good" Dunbrack rotamer, and the backbone has to be slightly tweaked in order for that residue to get a good dunbrack score. Also many atoms are slightly too close, and they give the fa_rep contribution.


# RNA threading/mutation

# Author
Rhiju Das, rhiju@stanford.edu

## Brief Description

Take an RNA template for homology modeling, and thread in a new sequence.

## Abstract

RNA homology modeling often involves taking a template coordinate file and threading on a new sequence. In the midst of the 2011 "RNA puzzles" community-wide blind trials, I hacked together some Rosetta code to do this.


# Running

```
rna_thread -in:file:fasta rosetta_inputs/3bo3_REARRANGE_to_GIR1.fasta -s rosetta_inputs/3bo3_REARRANGE.pdb  -o 3bo3_REARRANGE_to_GIR1_thread.pdb -seq_offset 63
```

Note that insertions will *not* be modeled -- there will just be a chainbreak at those residues, and the residue numbering will skip.

If you don't have any alignment gaps in the target model, you can skip the creation of the alignment.fasta file, and just supply -seq <target sequence> on command line. See demo for rna_mutate.


Floppy Tail Demo
================

This demo contains the starting structure for the original published use of 
FloppyTail, and is a copy of FloppyTail's integration test.  To run the demo 
(in the inputs subfolder):

    FloppyTail.linuxgccrelease @options

Read the options file first to set local paths as needed.

Please refer to FloppyTail's documentation (in the [online manual](https://www.rosettacommons.org/docs/latest/floppy-tail.html), or in your copy of the code at 
Rosetta/documentation/application_documentation/structure_prediction/floppy-tail.md, and its [publication](http://www.ncbi.nlm.nih.gov/pubmed/19945379) 
Kleiger G, Saha A, Lewis S, Kuhlman B, Deshaies RJ. Rapid E2-E3 assembly and 
disassembly enable processive ubiquitylation of cullin-RING ubiquitin ligase 
substrates. Cell. 2009 Nov 25;139(5):957-68. PubMed PMID: 19945379.
Note that the paper is primarily about biology, not modeling.

Note that the provided outputs are from the integration test, which runs in ~30 
s.  You will need to run FloppyTail for much longer (both longer individual 
runs, and many trajectories) to get scientifically useful results.  The options 
file and documentation provide details on how to do that.
the two scripts should be used in tandem. See RosettaScripts documentation for details on how to use this. The scripts are hard-coded for a particular example and should be changed to fit different examples.
this is a new (and not very well tested) strategy for doing multicriterion Monte Carlo. The idea is to optimize more than one criterion, e.g., the total score and the binding energy. You select positions in which to introduce random mutations and then you accept these mutations subject to MC on both criteria. Simple and effective. It's a good idea to be as restrictive with the designable positions (and their identities) as possible since sampling can never be complete.
# Zinc Site Redesign

This is a demo for the mononuclear zinc metalloenzyme redesign procedure described in Khare, Kipnis, Greisen et al. Nature Chemical Biology (2012).

Authors: Sagar Khare (khares@uw.edu), Per Jr Greisen (pgreisen@gmail.com)

## Summary

Starting from a TS ensemble model and a set of zinc-containing PDBs (only one - PDBid 1A4L - is included here, list of other PDB codes used in the paper is in list\_of\_input\_pdbs file) as inputs, we generate a design model and evaluate it. 

## Input

Requires an input TS model including the metal, and a (set of) PDB file(s). The TS ensemble is generated using a starting molfile, and converted to Rosetta parameters using the script `~/rosetta/rosetta_source/src/python/apps/public/molfile_to_params.py`

## Overview of the Procedure

Each Step is included as a directory with a separate README.

1.	Analysis of ZincSite
2. 	CleanPDB after running the zinc analysis. 
3.	Align Transition State model(s) to the Zinc
4.	Minimization in polyAla pocket (optional)
5.	Matching to introduce additional catalytic residues
6.	Design to maximize TS affinity
7.	Revert to Native - Semi-automated refinement
8. 	Dock (validation)

## Details

For most Python scripts, running with a -h option will give a list of available options.

1. Given a PDB file, obtain the co-ordination sphere details of metal (how many protein-metal, and HETATM-metal interactions, in what geometry in each chain).

    `python analyze_zinc_site.py -f 1A4L.pdb`

2. Given an input PDB, "clean" it i.e. keep one chain, change MSE->MET, KCX->LYS etc. Keep HETATMS.

    `python cleanPDBfile.py -f 1A4L.pdb`

3. Align transition state model(s) onto the zinc site according to co-ordination sphere details: 

	1. Given TS model (LG_0001.pdb) and "clean" PDBfile (1A4L_clean.pdb) align them.

	    `python align.py -f 1A4L_clean_A.pdb -l LG_0001.pdb`

	2. Generate Rosetta Inputs given the aligned ligand.

	    `python generate_metal_cstfile.py -f 1A4L_clean_A.pdb -m ZN -a aligned_ligand.pdb`

4. Minimize the constraints energy for the superimposed ligand using Rosetta. Except the metal-chelating protein residues, every other ligand-proximal protein residue is temporarily converted to Ala.

    `~/rosetta/bin/enzyme_design.static.linuxiccrelease @optcst.flags -linmem_ig 10 -in:file::s rosetta_cst.pdb`

5. To introduce additional catalytic residues, run a round of RosettaMatch.

    `~/rosetta/bin/match.linuxiccrelease @general_matching.flags @scaf.flags @subs.flags  -linmem_ig 10 -in:file::s 1A4LA_clean_r.pdb`

6. Design the rest of the pocket for maximizing TS affinity.

    `~/rosetta/bin/enzyme_design.linuxiccrelease @enzdes.flags -correct -linmem_ig 10 -in:file::s <file_from_matching.pdb> > design.log&`

7. For a semi-automated refinement step, generate a so-called resfile, so that the wildtype and any other user-specified residue can be introduced at any given position. Use this resfile with similar commandline as Step 6.

8. To validate the design, perform docking of the ligand in the designed pocket to check if there are any alternative binding modes with similar energies.
Next step is to optimize the interaction between TS model and protein
- design process - here new residues will be inserted in the protein
by Rosetta. 

Inputs:

1. A file from the matching step
2. Constraint file with the new constraint added
3. Rotamer file for the TS model
4. Parameter file for the TS model
5. Flag file for the design enzdes.flags

Usage:

./run.sh UM_1_H15H17H214D295Q58_1A4L_clean_A_r_1A4L_clean_A_1.pdb

Output:

UM_1_H15H17H214D295Q58_1A4L_clean_A_r_1A4L_clean_A_1__DE_1.pdb
enz_score.out

Here we search for additional catalytic residues to stabilize the TS using RosettaMatch. 
In this example, we are looking for a Q to make Hbond to the attacking nucleophile hydroxyl.
Match-style constraint block corresponding to this interaction is appended to the constraints file obtained in the previous step. 

The following files are needed

1. A protein scaffold file (1A4LA_clean_r.pdb) which includes a description of the pre-existing metal-protein interactions in the REMARK field. The protein co-ordinates could be the minimized ones from the previous CstOpt step or the starting scaffold co-ordinates without minimization (as in this example).
2. Flags for the matching: general_matching.flags subs.flags scaf.flags
3. all.pos: file describing what positions in the scaffold are to be used for which constraint (default all positions for each constraint)         
4. constraints.cst: geometric constraints including the pre-existing metal site constraints, and the "new" desired interaction.
5. Parameters for the ligand (LG.params) and optionally ligand rotamers as a concatanated PDB file - in this example, all atoms except the phosphorous and oxygens bonded to it are made virtual in the LG.params file to avoid using the rotamer ensemble for speed reasons.

Usage:
./match.sh 1A4LA_clean_r.pdb

Output:
If any matches are found new pdb files will be return with the residue inserted, named  UM_*pdb

Notes: For more details on the setup see the matcher documentation
The script will analyze the zinc site in the PDB file and print type
of metal site. 

Input parameter are PDB file and type of metal ion - it has mainly
been tested with Zinc ions. 

Example with 1A4L. 

python analyze_zinc_site.py -f 1A4L.pdb

Output from analysis

Number of chains in pdb file:  4
Number of interactions between metal and protein: 4
Number of interactions between hetero atom : 1
PDB chain D contains 5 coordinated metal site
Number of interactions between metal and protein: 4
Number of interactions between hetero atom : 1
PDB chain C contains 5 coordinated metal site
Number of interactions between metal and protein: 4
Number of interactions between hetero atom : 1
PDB chain A contains 5 coordinated metal site
Number of interactions between metal and protein: 4
Number of interactions between hetero atom : 1
PDB chain B contains 5 coordinated metal site

To summarize: 4 chains with a triangular bipyramidal zinc site.
Step 3.1:

Align the transition state model on to the zinc site in the "clean" protein PDB file. 

Input:
Protein PDB file: 1A4L_clean_A.pdb
Ligand PDB file: LG_0001.pdb

Usage: 

python align.py -f 1A4L_clean_A.pdb -l LG_0001.pdb

Output:
A file named aligned_ligand.pdb which is the transformed co-ordinates of the ligand after alignment (protein is held fixed).

Notes:

MANUALLY SET THE ATOMNAMES USED TO ALIGN CORRECTLY IN THE ALIGN.PY FILE 
ORDER DOES MATTER: ZINC SHOULD BE SPECIFIED FIRST.

------------------------------------------------------------------------------------------------------------------------------------------
Step 3.2

Generate files for next steps in Rosetta

Input:

Clean PDB file: 1A4L_clean_A.pdb
Aligned ligand: aligned_ligand.pdb
name of metal: ZN

Usage:

python generate_metal_cstfile.py -f 1A4L_clean_A.pdb -m ZN -a aligned_ligand.pdb 

Output:

Protein-ligand match contraints in a file (constraint.cst)
PDB for Rosetta: rosetta_cst.pdb

Other options:

-n to custom name the PDB for Rosetta 
Here other chains are removed  and residues like MSE are changed to
MET, remove structural metal sites, and change KCX to LYS. The default
chain that is kept is chain A. The script is executed by 

python cleanPDBfile.py -f 1A4L.pdb 

To keep chain B for example:

python cleanPDBfile.py -f 1A4L.pdb -c B

Other options: To remove a specific metal ion from PDB

-m: name of metal to remove from pdb
-n: residue number of metal to remove from pdb

python cleanPDBfile.py -f PDBFILE -m METAL -n NR

This is an optional step in the protocol where we optimize the protein-metal interaction in a polyAla context (i.e. everything in the pocket except the TS model and metal-chelating residues are trimmed back to Ala, and the constraints energy is optimized by minimization. See enzdes documentation for more details.)

Input files:

1. rosetta_cst.pdb: contains remark lines to specify catalytic residues, co-ordinates of the TS model
2. constraint.cst: contains match-style constraints between metal site and TS model
3. con_rotamers.pdb.gz : contains rotamers for the TS model
4. LG.params: parameter file for the TS model, which includes charges and topology
   of TS model. Here the following lines is added to the end of the file
   PDB_ROTAMERS con_rotamers.pdb.gz
   which ensures that the rotamers of the ligand are used for the optimization
5. optcst.flags
   Parameters to rosetta for the optimization

Usage:

./min.sh 

Output files:
PDB file for minimized interface:
e.g. rosetta_cst__DE_1.pdb

Typically this minimization is repeated multiple times to get slightly different conformations (number is controlled by the -nstruct option)

The resfile contains reversions back to the native sequence of the
protein. The script generates a template residue file where one can
change residue back to native or test certain residues at different
position
	python generate_residuefile.py rosetta_cst.pdb
returns a resfile with type NATAA, NATRO (keep rotamer) - one should
add PIKAA for picking new residue e.g.
    1 PIKAA WFY 
As one of the test to evaluate the design we dock the TS structure
into the protein. 

The docking is described in detail in Davis & Baker JMB 2009.

Here a new LG.params is needed as the zinc ion is
kept in the protein and not free to move. 

A constraint file is used for the zinc-ligand distance.

Number of alternative conformations is controlled by -nstruct. 
# StepWise Monte Carlo (examples for RNA)

## Author
Rhiju Das, rhiju@stanford.edu

## Brief Description

Solve structure of an RNA loop or motif in the context of a starting structure.

## Abstract

Ab initio and comparative modeling of biopolymers (RNA, protein, protein/RNA) often involves solving well-defined small puzzles (4 to 20 residues), like RNA aptamers, RNA tertiary contacts, and RNA/protein interactions. If these problems have torsional combinations that have not been seen previously or are not captured by coarse-grained potentials, most Rosetta approaches will fail to recover their structures.  This app implements a stepwise ansatz, originally developed as a 'stepwise assembly' enumeration that was not reliant on fragments or coarse-grained modeling stages, but was computationally expensive. The new mode is a stepwise monte carlo, a stochastic version of stepwise assembly. 


## Running

### Example Rosetta Command Line
```
stepwise -in:file:fasta rosetta_inputs/1zih.fasta -s rosetta_inputs/start_helix.pdb  -out:file:silent swm_rebuild.out -extra_min_res 4 9
```
Currently, we are mainly using a scorefunction with a more stringent torsional and repulsive potential, enabled by flags `-score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj`. These may become default soon (in 2015 onwards), after further tests.

To get out models:

```
extract_pdbs -silent swm_rebuild.out 
```

(Or use extract_lowscore_decoys.py which can be installed via tools/rna_tools/.)

### Example Rosetta Command Line for Design
Simply use a fasta file that has n's at positions you want to design.

```
stepwise -in:file:fasta rosetta_inputs/NNNN.fasta -s rosetta_inputs/start_helix.pdb  -out:file:silent swm_design.out -extra_min_res 4 9
```

# Application of the Stepwise Assembly method to RNA loop modeling (SWA_RNA_LOOP)


## Authors
- Parin Sripakdeevong (sripakpa [at] stanford.edu)
- Rhiju Das (rhiju [at] stanford.edu)

Last updated on March 18, 2012.
 
# Purpose and Algorithm

This demo illustrates a protocol to built single-stranded RNA loops using a deterministic, enumerative sampling method called Stepwise Assembly. The modeling situation considered here is the lock-and-key problem. Given a template PDB that contains nucleotides surrounding a missing RNA loop, the Stepwise Assembly method finds the loop conformation (the key) that best fits the surrounding structure (the lock). Details of this method is described in "An enumerative stepwise ansatz enables atomic-accuracy RNA loop modeling" by P. Sripakdeevong, W. Kladwang, and R. Das (2012), Proc Natl Acad Sci USA.
	As detailed in the paper, the Stepwise Assembly method constructs full-length RNA loops through the recursive building of each individual RNA nucleotides over multiple steps. The enumerative nature of the method makes the full-calculation quite computationally expensive, requiring for example 15,000 CPU hours to build a single 6-nucleotides RNA loop. While this full-calculation is now feasible on a high-performance computer clusters, performing the full-calculation in this demo would be too excessive. 
	Instead we will illustrate in this demo, the Stepwise Assembly protocol to build the first (5' most) nucleotide of a 6-nucleotides RNA loop. Performing the individual building step takes roughly 15 minutes on an Intel Core i7 2.66 GHz processor. The same Stepwise Assembly protocol can then be recursively applied to build the remaining 5 nucleotides of the loop (see SI of the referenced paper for details).

# Required Tools and Input files

There are two required files: 
The template_PDB file (template.pdb): A PDB file containing the coordinates of surrounding nucleotides in the vicinity of the missing RNA loop to be build. We recommend including all surrounding nucleotides within a 10-Angstrom vicinity of the missing RNA loop. Supplied PDB file must be in the Rosetta RNA PDB format.

The fasta file (fasta): this is the sequence file of the full-length RNA. The fasta file has the RNA name on the first line (after >), and the sequence on the second line. Valid letters are a, c, g and u. 

# Optional Tools and Input Files

## Optional additional files:
The native_PDB file (native.pdb): A PDB file containing the 'native' crystallographic or NMR structure. The PDB file should contain the coordinates of the native loop nucleotides plus the coordinates of the surrounding nucleotides inherited from template_PDB. The supplied native_PDB file is not used to guide the modeling process and only used for reporting the RMSD of the generated rosetta models to the native loop. Supplied PDB file must be in the Rosetta RNA PDB format.

# How to run the job

The SWA_RNA_python package located at `rosetta_tools/SWA_RNA_python/` contains the scripts necessary to setup and run the Stepwise Assembly protocol. Instructions are provided in steps 1)-4) below: 

1. Specify the location of the rosetta bin folder and rosetta database folder by editing the file `rosetta_tools/SWA_RNA_python/SWA_dagman_python/utility/USER_PATHS.py`

    For example, if the main rosetta folder is located at `~/rosetta/`, then the file should look as follow:
    ```
	    #!/usr/bin/python

	    USER_ROSETTA_BIN_FOLDER="~/rosetta/rosetta_source/bin/"
	    USER_ROSETTA_DATABASE_FOLDER="~/rosetta/rosetta_database/"
    ```

2. Add the SWA_RNA_python package location to the PYTHON path. For bash shell users, the location can be directly added to the `~/.bashrc` file:

    ```
	export PYTHONPATH=$PYTHONPATH:~/rosetta/rosetta_tools/SWA_RNA_python/
    ```

3. After the paths are correctly specified, the following command is used to setup everything needed run the Stepwise Assembly job:

    ```
	rosetta_tools/SWA_RNA_python/SWA_dagman_python/SWA_DAG/setup_SWA_RNA_dag_job_files.py -s template.pdb -fasta fasta -sample_res 3-8 -single_stranded_loop_mode True -local_demo True -native_pdb native.pdb
    ```

    The "-s" flag specifies the template_PDB file

    The "-fasta" flag specifies the fasta file

    The "-sample_res" flag specifies the sequence number of nucleotides in the missing loop. In this demo case, this correspond to nucleotides at sequence number 3 4 5 6 7 and 8.

    The "-single_stranded_loop_mode" flag specifies that the job involve modeling a single-stranded loop (i.e. the lock-and-key problem).

    The "-local_demo" flag indicate that this is demo to be run on a local laptop or desktop. The calculation perform here is to only build the first (5' most) nucleotide of the 6-nucleotides RNA loop.

    The "-native_pdb" flag specifies the native_PDB file and is optional. 

4. Type `source LOCAL_DEMO` to execute the Rosetta protocol.

The provided instruction will allow the user to build the first (5' most) nucleotide of a N-nucleotide loop. As previously stated, the full-calculation to build the full-length RNA loops is quite computationally expensive and is beyond the scope of this demo. The SWA_RNA_python package is, however, equipped to run this recursive full-calculation on a high-performance computer clusters. The package utilize concept familiar from the Map/Reduce Direct Acyclic Graph framework to order the calculation steps and allocate resources to recursive build the full-length RNA loop over multiple steps, one individual RNA nucleotide at a time. If any user is interested, please contact Parin Sripakdeevong (sripakpa [at] stanford.edu) and we will be happy to provide additional instructions.

# Expected Outputs

The expected outputs are two silent_files:
A) region_0_1_sample.out: This silent_file contain 108 structures, corresponding to the 108 lowest energy conformations.
B) region_0_1_sample.cluster.out: Same as A) but after clustering of the models to remove redundant conformations.

In both silent_files, the total energy score is found under the 'score' column. If the "native_pdb" flag was included, then the RMSD (in angstrom units) between the native_pdb and each Rosetta model is found under the 'NAT_rmsd' column.

Finally, use the following command to extract the top 5 energy cluster centers:

```
	rosetta_tools/SWA_RNA_python/SWA_dagman_python/misc/SWA_extract_pdb.py -tag S_0 S_1 S_2 S_3 S_4  -silent_file region_0_1_sample.cluster.out
```

After running the command, the extracted PDB files should appear in the pose_region_0_1_sample.cluster.out/ subfolder.

# StepWise Monte Carlo (examples for RNA)

## Author
Rhiju Das, rhiju@stanford.edu

## Brief Description

Solve structure of a mini-protein

## Abstract

Ab initio and comparative modeling of biopolymers (RNA, protein, protein/RNA) often involves solving well-defined small puzzles (4 to 20 residues), like RNA aptamers, RNA tertiary contacts, and RNA/protein interactions. If these problems have torsional combinations that have not been seen previously or are not captured by coarse-grained potentials, most Rosetta approaches will fail to recover their structures.  This app implements a stepwise ansatz, originally developed as a 'stepwise assembly' enumeration that was not reliant on fragments or coarse-grained modeling stages, but was computationally expensive. The new mode is a stepwise monte carlo, a stochastic version of stepwise assembly. 


## Running
### Example Rosetta Command Line

```
stepwise -fasta rosetta_inputs/2jof.fasta -native rosetta_inputs/2jof.pdb -score:weights stepwise/protein/protein_res_level_energy.wts -silent swm_rebuild.out -cycles 2000 -nstruct 50
```

Most of the simulation may be spent flickering bits of secondary structure -- in the future, we will probably setup some precomputation of these bits so that computation can be focused on build up of the complete mini-protein structure.



Clustering Demo
===============

The following command is a good place to start:

    /path/to/rosetta/bin/cluster.linuxgccrelease -database /path/to/database -in:file:s rosetta_inputs/*.pdb -in:file:fullatom -cluster:gdtmm -cluster:radius -1 -cluster:population_weight 0.0 -sort_groups_by_energy 

These flags may be also be useful:

    -in:file:silent <filename.out>
    -in:file:silent_struct_type <desired type but binary is best choice>
    -limit_cluster_size <n>
    -limit_clusters <n>
    -cluster:radius <for gdt you want some number between 10-50>
    -out:file:silent clustered.out -out:file:silent_struct_type binary

You can extract the silent file by running the mini app `extract_pdbs.linuxgccrelease`.

If you have thousands of structures it may be best to first do an energy cut to filter out high energy structures.
The pdb's should be first scored and stored in an out file.
Try using the following app

    score_jd2 -in:file:s *.pdb -out:file:silent scored.out -out:file:silent_struct_type binary -in:file:fullatom

After the pdbs have been scored use the included script to do an energy cut.

    make_sub_silent_file_percentile.py scored.out ecut_10.out -1 0.10

where 0.10 is a 10% energy cut.

 
# Nucleobase Sample Around

## Author
Rhiju Das, rhiju@stanford.edu

## Brief Description

Make tables of interaction energies between an adenosine nucleobase and, say, 
 a simple carbon atom or phosphate as a probe.

## Abstract

We wanted to compare potentials of mean force of various atoms like a single methyl probe, water, adenosine, etc. around a fixed adenosine to explicit molecular dynamics solutions.

Developed in summer 2012, Das Lab hackathon -- Kyle Beauchamp, Fang-Chieh Chou, Rhiju Das, Parin Sripakdeevong. 
Extended to include phosphate by Rhiju Das in Dec. 2014.

## Running
To sample a 'carbon' probe atom:
```
 nucleobase_sample_around   [-s a_RNA.pdb]
```

To sample a water (sampling all possible orientations and outputting Boltzmann summed free energies)
```
 nucleobase_sample_around   [-s a_RNA.pdb]  -sample_water  [ -extra_res ~/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/water/TP3.params ]
```
To sample a nucleobase
```
 nucleobase_sample_around   [-s a_RNA.pdb]  -sample_another_nucleobase   -copy_nucleobase_nucleobase_file double_A_ready_set.pdb
```
To sample an nucleobase, reading in a starting nucleobase-nucleobase pairing conformation.
```
 nucleobase_sample_around   [-s a_RNA.pdb]  -sample_another_nucleobase   -copy_nucleobase_nucleobase_file double_A_ready_set.pdb
```
Can now sample phosphates with the flags
```
 nucleobase_sample_around -sample_phosphate [-center_on_OP2]
```
The phosphate center is on the phosphorus atom, unless user specifies -center_on_OP2 . 
Note that due to some silliness in available variant types and the desire to use a phosphate from an actual nucleotide residue_type, the probe phosphate also has a floating C1'.

Recently added option
```
 -nucleobase g [or a,c,u]
```
which will use something other than adenosine as the central nucleotide to sample around

## Plotting Results

The plotting script is available in Rosetta/tools/rna_tools/pdb_util/plot_contour.py
```
plot_contour.py score_xy_0.table score_xy_0.png
plot_contour.py score_xy_1.5.table score_xy_1.5.png
plot_contour.py score_xy_3.table score_xy_3.png
plot_contour.py score_xy_4.table score_xy_4.png
```

Making non-canonical amino acid rotamer libraries
=================================================

This demo was written by P. Douglas Renfrew (renfrew at nyu dot edu).

This demo illustrates the make_rot_lib protocol. It was originally written as 
the 2nd (of 4) parts of a "protocol capture" acompanying the paper "Using 
Noncanonical Amino Acids in Computational Protein-Peptide Interface Design" by 
P. Douglas Renfrew, Eun Jung Choi, Brian Kuhlman (2011) PLoS One, from the 
RosettaCon2010 special PLoS ONE special issue. The four separate protocol 
captures that describe the four different ways Rosetta was used in the 
publication: creating Rosetta ResidueType parameter files for NCAAs, creating 
backbone dependent rotamer libraries for NCAAs, calculating explicit unfolded 
state reference energies for NCAAs, and the running of the 
DougsDockDesignMinimizeInterface protocol using NCAAs. The protocol caputres 
for creating ResidueType paramater files, rotamer libraries and explicite 
unfolded state energies describe the process used in the publication but are 
written from the standpoint of a researcher looking to add an additional NCAA 
to Rosetta. 

Creating a Noncanonical Amino Acid (NCAA) rotamer library is the second of two 
steps toward being able to use a NCAA in Rosetta. To add a new NCAA or to 
better understand how the NCAAs in the related publication were added one 
should have already completed or understand the steps in 
HowToMakeResidueTypeParamFiles.

Rotamer libraries are sets of common side chain conformations that generally 
correspond to local minima on the side chain conformational energy landscape. 
Side chain conformations are usually represented as a set of mean angles and a 
standard deviation to indicate variability. Rotamer libraries are used in for 
two main purposes in Rosetta: to provide starting points for side chain 
optimization routines, and the relative frequency is used as a pseudo-energy. 
Traditionally rotamer libraries are created by collecting statistics from 
protein structures. Rosetta uses the backbone dependent Drunbrack rotamer 
libraries. Since there are not enough structures containing NCAAs they must be 
generated.

Running the MakeRotLib protocol consists of four steps:

1. Creating and input template and generating the MakeRotLib options files.
2. Running the MakeRotLib protocol on each option file.
3. Assembling the individual rotamer libraries in a single file.
4. Modifying the ResidueType parameter file to be aware of our new rotamer 
   library.

Step 1: Making input files
--------------------------

Rosetta primarily uses backbone dependent rotamer libraries. Backbone-dependent 
rotamer libraries list provide side chain conformations sampled from residue 
positions whose backbone dihedral angles fall in particular bins. In the case 
of the Drunbrack rotamer libraries used by Rosetta the bins are in 10 degree 
intervals for for both phi and psi for a total of 1296 (36x36) phi/psi bins. 
To replicate this for the NCAAs we need to create a set of side chain rotamers 
for each member of a set of phi/psi bins.

The MakeRotLib protocol takes an option file as input. It requires an options 
file for each phi/psi bin. The first step in running it is creating these 1296 
options files. Continuing from the HowToMakeResidueTypeParamFiles protocol 
capture we are again using ornithine as an example. Ornithine has 3 sidechain 
dihedral angles (chi). We want to sample each chi angle from 0 to 360 degrees 
in 30 degree intervals, and based on the chemistry of the side chain we predict 
that were will probably be three preferred angles for each chi angle at 60, 
180, and 300 degrees for a total of 27 rotamers (3x3x3). We setup our 
MakeRotLib options file template as shown bellow.

    <<<<< C40_rot_lib_options_XXX_YYY.in start >>>>>
    AA_NAME C40
    PHI_RANGE XXX XXX 0
    PSI_RANGE YYY YYY 0
    NUM_CHI 3
    CHI_RANGE 1 0  330  30
    CHI_RANGE 2 0  330  30
    CHI_RANGE 3 0  330  30
    CENTROID 300 1 300 1 300 1
    CENTROID 300 1 300 1 180 2
    CENTROID 300 1 300 1  60 3
    CENTROID 300 1 180 2 300 1
    CENTROID 300 1 180 2 180 2
    CENTROID 300 1 180 2  60 3
    CENTROID 300 1  60 3 300 1
    CENTROID 300 1  60 3 180 2
    CENTROID 300 1  60 3  60 3
    CENTROID 180 2 300 1 300 1
    CENTROID 180 2 300 1 180 2
    CENTROID 180 2 300 1  60 3
    CENTROID 180 2 180 2 300 1
    CENTROID 180 2 180 2 180 2
    CENTROID 180 2 180 2  60 3
    CENTROID 180 2  60 3 300 1
    CENTROID 180 2  60 3 180 2
    CENTROID 180 2  60 3  60 3
    CENTROID  60 3 300 1 300 1
    CENTROID  60 3 300 1 180 2
    CENTROID  60 3 300 1  60 3
    CENTROID  60 3 180 2 300 1
    CENTROID  60 3 180 2 180 2
    CENTROID  60 3 180 2  60 3
    CENTROID  60 3  60 3 300 1
    CENTROID  60 3  60 3 180 2
    CENTROID  60 3  60 3  60 3
    <<<<< C40_rot_lib_options_XXX_YYY.in end >>>>>

Field definitions:

* AA_NAME \<three letter code for the amno acid> 

* PHI_RANGE \<phi value for this bin> \<phi value for this bin> 0 : The phi range 
  functionality is not functional. Both values need to be the same and the 
  interval set to 0

* PSI_RANGE \<psi value for this bin> \<psi value for this bin> 0 : The psi range 
  functionality is not functional. Both values need to be the same and the 
  interval set to 0

* NUM_CHI \<number side chain dihedral angles> : This should be the same as in 
  the parameter file.

* CHI_RANGE \<chi number> \<starting value> \<ending value> \<interval> : The 
  number of CHI_RANGE fields needs to equal the values specified for NUM_CHI.

* CENTROID \<Rotamer number for chi 1> \<starting value> {\<rotamer number for chi 
  2> \<starting value>}{etc.} : CENTROIDS specify the starting points for the 
  K-means clustering described in the related publication. A CENTROID field is 
  needs for each potential rotamer. The number of CENTROID fields defines the 
  number of rotamers listed in the resulting rotamer library.

To generate the 1296 input files we use a provided script that simply replace 
the XXX and YYY with the phi and psi values. The script is run as shown bellow.

    $ cd inputs
    $ ../scripts/make_inputs C40_rot_lib_options_XXX_YYY.in

The number of chi angles and the CHI_RANGE sampling interval are the primary 
determinants of the run time as they determine the number of rotamers that will 
be tested for each phi/psi bin. It is recommended to have at least 500 samples 
per chi. In the ornithine example we sample in 30 degree intervals for each of 
the 3 chi angles giving us a total of 1728 (12x12x12) conformations tested for 
each phi/psi bin. For a residue with a single chi 1 degree bins will suffice. 

Step 2: Running the MakeRotLib protocol
----------------------------------------

The next step is to run the MakeRotLib protocol on each of the input files we 
created in step one. This is the most time consuming portion of the process and 
should probably be done on a cluster. As cluster setups vary, an example for a 
single MakeRotLib options file is provided. The other 1295 should be run 
identically.

    $ cd outputs
    $ PATH/TO/ROSETTA/bin/make_rot_lib.macosgccrelease -database PATH/TO/rosetta_database -rot_lib_options_file ../inputs/C40_rot_lib_options_-60_-40.in >& C40_rot_lib_options_-60_-40.log &

NOTE: The extension on your executable maybe different.

The only options passed to the executable are the path to the database and the 
MakeRotLib options file. After the run completes a file called 
C40_-60_-40_bbdep.rotlib should be in the output directory. This is the 
backbone dependent rotamer library for a phi of -60 and a psi of -40.

The log file from the rosetta run in includes quite a bit of useful output. 
There are three main sections to the log output,  "ROTAMERS", "CENTROIDS" and 
"FINAL ROTAMERS" sections. Each one shows the following data: phi, psi, omega 
and epsilon backbone dihedral angles, probability, total energy, torsion 
energy, intra-residue repulsive, intra-residue attractive, the number of chi 
angles, the assigned cluster number, the set of input chi angles, the set of 
minimized chi angles, the standard deviation, and the distance from that point 
to each of the cluster centroids. The log file also displays the number of 
conformations per cluster, the average distance between the cluster center and 
the members of that cluster. Lack of conformations in a cluster and a large 
(>30) average cluster centroid distance suggests that that cluster is higher in 
energy. 

Step 3: Assembling the individual rotamer libraries in to a single file
-----------------------------------------------------------------------

After the MakeRotLib protocol has been run on all of the MakeRotLib options 
files the individual rotamer libraries for each phi psi bin need to be 
assembled in to a single file. This is accomplished with a provided script as 
shown bellow. 

    $ cd outputs
    $ ../scripts/make_final_from_phi_psi.pl C40

The single file rotamer library should be called C40.rotlib. The file should be 
placed in the ncaa_rotlibs directory in the database. 

    $ cp C40.rotlibs PATH/TO/DATABASE/rosetta_database/ncaa_rotlibs/

Step 4: Modifying the residue type PARAMS file
-----------------------------------------------

The last step is modifying the residue type parameter file to use the new 
rotamer library. To do this we need to add the name of the rotamer library, the 
number chi angles it describes, and how many bins there are for each chi angle 
to the Residue type parameter file. The ornithine rotamer library is called 
C40.rotlib and the rotamer library describes 3 chi angles and each of those 3 
chi angle has 3 rotamer numbers. So we would use the following commands to add 
that information to the file we created in HowToMakeResidueTypeParamFiles.

    $ echo "NCAA_ROTLIB_PATH C40.rotlib" >> PATH/TO/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/l-ncaa/ornithine.params
    $ echo "NCAA_ROTLIB_NUM_ROTAMER_BINS 3 3 3 3" >> PATH/TO/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/l-ncaa/ornithine.params

Calculate Protein Protein ΔΔG
=============================

A common computational problem involves finding the binding energy of a
protein-protein complex. In this tutorial, we will calculate the change in 
binding energy caused by point mutations in the complex. 

We will use the ubiquitin ligase E3a - ubiquitin conjugating enzyme (UBC)
complex as an example. For instance, when we introduce a mutation in the E3a 
enzyme, the binding reactions for the wild-type and the mutant are as the 
following:

    (1) E + P --> EP    Delta_G   (wild-type or WT)

    (2) E' + P --> E'P  Delta_G'  (mutant or WT)

where E represents the E3a ligase and P represents the UBC protein. 
The binding energy change due to the mutation can be obtained by: 

    (3) Delta_Delta_G = Delta_G' - Delta_G

To obtain Delta_G in Rosetta, we calculate the energy of the complex
using the Rosetta scoring function. We then score the protein and ligase
that have been pulled apart. Subtracting these two scores gives 
the binding energy of the complex, Delta_G.

To obtain a change in binding energy, we create point mutations using resfiles
(explained below). Subtraction of the wild-type Delta_G from the mutant
Delta_G' yields the change in binding energy due to the point mutations.

In this example we will use Rosetta Scripts to peform all of the operations in
a stand-alone application. Using a protocol defined in an XML script we can
perform operations on the complex in a stepwise manner. The protocol will:

1. Relax the input structure to relieve possible clashes in the PDB.
2. Repack the structure.
3. Calculate the Deleta_G of the wild type complex. 
4. Repack the structure with a resfile to make a point mutation.
5. Calculate the Delta_G' of the mutated complex.
6. Using Equation (3) to obtain the change in binding energy. The energy value
is in Rosetta energy unit. 

Resfiles contain information for how the packer should behave, such as
telling the packer to make a point mutation. For full Resfile documentaion, 
see:

http://graylab.jhu.edu/Rosetta.Developer.Documentation/all_else/d1/d97/resfiles.html

In root directory of this demo, you will find a sample resfile with the 
following content. It ask Rosetta to mutate residue 641 on chain A (E3a ligase) 
into a tryptophan.  The contents of this resfile are:

    NATAA
    USE_INPUT_SC
    EX1 EX23
    start 
    641 A PIKAA W

Running the demo
----------------

    rosetta_scripts.operatingsystem.release -parser:protocol mutation_script.xml -s ../starting_files/1C4Z.pdb -ignore_unrecognized_res -database /path/to/database -out:path:pdb ../output_files -out:path:score ../output_files -nstruct 1

The command line arguments are explained below. 

* operatingsystem:
  indicates the platform of the user. 

* `-parser:protocol`:
  this flag indicates the XML file that contains the point mutation ddg protocol.

* `-s`:
  input PDB file.

* `-ignore_unrecognized_res`:
  this flag ignores lines in the PDB file that Rosetta doesn't recognize.

* `-database`:
  a path to the Rosetta database shipped with the release.

* `-out:path:pdb`:
  indicates the output directory for pdbs.

* `-out:path:score`:
  the output directory for the score file.

* `-nstruct`:
  the number of models to output.
  To get an accurate representation of the energy landscape, many models should be created (e.g., >1000).
  During testing, you may set nstruct to 1. 

The output file will be a repacked and mutated PDB file with the dg_wt and
dg_mut lines at the end of the file. Subtration of these two numbers yields
the change in binding energy. 




# Rosetta VIP

This README was written in Feb. 2012, by Jim Havranek (havranek@genetics.wustl.edu).

This demo illustrates a protocol to identify candidate mutations for stabilizing proteins that have sub-optimally packed cores.

It has been published in the paper "Automated selection of stabilizing mutations in designed and natural proteins" by B. Borgo J.J. Havranek (2012), Proc. Natl. Acad. Sci. USA 109(5) pp 1494-99.

The source code is found at:  `(rosetta_source)/src/apps/public/vip.cc`

The app can be run as follows:

```
./(executables_path)/vip.(arch)(mode) -database (database_path) -s input.pdb -cp:ncycles 3 -cp:cutoff 6.0 -sasa_calculator_probe_radius 1.0 -run:silent
```

options are:

```
-cp:ncycles (Size)
	This will run the iterative protocol (find point mutations, relax, output best relaxed pose) a fixed number of times. If you don't use this option it will continue to run until it no longer finds favorable mutations. The latter option can take a while for large proteins.

-cp:max_failures (Size)
	This allows you to try each iteration multiple times before deciding it is pointless.  The reason for this option is that, due to RosettaHoles, the process is stochastic.  Whether a good mutation can be found can depend on the pseudorandom numbers pulled during the void identification steps.  The highest this should be set is around five.  The default is 1.

-cp:cutoff (real)
	This is the cutoff for choosing mutatable residues (ie distance from cavity ball to a non-bb, non-surface atom on the residue). The smallest cutoff you can use is best since that will mutated the smallest # of residues that line the cavity.

-sasa_calculator_probe_radius (real) Increasing this will likely give you surface clefts in addition to buried voids.
```

further options for the fast relax mover, sasa metric options, etc.:

```
-cp:pack_sfxn (score function) (e.g. -cp:pack_sfxn score12_full)
	Allows you to use a different score for the point mutant trials

-cp:relax_sfxn (score function) to use a different score for the relax stage

-cp:skip_relax
This causes the protocol to skip the relax step, which is a quick fix until I can exchange the relax mover for something more efficient. I'd only caution that scoring of the fixed bb point mutations isn't quantitative wrt to magnitude of the ddE. Rather, it seems to be very accurate in terms of favorable vs. unfavorable.

-cp:relax_mover
This controls the mover used for the relaxation step after inserting mutations.  The default is "relax", which performs the fast_relax protocol.  As available is "classic_relax", which can be a little slower.

-cp:local_relax
This restricts the scope of the relaxation protocol for possible mutations to 'neighbor' residues, which are defined as those residues within a X Ang. CB-CB distance.

-cp:print_reports
This generates a file (reports.txt ) that tells you which mutations were identified by the simple geometric screen, and which mutations still look good after relaxation.

-cp:print_intermediate_pdbs
This option outputs a fullatom pdb for each accepted mutation, named "vip_iter_1.pdb", with the number incremented for each pass through.  The pdb for the final iteration will be the same as the final pdb.

-cp:exclude_file
This allows you to specify a file that contains positions that should not be allowed to mutate.  The format is one position per line, with pdb number and chain separated by a space ("128 A").
```

a couple of other notes:

- It is preferable to let the application find as many mutations as possible.  This is accomplished by _not_ using the -cp:ncycles option.  If the application takes too long to run, the -cp:ncycles option can be used to limit the run time to an acceptable amount.

- Extra rotamers are automatically included in the point mutant trials, so if you use -ex1 -ex2 etc flags, these will be applied in the relax step which slows things down quite a bit.

- Cavity finding is done via RosettaHoles1, which is stochastic, so:  1) if it finds no cavities, run it again and it probably will and 2) separate runs can result in a different sequence of mutations.



# Instructions for Calculating Explicit Unfolded State Energies for Noncanonical Amino Acids

This README was written by P. Douglas Renfrew (renfrew@nyu.edu).

This demo illustrates the UnfoldedStateEnergyCalculator protocol. It was originally written as the 3rd (of 4) parts of a "protocol capture" accompanying the paper "Using Noncanonical Amino Acids in Computational Protein-Peptide Interface Design" by P. Douglas Renfrew, Eun Jung Choi, Brian Kuhlman (2011) PLoS One, from the RosettaCon2010 special PLoS ONE special issue. The four separate protocol captures that describe the four different ways Rosetta was used in the publication: creating Rosetta ResidueType parameter files for NCAAs, creating backbone dependent rotamer libraries for NCAAs, calculating explicit unfolded state reference energies for NCAAs, and the running of the DougsDockDesignMinimizeInterface protocol using NCAAs. The protocol captures for creating ResidueType parameter files, rotamer libraries and explicit unfolded state energies describe the process used in the publication but are written from the standpoint of a researcher looking to add an additional NCAA to Rosetta. 

Calculating the explicit unfolded state energies is the third of three steps toward being able to use a noncanonical amino acid (NCAA) in Rosetta. To add a new NCAA or to better understand how the NCAAs in the related publication were added one should have already completed or understand the steps in HowToMakeResidueTypeParamFiles and HowToMakeRotamerLibraries. 

The explicit unfolded state energies of an amino acid represent the energy of an amino acid in the unfolded state of a protein and is used to replace the reference energies in Rosetta. The UnfoldedStateEnergyCalculator uses a fragment based method to calculate the average unfolded state energies for each ResidueType. The protocols works on a large set of protein structures that are split in to randomly generated fragments. The central residue of each fragment is mutated to the residue of interest. The fragment is repacked. The unweighted energy for each energy method in the scoring function is recorded for the central residue. After the energies for all fragment central residues are collected, a boltzmann-weighted-average average energy is calculated for each term. 

Calculation explicit unfolded state energies for a NCAA requires three steps:
 - Obtaining a set of input pdbs
 - Running the UnfoldedStateEnergyCalculator protocol on the set of pdbs
 - Modifying the unfolded state energies file in the database

# Step 1: Obtaining a Set of Input PDBs

Since the UnfoldedStateEnergyCalculator protocol uses fragments from protein structures, we need a set of high quality structures to work with. Through their PISCES server, the Dunbrack laboratory maintains lists of structures in the Protein Data Bank organized based on xray resolution, precent sequence similarity, and r-factors\*\*. These lists are a convenient way to get a set of high quality structures. In this example we will use a list culled on May 20, 2011. It contains 1801 pdb files that have an xray resolution of at least 1.6 angstroms, less than 20% sequence identity, and r-factors of less than 0.25. To get the pdbs simply use a supplied script to download the pdbs from the Protein Data Bank ftp servers. 

```
$ cd inputs
$ ../scripts/get_pdbs.bash cullpdb_pc20_res1.6_R0.25_d110520_chains1859
```

There should be 1801 gzipped pdb files and a text file containing a list of them called cullpdb_pc20_res1.6_R0.25_d110520_chains1859_list_pruned in the inputs directory. Rosetta will sometimes fail to correctly read in particular pdbs files. The cullpdb_pc20_res1.6_R0.25_d110520_chains1859_list_pruned file is a list of the pdbs which have been screened to be read successfully by Rosetta. 

\*\*Citation: G. Wang and R. L. Dunbrack, Jr. PISCES: a protein sequence culling server. Bioinformatics, 19:1589-1591, 2003. 

# Step 2: Running the UnfoldedStateEnergyCaclulator Protocol

The UnfoldedStateEnergyCalculator is relatively easy to run. The command line options are described below:

- frag_size: single integer value, sets the number of residues in each fragment, should be an odd number and has a default of 5 which is what was used in the accompanying publication
- residue_name: string value, sets the three letter code of the residue type which the central residue will be mutated to
- repack_fragments: boolean value, controls if the fragments will be repacked before scoring and defaults to true
- native_sequence: boolean value, controls if the central residue will be mutated before scoring and defaults to false

Additionally it is strongly recommended to add the following flags as they will make Rosetta handle more pdb files and improves runtime by disabling default features that will be negated by the fragmenting and prepacking:

- ignore_unrecognized_res: causes Rosetta to ignore unrecognized residue types and 
- ex1 and ex2 and extrachi_cutoff 0: force rosetta to use additional rotamer during the fragment repacking
- mute all and unmute devel.UnfoldedStateEnergyCalculator and unmute protocols.jd2.PDBJobInputer: reduces the size of the log file significantly by turning off unnecessary output
- no_optH true: turns off the hydrogen optimization done when the protein is first read in 
- detect_disulf false: turns off disulfide detection

Continuing the ornithine example we have used in the two previous protocol captures, to calculate the unfolded state energies one would run the following command.

```
$ cd outputs
$ PATH/TO/bin/UnfoldedStateEnergyCalculator.macosgccrelease -database PATH/TO/rosetta_database -ignore_unrecognized_res -ex1 -ex2 -extrachi_cutoff 0 -l ../inputs/cullpdb_pc20_res1.6_R0.25_d110520_chains1859_list_pruned -residue_name C40 -mute all -unmute devel.UnfoldedStateEnergyCalculator -unmute protocols.jd2.PDBJobInputer -no_optH true -detect_disulf false >& ufsec_log_c40 &
```

**Note:** The extension on your executable my be different.

The run will take between 30-60 seconds per pdb file.

The log file contains lots of useful information. It contains the unweighted energies for each of the energy methods for each of the individual fragments. At the end it will print the average unweighted energies for each ResidueType as well as the Boltzmann weighted average unweighted energies. Boltzmann weighted average unweighted energies are used because some backbones just can't tolerate a mutation to a particular ResidueType and there are extremely high repulsive energies for some fragments that skew the average value. Using the Boltzmann weighting removes the higher energy outliers in a more elegant fashion than a hard energy cutoff.

# Step 3: Modify the Unfolded State Energies File

Once the UnfoldedStateEnergyCalculator has finished running the Boltzmann weighted average unweighted energies need to be added to the database. The line you want is the "BOLZMANN UNFOLDED ENERGIES". These are the Boltzmann weighted average unfolded energies for each energy method. The file you need to modify is unfolded_state_residue_energies_mm_std.

Using the ornithine line as an example, the line form the log file is... 

```
BOLZMANN UNFOLDED ENERGIES:  fa_atr:    -2.462 fa_rep:     1.545 fa_sol:     1.166 mm_lj_intra_rep:     1.933 mm_lj_intra_atr:    -1.997 mm_twist:     2.733 pro_close:     0.009 hbond_sr_bb:    -0.006 hbond_lr_bb:     0.000 hbond_bb_sc:    -0.001 hbond_sc:     0.000 dslf_ss_dst:     0.000 dslf_cs_ang:     0.000 dslf_ss_dih:     0.000 dslf_ca_dih:     0.000
```

We could add the following to the unfolded_state_residue_energies_mm_std file in the database using the command bellow.

```
$ echo "C40 -2.462 1.545 1.166 1.933 -1.997 2.733 0.009 -0.006 0.000 -0.001  0.000" >> minirosetta_database/scoring/score_functions/unfolded/unfolded_state_residue_energies_mm_std 
```

The ResidueType can now be used in almost any Rosetta protocol that is compatible with the MM_STD scoring function.
Doug's DockDesignMinimize Demo
==============================

This demo was written by P. Douglas Renfrew (renfrew at nyu dot edu)

This demo illustrates the DougsDockDesignMinimizeInterface protocol. It was 
originally written as the 4th (of 4) parts of a "protocol capture" accompanying 
the paper "Using Noncanonical Amino Acids in Computational Protein-Peptide 
Interface Design" by P. Douglas Renfrew, Eun Jung Choi, Brian Kuhlman (2011) 
PLoS One, from the RosettaCon2010 special PLoS ONE special issue. The four 
separate protocol captures that describe the four different ways Rosetta was 
used in the publication: creating Rosetta ResidueType parameter files for 
NCAAs, creating backbone dependent rotamer libraries for NCAAs, calculating 
explicit unfolded state reference energies for NCAAs, and the running of the 
DougsDockDesignMinimizeInterface protocol using NCAAs. The protocol captures 
for creating ResidueType parameter files, rotamer libraries and explicit 
unfolded state energies describe the process used in the publication but are 
written from the standpoint of a researcher looking to add an additional NCAA 
to Rosetta.

The DougsDockDesignMinimize (DDDM) protocol was used in the accompanying 
manuscript to redesign the protein/peptide interface of Calpain and a fragment 
of its inhibitory peptide calpastatin. The protocol was written for this 
specific protein/peptide interaction and modifications to the code will be 
necessary to run the protocol on a different system. A modified form was used 
as the example protocol in the advanced section of the Rosetta 3.0 release 
manual.

Minor modifications to the protocol have been made from the version used to 
produce the designs in the accompanying publication. The interfaces of the 
Rosetta libraries have changed since the initial implementation of the protocol 
and these modifications were necessary to allow the protocol to work with the 
current release. 

The protocol has two loops referred to as the inner and outer loops. The outer 
loop controls the number of structures generated, outputting structures only if 
certain filters are met. The inner loop iterates between two phases: a 
perturbation phase and a design phase. During the perturbation phase three 
types of perturbations are used. Two types of perturbations are applied to the 
peptide: rotational and translational rigid body permutation of the peptide in 
the binding pocket of the protein, and small and shear perturbations of the 
backbone dihedral angles in all residues of the peptide. The third type of 
perturbation is preformed on residues 1-45 of the N-terminus which comprise the 
peptide binding site and surrounding residues of the protein are perturbed 
using small and shear perturbations of the backbone dihedral angles. One of the 
three perturbations is randomly chosen to be used each time the function in 
called and is followed by a coarse side chain optimization (RotamerTrials). 
During the design phase of the protocol the side chains positions are optimized 
using a more intensive side chain optimization routine (PackRotamers) followed 
by minimization of the backbone and side chain dihedrals angles as well as the 
jump between the peptide and protein. The number of perturbations per design as 
well as the magnitude of perturbations are controlled using command line 
options described bellow.

There are three basic steps to running the protocol as it was used in the 
accompanying manuscript.

 - Generating the input resfiles and folders
 - Running the DougsDocDesignMinimize protocol for each mutation at each position 
 - Running the analysis scripts

Example output is provided for a small version of the full run.

Generating input resfiles and folders
-------------------------------------

The Protein Databank (PDB) code for the Calpain/Calpastatin structure used in 
the designs is 1NX1. The files contain two copies of calpain (chains A and B) 
and two copies of the calpastatin inhibitory peptide (chains C and D). To 
reduce the size of the system, designs were done on chain A and chain C only 
because chains A and C have lower B-factors than chains B and D. The 
calpain/calpastatin interface is distal to the calpain/calpain interface and 
the resides that make up the calpain/calpain interface were held fix during the 
protocol. Additionally the calcium atoms and water molecules were removed. 
Rosetta doesn't work with with water molecules and the calcium ions are not 
near the protein peptide interface. This input pdb was repacked and both the 
side chain and backbone dihedrals minimized using a modified version of the 
fixbb app. The input PDB file is called 1NX1_clean_repack_min_all.pdb

Each of the NCAAs added in the accompanying publication was tried at each 
position in the peptide. To do this and to keep all of the output pdb organized 
a script is provided that creates folders and generates resfiles based on 
templates for each sequence. To generate the resfiles and folders run the 
following commands...

    $ cd run_dir
    $ ../scripts/make_folders_resfiles.tcsh

NOTE: To save space, the script has been modified to only produce resfiles and 
folders for position 610 in the peptide and residue type MPA 
(4-methyl-phenylalanine). To modify the script to produce folders and resfiles 
for each position simply uncomment the lines that read `#foreach i ( 601 602 
603 604 605 606 607 608 609 610 611 )` and `#foreach j ( ABA APA HLU ..... C92 
C93 C94 )` comment out the line that reads `foreach i ( 610 )` and `foreach j ( 
MPA )`.

Using the templates in the run_dir the make_folders_resfiles.tcsh script makes 
folders and resfiles to preform all of the DDDMI runs.


Running DougsDockDesignMinimizeInterface protocol
-------------------------------------------------

To run the protocol modifications need to made to files in the database. In the 
file rosetta_database/chemical/residue_type_sets/fa_standard/residue_types.txt 
all of the paths to the residue type parameter files under the L-NCAA heading 
need to be uncommented by removing the "#" from the front of the line. 
Additionally the rotamer libraries for the NCAA are not provided in the default 
Rosetta database because they are more than 400MB. The rotamer libraries for 
the NCAAs added in the accompanying publication are provided as supplemental 
information. 

NOTE: Turning on all of the additional residue types dramatically increases the 
number of residue types and the memory footprint of Rosetta. The memory foot 
print can be reduced by commenting out unnecessary patches in the 
rosetta_database/chemical/residue_type_sets/fa_standard/patches.txt file. For 
the DougsDockDesignMinimizeProtocol all but the NtermProteinFull.txt and 
CtermProteinFull.txt can be safely commented out by placing "#" symbols at the 
beginning of each line of the patches.txt except for the lines that say 
"NtermProteinFull.txt" and "CtermProteinFull.txt".

To run the protocol as in the accompanying publication preform the following 
commands starting at the HowToRunDougsDockDesignMinimizeProtocol directory.

    $ cd run_dir
    $ cd pos_610_MPA
    $ /PATH/TO/ROSETTA/bin/doug_dock_design_min_mod2_cal_cal.macosgccrelease -database /PATH/TO/rosetta_database -s ../../inputs/1NX1_clean_repack_min_all.pdb -resfile ../resfile_pos_603_MPA -nstruct 255 -inner_num 45 -pert_num 25 -ia_ener 100 -use_input_sc -pdb_gz

The above command generates 255 structures and will take approximately 5 
minutes per structure depending on your hardware.

NOTE: The extension of your executable maybe different than the above. Also in 
the publication the nstruct command line option was 255. However to save space 
for the protocol capture an nstruct of 10 was used.

A script is provided that will preform the above command for each folder 
created by the make_folders_resfiles.tcsh script.

    $ cd run_dir
    $ ../scripts/run_script.bash

NOTE: You will need to set the path to your database and executable in the 
run_script.bash.

There are a number command line options to control the running of the 
application.

* `-pert_mc_temp`: The MC temperature to use for the perturbation phase of the DDDM protocol, defaults to 0.8 kT.
* `-pert_dock_rot_mag`: The rotation magnitude for the ridged body perturbation in the perturbation phase of the DDDM protocol, defaults to 0.5. 
* `-pert_dock_trans_mag`: The translation magnitude for the ridged body perturbation in the perturbation phase of the DDDM protocol, defaults to 0.25. 
* `-pert_pep_small_temp`: The temperature for the internal MC object in the small mover of the peptide perturbations, defaults to 0.8 kT. 
* `-pert_pep_shear_temp`: The temperature for the internal MC object in the shear mover of the peptide perturbations, defaults to 0.8 kT. 
* `-pert_ter_small_temp`: The temperature for the internal MC object in the small mover of the termini perturbations, defaults to 0.8 kT. 
* `-pert_ter_shear_temp`: The temperature for the internal MC object in the shear mover of the termini perturbations, defaults to 0.8 kT. 
* `-pert_pep_small_H`: The maximum angle of perturbation for helical secondary structure for the peptide small mover, defaults to 1 degree.
* `-pert_pep_small_L`: The maximum angle of perturbation for loop secondary structure for the peptide small mover, defaults to 1 degree.
* `-pert_pep_small_E`: The maximum angle of perturbation for strand secondary structure for the peptide small mover, defaults to 1 degree.
* `-pert_pep_shear_H`: The maximum angle of perturbation for helical secondary structure for the peptide shear mover, defaults to 1 degree
* `-pert_pep_shear_L`: The maximum angle of perturbation for loop secondary structure for the peptide shear mover, defaults to 1 degree.  
* `-pert_pep_shear_E`: The maximum angle of perturbation for strand secondary structure for the peptide shear mover, defaults to 1 degree.
* `-pert_ter_small_H`: The maximum angle of perturbation for helical secondary structure for the termini small mover, defaults to 1 degree
* `-pert_ter_small_L`: The maximum angle of perturbation for loop secondary structure for the termini small mover, defaults to 1 degree.  
* `-pert_ter_small_E`: The maximum angle of perturbation for strand secondary structure for the termini small mover, defaults to 1 degree.
* `-pert_ter_shear_H`: The maximum angle of perturbation for helical secondary structure for the termini shear mover, defaults to 1 degree
* `-pert_ter_shear_L`: The maximum angle of perturbation for loop secondary structure for the termini shear mover, defaults to 1 degree.  
* `-pert_ter_shear_E`: The maximum angle of perturbation for strand secondary structure for the termini shear mover, defaults to 1 degree.
* `-pert_pep_num_rep`: Number of small and shear iterations for the peptide, defaults to 100. 
* `-pert_ter_num_rep`: Number of small and shear iterations for the terminus, defaults to 100. 
* `-pert_num`: Number of iterations of perturbation loop per design, defaults to 100.
* `-inner_num`: Number of iterations of the inner loop, defaults to 100. 
* `-ia_ener`: Upper energy limit for final design/interface analysis checkpoint, defaults to 0.0. 
* `-desn_mc_temp`: The temperature to use for the design/minimization phase of the DDDM protocol, defaults to 0.8 kT.

For the most part the defaults should suffice and are what was used in the 
paper. The "ia_ener" command line option is dependent on the complex that the 
protocol is run on. Setting it unreasonably low will cause the protocol to run 
forever and never output any structures. Setting it to 100 above allows all but 
the most aggressions structures to be output. The nstruct command line option 
controls the outer loop. 

Analyzing the results
---------------------

Each of the output pdb files contains information about the protein peptide 
complex that can used to evaluate the designs.  For example at the end of the 
1NX1_clean_repack_min_all_0004.pdb is shown bellow. For each filter the value 
is calculated for the protein and peptide together (COMPLEX) and separated by 
1000 angstroms (SEPARATE) and the difference between the two (DIFF). ENERGY is 
the Rosetta energy. The ENERGY_COMPLEX is the primary determinant to how good a 
design is and the ENERGY_DIFF can give an estimate for the binding energy. SASA 
is the solvent accessible surface area. SASA_DIFF is indicative of sequences 
that make a more protein-peptide contacts and can be used for screening designs 
for example placing a very large side chain at a constrained interface position 
can cause the peptide to be pushed out of the binding pocket which would be 
reflected in a smaller magnitude SASA_DIFF. HB_ENER is the hydrogen bonding 
component of the Rosetta energy. Larger HB_ENER_DIFF values indicate that the 
design is making more or better hydrogen bonds across the protein peptide 
interface. PACK is the RosettaHoles score and is a measurement of how well the 
protein is packed. A PACK_COMPLEX that is larger than the PACK_SEPARATE is 
favorable and suggests that the complex is better packed than the protein 
alone. Additionally the RosettaHoles score penalizes holes that cannot be 
occupied by solvent so larger PACK_DIFF score indicate that the designed 
peptide is capable of filling cavities in the protein that are inaccessible to 
solvent.

    ENERGY_COMPLEX:	   -63.3501
    ENERGY_SEPERATE:   -48.599
    ENERGY_DIFF:	   -14.7512
    SASA_COMPLEX:	   11774.1
    SASA_SEPERATE:	   13092.3
    SASA_DIFF:	   -1318.15
    HB_ENER_COMPLEX:   -146.728
    HB_ENER_SEPERATE:  -145.037
    HB_ENER_DIFF:	   -1.69085
    PACK_COMPLEX:	   0.515277
    PACK_SEPERATE:	   0.466995
    PACK_DIFF:	   0.0482824

A script is provided that pulls the information out of a set of pdb files and 
sorts it based on the ENERGY_COMPLEX metric. 

    $ cd run_dir
    $ cd pos_610_MPA
    $ ../../scripts/get_interface_data.tcsh

The script produces a file call out.ALL that contains a single line for each 
pdb file with the metrics in the above order.
This folder contains material from a tutorial given at Baylor in January 2010,
on the use of Rosetta for structure modeling using low to medium resolution
cryoEM data.  This folder also contains a PDF distributed during the tutorial.

All the examples from the tutorial work with the current release, however, the
executible names and a few of the options have changed since the tutorial was
given.  Always refer to the scripts themselves rather than the command lines
given in the PDF file.

AbInitio Demo
=============

Run on linux
------------
    Rosetta/main/source/bin/AbinitioRelax.linuxgccrelease @flags

Run on macs
-----------
    Rosetta/main/source/bin/AbinitioRelax.macosgccrelease @flags

Flags
-----
    -in:file:fasta ./input_files/1elwA.fasta
    -in:file:frag3 ./input_files/aa1elwA03_05.200_v1_3
    -in:file:frag9 ./input_files/aa1elwA09_05.200_v1_3
    -in:file:native ./input_files/1elw.pdb

    -abinitio:relax
    -nstruct 1
    -out:pdb

    -use_filters true
    -psipred_ss2 ./input_files/1elwA.psipred_ss2
    -abinitio::increase_cycles 10
    -abinitio::rg_reweight 0.5
    -abinitio::rsd_wt_helix 0.5
    -abinitio::rsd_wt_loop 0.5
    -relax::fast

Analyze Output
--------------
The output_files directory contains example output.

In `score.fsc` get a score and RMS for each model.
Typical analysis makes a scatter plot of these with RMS on the x-axis and score on the y-axis.
Look for a "funnel" to low energies and RMS in a successful ab initio prediction.
A failed prediction will not have low RMS/energy structures.
For the following arguments for full production run (using the minirosetta compile):

    -abinitio::fastrelax
    -abinitio::increase_cycles 10
    -abinitio::rg_reweight 0.5
    -abinitio::rsd_wt_helix 0.5
    -abinitio::rsd_wt_loop 0.5
    -abinitio::use_filters false
    -database minirosetta_database
    -ex1
    -ex2aro
    -frag3 aa0000103_05.200_v1_3
    -frag9 aa0000109_05.200_v1_3
    -in:file:boinc_wu_zip ploop23_control_fold_data.zip
    -in:file:native 00001.pdb
    -mute all
    -mute all
    -out:file:silent default.out
    -relax::default_repeats 15
    -silent_gz

    resultfiles = default.out.gz queue = 3000

And run a relax:

    -database
    -ex1
    -ex2aro
    -frag3 aa0000103_05.200_v1_3
    -frag9 aa0000109_05.200_v1_3
    -in:file:boinc_wu_zip ploop23_control_fold_data.zip 
    -in:file:fullatom
    -in:file:native 00001.pdb
    -in:file:s 00001.pdb
     minirosetta_database
    -mute all
    -out:file:silent default.out
    -relax::default_repeats 15
    -run:protocol relax
    -silent_gz

    resultfiles = default.out.gz
Analyzing Interface Quality
===========================

Outline
-------
You have been provided with two PDBs to use for demonstration: 1U6E, a simple homodimer, and 3R2X, Sarel's hemaglutinin/designed protein structure.  In the latter case, we will calculate the interface between the designed protein and the HA - the interface between HA chains is unimportant.

Preparing the inputs
--------------------
We usually score the input pdbs to make sure they are able to be read by Rosetta and to replace any missing sidechain atoms.
From the main demo directory run the following command to score 

    cd rosetta_inputs
    /path/to/rosetta/bin/score_jd2.default.macosgccrelease -s ../starting_files/*.pdb.gz -no_optH false -database /path//to/rosetta_database/ -ignore_unrecognized_res -out:pdb

There are a few options here that need to be described:
* `-no_optH false`: This will make rosetta consider Q,N,H ring flips (this helps remove some buried unsatisfied polar atoms).
* `-ignore_unrecognized_res`: This will drop any residues and waters from the PDB that rosetta does not recognize.
* `-out:pdb`: This forces the score application to output a pdb to use later

This will also output file called score.sc which is the scored input structures.
Now we'll rename these PDBs to make this a bit easier to keep track of:

    mv 1u6e_0001.pdb 1u6e_scored.pdb; mv 3r2x_0001.pdb 3r2x_scored.pdb

Running InterfaceAnalyzer
-------------------------
Now we have all of our inputs set to go to run InterfaceAnalyzer.
If you are interested you should read the full documentation for this application located in `/path/to/rosetta/doc/apps/public/analysis/interface_analyzer.dox`:

    cd example_output

First we are going to look at the 1U6E homodimer interface.  It only has two chains so the default way of defining an interface will work. 
You will need to modify the options files to contain the path to your rosetta database.
To analyze this interface run the following commands:
This one will use only the input sidechains from the input structure and not repack them.  This probably only a good idea if you have already designed/minimized this structure, but is provided as an example here anyway

    /path/to/rosetta/mini/bin/InterfaceAnalyzer.default.macosgccrelease -s ../rosetta_inputs/1u6e_scored.pdb @../rosetta_inputs/no_pack_input_options.txt

This one will repack the input sidechains from the input structure before analyzing anything about the interface. This is a better idea in this case because we are using a raw pdb as input.

    /path/to/rosetta/mini/bin/InterfaceAnalyzer.default.macosgccrelease -s ../rosetta_inputs/1u6e_scored.pdb @../rosetta_inputs/pack_input_options.txt

You can look in the options files for a full description of each of the command line flags.

Now to analyze the 3R2X interface between the designed protein (Chain C) and Hemaglutinin (Chains A & B)  we need to consider chains A and B as one monomer and chain C as the binding partner. To do this we have the option available to keep any given chains together in this case the option needed is 

    -fixedchains A B

We will add this to the command lines given above; so now run these:

    /path/to/rosetta/mini/bin/InterfaceAnalyzer.default.macosgccrelease -s ../rosetta_inputs/3r2x_scored.pdb -fixedchains A B @../rosetta_inputs/no_pack_input_options.txt
    /path/to/rosetta/mini/bin/InterfaceAnalyzer.default.macosgccrelease -s ../rosetta_inputs/3r2x_scored.pdb -fixedchains A B @../rosetta_inputs/pack_input_options.txt

Now in the current directory (`example_output`) we have two different score files, one for when we repacked the interface (`pack_input_score.sc`), and one for when we kept the input rotamers fixed (`no_pack_input_score.sc`)
From here we can get some important info.  We will only concentrate on a few here. For a full description of what everything means see the documentiation page for this application: `/path/to/rosetta/doc/apps/public/analysis/interface_analyzer.dox`

Looking at the results
----------------------
Open the output files in some form of text editor or spreadsheet application.  Below we describe how to find a few important output values.
Binding energy (dG_separated): this is computed deltaG of binding. Notice that in no_pack_input the value for 3R2X is unrealistic.  This is probably due to clashes in the input structure, which is why it is a good idea to relax your structure somehow before calculating these values.
Number of buried unsatisfied polar atoms (delta_unsatHbonds): note that the number is higher for 1U6E when you pack the input, this sometimes happens when rosetta tries to relieve clashes and thus leaves some polars without a hydrogen bond partner
Packing score (packstat): This is a measure of how well packed the interface is with 0.0 being as poor as possible and 1.0 being perfect shape complementarity (usually values above 0.65 are good)
Buried Surface Area (dSASA_int): change in exposed surface area upon formation of an interface.
Binding energy per unit area (dG_separated/dSASAx100): This is the dG_separated binding energy divided by the total interface surface area (dSASA_int). We multiply by 100 to scale it up the value so it fits better in the score file.  We like using this to make sure that rosetta is making high quality contacts instead of making a lot of low quality contacts across the interface.  Usually values below -1.5 are pretty good. 

Molecular Replacement Demo
==========================

Written October 26, 2010  
Modified August 27, 2013

---

This document briefly walks through the use of Rosetta to solve difficult 
molecular replacement problems.  These tools assume that the user has access to 
the Phenix suite of crystallographic software (in particular, phaser and the 
mapbuilding script mtz2map); however, all intermediate files are included so 
that if the user does not, most of the demo may still be run.

The basic protocol is done in 5 steps; each step has a corresponding script in 
the folder:

1. Using HHSearch, find potential homology to the target sequence.  Use a 
   Rosetta "helper script" to prepare templates (and Rosetta inputs for 
   subsequent computations).

2. Use PHASER to search for placement of the trimmed templates within the unit 
   cell.

3. Generate a map correspoding to each putative MR solution.

4. Using Rosetta, rebuild gaps and refine each template/orientation in Rosetta, 
   constrained by the density of each solution.  After rescoring with PHASER, 
   the best template/orientation should be clear (if the correct solution was 
   among the starting models).

Step 1: prepare_template_for_MR.sh
----------------------------------

This command-line illustrates the use of my script for preparing templates for 
an initial phaser run.  Functionally, it's doing the same thing as the 
crystallographic software 'Sculptor' but it doesn't remap the residues as 
sculptor does (and makes it easier to run with different alignments).  The 
script takes just one arguments: an HHR format alignment file.

Alignments generally come from HHsearch's web interface 
(http://toolkit.tuebingen.mpg.de/hhpred).  After submitting the sequence 
through their website, export the results to a .hhr file.  Results may be 
trimmed so only alignments with a reasonable e-value and sequence coverage are 
included.

The script parses the .hhr file, downloads each template PDB, and trims the PDB 
to the aligned residues.  In addition, the script produces a 'rosetta-style' 
alignment file; the format is briefly introduced below.  These alignment files 
are used in Rosetta model-building.

    ## 1CRB_ 2qo4b
    # hhsearch
    scores_from_program: 0 1.00
    2 DFNGYWKMLSNENFEEYLRALDVNVALRKIANLLKPDKEIVQDGDHMIIRTLSTFRNYIMDFQVGKEFEEDLTGIDDRKCMTTVSWDGDKLQCVQKGEKEGRGWTQWIEGDELHLEMRAEGVTCKQVFKKV
    0 AFSGTWQVYAQENYEEFLRAISLPEEVIKLAKDVKPVTEIQQNGSDFTITSKTPGKTVTNSFTIGKEAEIT--TMDGKKLKCIVKLDGGKLVCRTD----RFSHIQEIKAGEMVETLTVGGTTMIRKSKKI
    --

The first line is '##' followed by a code for the target and one for the 
template.  The second line identifies the source of the alignment; the third 
just keep as it is.  The fourth line is the target sequence and the fifth is 
the template ... the number is an 'offset', identifying where the sequence 
starts.  However, the number doesn't use the PDB resid but just counds residues 
_starting at 0_.  The sixth line is '--'.

The results for this demo appear in the folder 'templates'.  For each 
alignement in the starting .hhr file, 3 files are produced.

Steps 2 & 3: run_phaser.sh and make_maps.sh
-------------------------------------------

This command line shows the use of Phaser to generate initial molecular 
replacement solutions.  For each template we run phaser to find potential 
placements of each template in the unit cell.

The example scripts here only generate a single model from a single template, 
but for a real-world case, one will often want to use many different templates 
and may want to generate more than one possible solution using 'TOPFILES n'.  
In general, though, we have found it is better to use fewer potential solutions 
from more templates than many solutions from few templates.

Sometimes weak hits may be found by lowering the rotation function cutoff in 
phaser by adding the line 'SELECT ROT FINAL PERCENT 0.65' (or even 0.5) to the 
phaser script.  Increasing the packing function threshold (with PACK 10) may 
also help in some cases.

Finally, for each template/orientation, we generate the 2mfo-dfc map for input 
to Rosetta in the next step.

Steps 4A & 4B: run_rosetta_mr.sh
--------------------------------

The final step illustrate the use of rosetta's comparative modeling into 
density.  After running the script and an initial phaser run, density maps are 
generated from each phaser hit, and cm-into-density is done.  The flag 
-MR::mode cm is used to run this mode.  This first application does not try to 
rebuild gaps in the alignment, it just performs the threading and runs relax 
into density.  Thus, the only inputs needed are: the target fasta file, the 
rosetta-style ali file, and the template pdb.  Because there is no rebuilding, 
not many models are needed to adequately cover conformational space, generally 
10-20 is sufficient.

This script is the same as above, but also rebuilds gaps in the alignment.  The 
main difference is that a non-zero value is given for 
`-MR::max_gaplength_to_model`; additionally, some flags must be given that 
describe how rosetta should rebuild gaps.

Several additional input files must be provided as well.  Rebuilding of gaps is 
done by fragment insertion (as in Rosetta ab initio); thus two backbone 
fragment files (3-mers and 9-mers) must be given.  The application for building 
these is included with rosetta but requires a bunch of external 
tools/databases.  The easiest way to generate fragments is to use the Robetta 
server (http://robetta.bakerlab.org/fragmentsubmit.jsp).  The fragment files 
should be built with the full-length sequence; rosetta handles remapping the 
fragments if not all gaps are rebuilt.

A brief overview of flags is given below:

    -in::file::fasta inputs/1crb.fasta
    -in::file::alignment inputs/1crb_2qo4.ali
    -in::file::template_pdb inputs/2qo4.PHASER.1.pdb
        The fasta, alignment and template PDBs.  See section 1 for the input file format if it needs to be hand-edited.

    -edensity:mapfile inputs/sculpt_2QO4_A.PHASER.1.map
    -edensity:mapreso 3.0
    -edensity:grid_spacing 1.5
        This is how the density map and scorefunction parameters are given to Rosetta.  The input map (-edensity:mapfile) is CCP4 format.  The flags 'mapreso' defines the resolution of the calculated density; If the data is high-resolution it is often good to limit this to 2.5 or 3.  The grid spacing should be no more than 1/2 the map resolution; if -MR::fast is used (see below), then the grid_spacing flag may be omitted (since more finely sampled grids will have much less of a speed penalty).

    -MR::cen_dens_wt 4.0
    -MR::fa_dens_wt 1.0
        This controls the weight on the experimental density data furing the two stages of refinement.  The second flag (fa_dens_wt) has the greatest impact on the final models.  If, after model generation and visual inspection, models don't seem to be fitting the density well (or overfitting to the density), this may be adusted accordingly.  If omitted, the values shown abve are the defaults that are used; generally, these defaults are sufficient for many cases.

    -MR::fast
        A special faster density scoring formulation is used.  Off by default, but it is recommended.

    -MR::max_gaplength_to_model 8
        Rosetta will close gaps up to this width; the larger this value is, the more sampling is required.  Values higher than 10 will often return incorrect loop conformations, although for very restrained segments, or largely helical segments, large insertions may be successfully modeled.

    -nstruct 20
        The number of output structures.  Generally 10-20 is sufficient, unless a large 'max_gaplength_to_model' is given.

    -ignore_unrecognized_res
        If the template contains nonstandard residues/ligands/waters, this tells Rosetta to ignore them.  This flag is recommended.

    -loops::frag_files inputs/aa1crb_09_05.200_v1_3.gz inputs/aa1crb_03_05.200_v1_3.gz none
        (Optional) 	Fragment files from Robetta.  If omitted, MR-Rosetta will automatically generate fragments for the input structure; this may slightly reduce final model accuracy.  (The two separate command lines illustrate using and omitting this flag).


Since each model is independently generated, multiple processes may be used to 
produce all the necessary models.  To manage the output, either each process 
can be run from a separate directory, or '-out:prefix <prefix>' can be used to 
keep jobs from overwriting each other's structures.  Rosetta workloads may also 
be split using MPI; see the rosetta documentation for more details.

Alternatively, there is a compact output format, 'silent files' that can be used to dump structures to.  Simply add the flags '-out:file:silent <silent_filename> -out:file:silent_struct_type binary' and all structures from one process will be written to this compact file.  Then the rosetta program 'extract_pdb' can be used to extract:

    $ bin/extract_pdbs.default.linuxgccrelease -database $DB -in:file:silent <silent_filename> -silent_struct_type binary

Enzyme Design Demo
==================

Tutorial for a complete *de novo* enzyme design run, using the TIM reaction as an 
example, as published in 

* Richter F, Leaver-Fay A, Khare SD, Bjelic S, Baker D (2011) De Novo Enzyme 
  Design Using Rosetta3. PLoS ONE 6(5): e19230. 
  doi:10.1371/journal.pone.0019230

Tutorial written at RosettaCon2011 by Florian Richter (floric at uw dot edu), 
with help from Patrick Conway, Amanda Loshbaugh, Neil King, and Gert Kiss. 

The contents of the demo directory should be:

    starting_files/
     -- directory in which the raw input files are given - these
        are provided for you and serve as the base for your
        tutorial

    rosetta_inputs/
     -- directory in which the modified starting files should
        be placed which will be used as inputs to Rosetta.
        You may need to make modifications like stripping
        extra chains from the input PDB; store the modified
        PDB here and leave the unaltered one in starting_files

    scripts/
     -- python scripts, shell scripts, used in the workflow
     -- awk, grep, sed, etc. command lines

    README.dox
     -- A prose or list description of how to perform the protocol

    FOR_AUTHORS.txt
     -- A description for the demo creators of what their demo
        should achieve.
     -- Most of what starts in this file should end up in the
        README file as well.

Relevant documentation
----------------------

1. The above cited PLoS ONE paper,
2. Documentation for the enzyme design app  
   https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d6/dbc/enzyme_design.html
3. Documentation for the match app  
   https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d7/dfc/match.html)
4. Documentation about the enzdes cstfile format used for both matching and 
   enforcing catalytic geometries during design  
   https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d5/dd4/match_cstfile_format.html)
5. Familiarize yourself with how to generate a .params file and rotamer library 
   for your ligand of interest, as described in the ligand docking app 
   documentation.  
   https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d4/d47/ligand_dock.html  
   (in particular the section "Preparing the small-molecule ligand for docking")

The purpose of this tutorial is to reenact all steps described in the
PLoS ONE paper.

Step 1: Defining a theozyme in Rosetta format
---------------------------------------------
1.  Inputs required:

    A .params file for all residues used in the theozyme. The Rosetta3 database 
    has .params files for all amino acids, but they will have to be generated 
    for the ligand/reaction substrate. Refer to ligand docking documentation 
    under the link above.

2.  Outputs generated:

    A .cst (constraints) file. This file is written by the user and defines the 
    geometry of the theozyme in Rosetta format (see enzdes cstfile 
    documentation linked above) for use in subsequent steps of the design 
    process. 

3.  Defining a theozyme:

    A theozyme is not defined using Rosetta, but usually in one of the
    following ways: 1) through quantum mechanical calculations (see
    references of Ken Houk's work in the PLoS ONE paper), 2) from chemical 
    intuition, or 3) by stealing or using as inspiration a naturally occurring 
    enzyme’s active site. In this tutorial, we will design a novel triose 
    phosphate isomerase (TIM) based upon a naturally occurring TIM active site. 
    Our theozyme will consist of three catalytic residues and the DHAP ligand 
    from the S. cerevisiae TIM (PDB code 1ney).

4.  Making and checking a .cst file:

    Now that we have dreamed up our theozyme, it needs to be expressed in a 
    format that Rosetta can read. In general, the most unambiguous or precisely 
    defined interaction should come first in the .cst file. We use the Rosetta 
    executable CstfileToTheozymePDB.<extension> to generate .pdb format models 
    from our .cst file so that we can visually check that it defines our 
    theozyme correctly. The command:

        <path_to_mini>/bin/CstfileToTheozymePDB.<extension> -database <path_to_database> -extra_res_fa
        rosetta_inputs/1n1.params -match:geometric_constraint_file rosetta_inputs/mocktim_first_2interactions_example.cst

    produces a file called 
    `PDB_Model_mocktim_first_2interactions_example.cst.pdb` in the working 
    directory which can then be visualized with PyMOL. Inspecting this output 
    .pdb file will ensure that the theozyme geometry that is given to Rosetta 
    is as the user intends. See the additional comments on running 
    CstfileToTheozymePDB in the .cst file documentation linked above.

Step 2: Matching
----------------
(i.e. finding suitable sites for the active site theozyme in a library of scaffold proteins)

1.  Inputs required

    * The .params file from Step 1.

    * The .cst file from Step 1.

    * A library of scaffold protein structures in .pdb format. The scaffold 
      library should be as big as possible. Refer to the matcher documentation 
      linked above for instruction on how to prepare a scaffold for matching 
      (which is mainly deciding where the binding site is and what positions 
      will be considered for theozyme residue attachment). In this tutorial, we 
      will match our theozyme into one scaffold, PDB code 1tml, which can be 
      found here: rosetta_inputs/scaffolds/1tml_11.pdb.

    * A scaffold position file for each scaffold that defines which residues in 
      the scaffold structure will be considered for theozyme residue 
      attachment. The format of scaffold position files, and instructions on 
      how to generate them, can be found in the matcher documentation. In this 
      tutorial, our scaffold position file can be found here: 
      rosetta_inputs/scaffolds/1tml_11.pos.

    * Optionally, a ligand grid file that defines where in three-dimensional 
      space the ligand should be placed during matching. In this tutorial we 
      will not use a ligand grid file, but in case one wants to make sure that 
      the ligand is confined to a certain region of space, it is recommended 
      that one be used. More information on ligand grid files can be found in 
      the matcher documentation linked above.

2.  Outputs generated

    * A .pdb file for each “match”. Each match file contains the 
      three-dimensional coordinates of both the scaffold protein and the 
      theozyme, including the ligand, and will be used as an input in step

    * The number of matches found in the scaffold depends on the complexity of 
      the theozyme, and be anywhere between 0 and hundreds. 

3.  Performing matching

    The command line

        <path_to_mini>/bin/match.<extension> -database <path_to_database> @rosetta_inputs/general_match.flags @rosetta_inputs/1tml_sys.flags

    finds a bunch of matches (~11 in this tutorial) and writes them to the 
    working directory. In a real-life enzdes project, one should look at a few 
    of the matches in PyMOL to make sure that they look roughly as envisioned. 
    Alternatively, the matches could be ranked by similarity to the ideal 
    theozyme (watch for Scott Johnson’s / Ken Houk's EDGE publication).

Step 3: Design
--------------

1.  Inputs required
    * The .params file from step 1.
    * The .cst file from step 1.
    * The matches from step 2.

2.  Outputs generated
    * Designed 'enzymes' in .pdb format.
    * A score file that contains scoring information for each designed enzyme, 
      which will be used in step 4 to evaluate and rank the designs. In this 
      tutorial, our score file generated can be found here: scorefile.txt, 
      while an example score file from a larger design run is provided here: 
      rosetta_inputs/mocktim_all_design_scores.out.

3.  Performing design

    Usually, every match from step 2 is designed several times (between 10-100, 
    depending on computational resources) with the Rosetta enzyme_design 
    executable. For a detailed explanation of what this executable does (and 
    potential options), refer to documentation/paper linked above. Briefly, all 
    residues that are within a given (8A) sphere of the ligand (excepting 
    theozyme residues) are identified and considered changeable in design. 
    These residues are mutated to alanine, and the resulting structure 
    (scaffold with theozyme residues and ligand in a poly-ala cavity) is 
    minimized with respect to the theozyme geometry as specified in the .cst 
    file. Then, 2-3 cycles of constrained sequence design and minimization are 
    carried out, and the final sequence is subjected to an unconstrained repack 
    and minimization. In this tutorial, we will run this stage for one of the 
    matches (rosetta_inputs/UM_1_D41H116K189_1tml_11_mocktim_1.pdb). The 
    command line:

        <path_to_mini>/bin/enzyme_design.<extension> -database <path_to_database> @rosetta_inputs/general_design.flags -s rosetta_inputs/UM_1_D41H116K189_1tml_11_mocktim_1.pdb -out:file:o scorefile.txt

    generates a designed protein in .pdb format and a score file that has one 
    line of values for several score terms and other metrics. In a real life 
    example, there would be a .pdb file for every design as well as a score 
    file that contains one line of values for each design. An example of such a 
    score file with information about all designs is 
    rosetta_inputs/mocktim_all_design_scores.out, which was taken from the runs 
    done for the calculations in the PLoS paper.


Step 4: Evaluating and ranking designs
--------------------------------------

The purpose of this stage is to reduce the number of candidate designs to an 
amount tractable for visual inspection, and get rid of designs that have 
obvious defects. In short, one is looking for designs that have 1) good 
catalytic geometry, 2) good ligand binding score, 3) a preformed/preorganized 
active site and 4) a well behaved protein (expressible, soluble, stably folded, 
etc).

For 1), the constraint energy is taken as a metric, for 2) there is
the straight rosetta ligand binding score, for 3), the designed site
was repacked without the ligand in the design calculation, and for 4)
the designed protein is compared to the scaffold it came from at the
end of the design calculation. For each of these 4 criteria, there are
terms in the scorefile (refer to documentation linked above). The script 
DesignSelect.pl, which is part of the Rosetta3 distribution, can read the 
output scorefile, as well as a file that specifies required values for certain 
columns, and will then output only those designs in the scorefile that satisfy 
all required values:

    perl <rosetta_location>/src/apps/public/enzdes/DesignSelect.pl -d
    rosetta_inputs/mocktim_all_design_scores.out -c
    rosetta_inputs/mocktim_design_selectreqs.txt > selected_designs.txt

These commands will output 44 designs from the 3720 produced for the PLoS ONE 
paper. These 44 would then be visually examined for whether any of them look 
promising enough to be expressed and characterized.
# Ubiquitin Fasta Fragment Generation

The entire workflow for this demo should be described in a file named README.md.
It should describe an entire work flow, with command lines, tested if possible.

The contents of each demo directory should be:

```
starting_files/
  -- directory in which the raw input files are given
rosetta_inputs/
  -- directory in which the modified starting files should
     be placed which will be used as inputs to Rosetta.
scripts/
  -- python scripts, shell scripts, used in the workflow
  -- awk, grep, sed, etc. command lines

README.md
  -- A prose

FOR_AUTHORS.txt
  -- A description for the demo creators of what their demo
     should achieve.
  -- Most of what starts in this file should end up in the
     README file as well.
```
     
## Introduction and Preliminary Installation Steps

The following tutorial will allow you to generate protein fragment files using the `tools/fragment_tools/make_fragments.pl` script.
The directory provides also other necessary files, such as fragment database (vall.jul19.2011) or a script useful to convert between secondary structure prediction file formats (ss_pred_converter.py, also in `tools/fragment_tools/` dir).

**Note, that:**
- All commands listed below assume use of the bash shell.
- It requires the use of third party software that is not included as part of the Rosetta distribution. In particular, you will need to use secondary structure prediction software such as PsiPred, SAM, Jufo, or Porter.

### Instructions for Installation of PsiPred

- Download PsiPred into using the following command:
	`curl –O http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/psipred3.5.tar.gz`

- PsiPred should then be extracted with the following command:
	`tar –xzf psipred3.5.tar.gz`

- You then need to compile PsiPred using the following commands, entered sequentially:
	```
	cd psipred/src/

	make

	make install
	```


If you would like to install other secondary structure prediction software, you will need to obtain and
install the standalone versions of the software yourself. This tutorial does not include information on
the use of these secondary structure prediction programs, but information about them can be found at
the following locations:

**Porter:** To obtain this software, you should send an email to gianluca.pollastri@ucd.ie inquiring about
the standalone version.

**Sam:** Sam can be obtained at: compbio.soe.ucsc.edu/sam2src/

**Jufo:** Information about Jufo can be found here: http://www.meilerlab.org/index.php/bclcommons/
show/b_apps_id/1

**Note:** these may be out of date.

### PsiBlast Installation

You will then need to obtain and install PsiBlast. This can be done using the following commands:

`curl –O ftp://ftp.ncbi.nlm.gov/blast/executables/release/LATEST/blast-*.tar.gz`

It’s important that you download an older version of Blast called “legacy Blast” (not Blast+) because the
current version doesn’t write files in a format Rosetta can understand.

You should then also obtain and install the NCBI non-redundant (nr) database with the following
command:

`curl –O ftp://ftp.ncbi.nih.gov/blast/db/nr0?.tar.gz`

If prompted for a username and password, the username is anonymous and password is any email address.

Within the db dir (which will contain the nr.0x.tar.gz files) run the following command to uncompress the archives: 
`for i in *.tar; do tar –xzf $1; done`

## Preparations to the `make_fragments.pl` script

Once the above has been done, you should open the script in order to modify its contents. The section
that must be modified is labeled “USER CONFIGURATION”, and points the script to the pertinent paths
to the database and programs you have just installed. You should change the paths in this section to
reflect the directories in which you have installed all pertinent programs / databases.

## Actual Fragment Picking
All the steps listed above have to be done only once for a given machine, which now is ready to 
pick fragments. This demo will find fragments from the protein ubiquitin (PDB: 1ubq). The commandline entry for the
script is as follows, and should be run from the directory containing the script.

```
./make_fragments.pl -id 1ubqA -nosam -nojufo -noprof -nonnmake -picker -ref_struc ../starting_files/
1ubq.pdb.gz -nohoms ../starting_files/1ubq.fasta
```

The flags `-nosam`, and `-nojufo`, ensure that the script does not attempt to run SAM or Jufo. The -
picker flag signifies use of the fragment picker. The `-ref_struc` flag signifies the input .pdb file to be
used for analysis, and should include a correct path to the file. You will also need as an input a .fasta
corresponding to the sequence of your protein of interest. The `-nohoms` flag ensures that homologues
to the query sequence are not considered.


Basic Homology Modeling Demo
============================

Conceptual steps
----------------

1. Create fragments for target sequence
2. Find template structure
3. Make a sequence alignment between template sequence and target sequence
4. Download and clean the template pdb
5. Make Rosetta flags file
6. Execute Rosetta
7. Analyze results

Detailed steps
--------------

1.  Create fragments for target sequence

    Obtain a FASTA file of your target sequence. You can get this from 
    [[NCBI|http://www.ncbi.nlm.nih.gov/guide/proteins/]] or you can make the 
    file manually using the following format:

        >template
        AAAAKDLLEKDF

    Open a web browser and go to 
    http://robetta.bakerlab.org/fragmentsubmit.jsp. Paste the target sequence 
    into the sequence input box and click submit. Wait until the server 
    finishes and download all the files produced to your working directory

2.  Find template structure

    Finding a template can be done in multiple ways. One suggestion is to 
    submit the target sequence to HHpred or similar. HHpred can be found here: 
    http://toolkit.tuebingen.mpg.de/hhpred.

3.  Make a sequence alignment between the template sequence and the target 
    sequence

    In order to generate the alignment, you need the primary sequences of both 
    the target and the template. For the template, go to www.pdb.org and search 
    for the PDB code of your template (obtained in Step 2). In this case, we 
    are using 2AST.  In fact, we are using part of chain B from 2AST. On the 
    top right, click the down arrow by "Download Files" and click on "FASTA 
    sequence". Alternatively, you can get the FASTA sequence from NCBI (see 
    Step 1) or make it manually.

    Align the sequences using ClustalW  or the alignment tool of your choice 
    (e.g., EMBOSS, FASTA, Nexus, etc.). Create a FASTA file 
    (template_target.fasta) containing both primary sequences. Go to 
    http://www.ebi.ac.uk/Tools/msa/clustalw3/ and upload the newly created 
    fasta file. The default settings are a good place to start to generate a 
    generally good sequence alignment.  Use the slow alignment for better 
    results Enter your email so you can results. Click submit at the bottom of 
    the page. Depending on the size of the sequences being aligned, it will 
    take anywhere from 5 minutes to some hours or days to run

    After you generate the alignment file, put it in a format like this:

        score 123.456
        t000_			1 VIAFRCPRSFMDQPLAEHFSPFRVQHMDLSNS------VIEVSTL
        2astB_66-105_renumbered	1 ILSLRRSLSYVIQGMANIESLNLSGCYNLTDNGLGHAFVQEIGSL

    where the score is a floating point number, column 1 starting on line 2 is 
    the name of the target (or template, on line 3), column 2 is the beginning 
    sequence position for the threading, and the rest of the line is the 
    sequence.

4.  Download and clean the template PDB

    Once you have decided upon a template, search for the PDB ID at www.pdb.org 
    and go to "Download Files" and select to download the PDB.  This file then 
    needs to be "cleaned" to run properly in Rosetta. To avoid errors when 
    Rosetta reads in the PDB file, the protein must be formatted correctly or 
    “cleaned”. A correctly formatted PDB file includes removed non-ATOM 
    records, renumbered residues from 1, renumbered atoms from 1, and corrected 
    chain ID inconsistencies. The script clean_pdb.py located in the scripts 
    directory will be used to format the template PDB file.

    The clean_pdb.py script requires that you have python2.2 (go to 
    www.python.org for instructions on download and installation) or higher 
    installed, as well as biopython (biopython.org)..  

    Execute the script by typing:

        ./scripts/clean_pdb.py template.pdb A

    The clean_pdb.py script will output two files:

        template_A.pdb
        template_A.fasta

5.  Make the Rosetta flags file

    Rosetta supports a threading protocol (minirosetta application in the bin) 
    which needs a set of input files and specific commandline flags. The 
    commandline flags can be stated in a file in which case the commandline 
    reduces to 

        minirosetta @flags

    The format of the flagsfile is key-value pairs and an example of a working 
    flagsfile to perform thrading looks like:

        -run:protocol threading # run the rosetta threading, loopbuilding, and refinement protocol
        -in:file:alignment ./starting_files/template_target_short.aln # path to alignment file
        -cm:aln_format general # format of alignment file
        -frag3 ./starting_files/fragments/aat000_03_05.200_v1_3.gz # path in 9mer fragment file
        -frag9 ./starting_files/fragments/aat000_09_05.200_v1_3.gz # path to 3mer fragment file
        -in:file:fasta ./starting_files/fragments/t000_.fasta # path to target fasta file
        -in:file:fullatom # we have fullatom format for the input template
        -loops:frag_sizes 9 3 1 # which size fragments are you using?
        -loops:frag_files ./starting_files/fragments/aat000_09_05.200_v1_3.gz ./starting_files/fragments/aat000_03_05.200_v1_3.gz none  # paths to 9mer file, 3mer file, and say "none" for the 1mer file
        -in:file:psipred_ss2 ./starting_files/t000_.psipred_ss2 # path to psipred secondary structure prediction file
        -in:file:fullatom
        -out:nstruct 1 # number of structures you want to build.  Should build at least 1000, but 10,000 would be better
        -in:file:template_pdb ./starting_files/2astB_66-105_renumbered.pdb # path to template pdb
        -database ./rosetta-3.3/rosetta_database/ # path to rosetta database
        -loops:extended # Force extended on loops, independent of loop input file
        -loops:build_initial # Precede loop-modeling with an initial round of just removing the missing densities and building simple loops
        -loops:remodel quick_ccd # closing loops by quick_ccd
        -loops:refine refine_ccd # small movements to remodel loop
        -silent_decoytime # Add time since the last structure was written to score line
        -random_grow_loops_by 4 # randomly grow loops by up to this number of residues on either side.
        -select_best_loop_from 1 # Keep building loops until N and choose best (by score)
        -out:file:fullatom # output in fullatom mode
        -out:output
        -out:file:silent threaded_model.out # silent file stores internal coordinates of the PDB
        -out:file:silent_struct_type binary # output the silent file in binary format
        -out:file:scorefile threaded_model.fasc # output a table of Rosetta scores
        -run:constant_seed # Use a constant seed (1111111 unless specified)
        -run:jran 1111111 # this is good for testing since you should always get the same result
        -overwrite # overwrite any already-existing results having the same name

6.  Execute Rosetta

        ./rosetta-3.3/rosetta_source/bin/minirosetta.macosgccrelease @flags >& test.log &

7.  Analyze results

    There are several options for things to do next, see the following demos:

    * [[analyzing_structure_quality|public/analyzing_structure_quality/readme]]
    * [[clustering|public/clustering/readme]]
    * [[homology_modeling_with_end_extension|public/homology_modeling_with_end_extension/readme]]
    * [[relax_a_large_structure|public/relax_a_large_structure/readme]]

For more information, go to:  
https://wiki.rosettacommons.org/index.php/RosettaCon2011_Documentation_Event
The scripts and input files that accompany this demo can be found in the 
`demos/` directory of the Rosetta weekly releases.
# Symmetry Examples
## Authors
Frank Dimaio and Ingemar André

dimaio@u.washington.edu and ingemar.andre@biochemistry.lu.se

## Brief Description
Examples of how to run symmetry-enabled rosetta protocols.

## Abstract

Symmetric protein assemblies play important roles in many biochemical processes. However, the large size of such systems is challenging for traditional structure modeling methods. This paper describes the implementation of a general framework for modeling arbitrary complex symmetries in Rosetta3.  We describe the various types of symmetries relevant to the study of protein structure that may be modeled using Rosetta’s symmetric framework.  We then describe how this symmetric framework is efficiently implemented within Rosetta, which restricts the conformational search space by sampling only symmetric degrees of freedom, and explicitly simulates only a subset of the interacting monomers.  Finally, we describe structure prediction and design applications that utilize the Rosetta3 symmetric modeling capabilities, and provide a guide to running simulations on symmetric systems.

## Software

Rosetta can be downloaded at www.rosettacommons.org/software

## Documentation

The applications that are exemplified are fully documented in the regular Rosetta documentation. For a general description see http://www.rosettacommons.org/manuals/archive/rosetta3.2_user_guide/.

Documentation for symmetric docking: http://www.rosettacommons.org/manuals/archive/rosetta3.2_user_guide/symmetric_docking.html
Documentation for fold-and-dock: http://www.rosettacommons.org/manuals/archive/rosetta3.2_user_guide/fold_and_dock.html
Documentation for fixed backbone design: http://www.rosettacommons.org/manuals/archive/rosetta3.2_user_guide/fixed_backbone.html
Documentation for comparative modeling: http://www.rosettacommons.org/manuals/archive/rosetta3.2_user_guide/comparative_modeling.html
Backrub Sequence Tolerance
==========================

Author: Colin A. Smith  
Protocol backrub_seqtol

---

This protocol is designed to predict the tolerated sequence space for a
given protein-protein interface or protein domain. It involves generating an
ensemble of backbone structures by running backrub simulations on an input
structure. For each member of the ensemble, a large number of sequences are 
scored and then Boltzmann weighted to generate a position weight matrix for
the specified sequence positions. Interactions within and between different
parts of the structure can be individually reweighted, depending on the
desired objective.

Updates to this protocol capture can be found at:
http://kortemmelab.ucsf.edu/data/

Running the Protocol
--------------------
To run this protocol, the backrub app is used to generate an ensemble of
structures. After that, sequence_tolerance is used to sample a large number
of sequence for each member of the ensemble. Finally, an R script is used to
process the results. A python script is included which handles generation of
single ensemble member using backrub and sequence scoring using
sequence_tolerance. It includes a few paths which must be customized to run
correctly on a user's system. Please note that 20 backbones is the minimum
suggested to get acceptable output. The more backbone structures that are
generated, the less prone the results will be to stochastic fluctuations.

Common flags:

    -s
      This flag specifies the starting structure.
    -resfile
      This is used in backrub and sequence_tolerance to specify mutations and 
      control sequence sampling. It is required for sequence_tolerance.
    -score:weights
      This flag is used to specify a weights file that disables environment 
      dependent hydrogen bonds.
    -score:patch
      This flag must be used to reapply the score12 patch to the standard scoring
      function.
    -ex1 -ex2 -extrachi_cutoff
      These flags enable higher resolution rotamer librares for mutation and
      sequence redesign.

Backrub flags:

    -backrub:ntrials
      This flag is used to increase the number of Monte Carlo steps above the
      default of 1000.
    -backrub:minimize_movemap
      If mutations are specified in the resfile, this movemap is used to 
      specify degrees of freedom to be minimized in a three stage process:
      CHI, CHI+BB, CHI+BB+JUMP.
    -in:file:movemap -sm_prob
      Both of these flags are required to enable small phi/psi moves during
      backrub sampling.

Sequence tolerance flags:

    -ms:checkpoint:prefix -ms:checkpoint:interval
      Both of these flags must be specified to get output of the scored sequences.
    -ms:generations -ms:pop_size -ms:pop_from_ss
      These flags affect the genetic algorithm used for sequence sampling.
    -score:ref_offsets
      This flag is used to reweight the reference energies for given residues.
    -seq_tol:fitness_master_weights
      This flag controls the fitness function used for the genetic algorithm.

Example Rosetta command-line:

    rosetta-3.2/rosetta_source/bin/backrub.linuxgccrelease -database rosetta-3.2/rosetta_database -s input_files/2I0L_A_C_V2006/2I0L_A_C_V2006.pdb -ex1 -ex2 -extrachi_cutoff 0 -mute core.io.pdb.file_data -backrub:ntrials 10000 -score:weights input_files/standard_NO_HB_ENV_DEP.wts -score:patch score12
    rosetta-3.2/rosetta_source/bin/sequence_tolerance.linuxgccrelease -database rosetta-3.2/rosetta_database -s 2I0L_A_C_V2006_0001_low.pdb -ex1 -ex2 -extrachi_cutoff 0 -score:ref_offsets HIS 1.2 -seq_tol:fitness_master_weights 1 1 1 2 -ms:generations 5 -ms:pop_size 2000 -ms:pop_from_ss 1 -ms:checkpoint:prefix 2I0L_A_C_V2006_0001 -ms:checkpoint:interval 200 -ms:checkpoint:gz -score:weights input_files/standard_NO_HB_ENV_DEP.wts -out:prefix 2I0L_A_C_V2006_0001 -score:patch score12 -resfile input_files/2I0L_A_C_V2006/2I0L_A_C_V2006_seqtol.resfile

Using the seqtol_resfile.py python script
-----------------------------------------

The seqtol_resfile.py takes as input a PDB file and generates a resfile for
use with the sequence_tolerance app. It takes at least two other required
arguments. The first is the command used for making residues designable. This
is usually either "ALLAA" for all amino acids, or "PIKAA ..." for a 
restricted set of amino acids. The next arguments are the residues which 
should be designable, with the chain and residue number separated by a
colon.

Example seqtol_resfile.py command-line:

    scripts/seqtol_resfile.py input_files/2I0L_A_C_V2006/2I0L_A_C_V2006.pdb "PIKAA ADEFGHIKLMNPQRSTVWY" B:2002 B:2003 B:2004 B:2005 B:2006

Using the backrub_seqtol.py python script
-----------------------------------------

The backrub_seqtol.py script takes as input a PDB file and other similarly
named configuration files, and produces a single backrub ensemble member
along with approximately 10,000 scored sequences on that member. All of the
input files use a base name derived from removing the ".pdb" extension from
the PDB file. For instance, the base name of 1MGF.pdb would be 1MFG.

If you want to use one PDB file with many different input files you can 
specify a different path from which to get the input files.

Required input files:

    <base name>_seqtol.resfile
      This resfile specifies which sequence positions to sample, along with the
      residue positions that should be repacked.

Optional input files:

    <base name>_backrub.resfile
      This resfile specifies which residues should have flexible side chains
      during the backrub run. By default, all side chains are flexible. This file
      can also define mutations that should be made to the input structure prior
      to the backrub simulation.
    <base name>_minimize.movemap
      This file is passed to the -backrub:minimize_movemap flag (see above).
    <base name>_perturb.movemap
      This file is passed to the -in:file:movemap flag (see above). It also sets
      -sm_prob flag to 0.1.

Example overall command-line:

    scripts/backrub_seqtol.py input_files/2I0L_A_C_V2006/2I0L_A_C_V2006.pdb 1

Post-processing with R:

    $ cd output/2I0L_A_C_V2006
    $ R
    > source("rosetta-3.2/rosetta_source/analysis/apps/sequence_tolerance.R")
    > process_specificity()

Generating Smith & Kortemme PLoS One 2011 figures:

    $ cd output
    $ R
    > rosetta_analysis_dir <- "rosetta-3.2/rosetta_source/analysis"
    > source("../scripts/figures.R")

Versions
--------

This protocol used the Rosetta 3.2 release version

Several other scripting/analysis tools were used:

* Python 2.4.3
* R 2.12.1

References
----------

* Smith, C. A. & Kortemme, T. (2010) Structure-Based Prediction of the Peptide 
  Sequence Space Recognized by Natural and Synthetic PDZ Domains. J Mol Biol 
  402, 460-474.
* Smith, C. A. & Kortemme, T. (2011) Predicting the Tolerated Sequences for 
  Proteins and Protein Interfaces Using RosettaBackrub Flexible Backbone 
  Design. PLoS One 10.1371/journal.pone.0020451
Canonical Sampling for Protein-Protein Docking Refinement
=========================================================

Author: Zhe Zhang (zhezhang1986 at gmail dot com)  
Corresponding PI: Martin Zacharias (martin.zacharias at ph dot tum dot de)  
Last Updated: 06/01/2015  

Reference: Zhang Z, Schindler C, Lange OF, Zacharias M (2015): Comparison of
Replica-Exchange Approaches for Protein-Protein Docking Refinement in Rosetta

---

The different protocols described in this paper are tested on unbound docking
targets selected from docking benchmark4.0. In docking refinement practice,
approximate interaction site is known. Thus the initial start conformation is
generated thus from the superimposed unbound structure, by first translating
the second binding partner 15Å, then rotating 60°. This gives the initial
Ligand RMSD approximately sqrt( 15*15 + (2*pi*r*60/2/360)*(2*pi*r*60/2/360) ),
which r denoting the radii of the second binding partner. The translation
direction and rotation axis are both drawn from uniformly distributed vectors
on unit sphere.

Generating an initial conformation
----------------------------------

### Executable/Script:

    Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease

### PDB:

1. Download cleaned up pdb files from http://zlab.umassmed.edu/benchmark/, or 
   directly from http://www.rcsb.org and clean with clean_pdb.py in folder 
   Rosetta/tools/protein_tools/scripts/
2. Replace chain id if neccessary using scripts/replace_chain.py
       $ cat 1PPE_r_b.pdb 1PPE_l_b.pdb > native.pdb
       $ scripts/replace_chain.py 1PPE_r_u.pdb A > 1PPE_r_u.pdbA.pdb
       $ scripts/replace_chain.py 1PPE_l_u.pdb B > 1PPE_l_u.pdbB.pdb
       $ cat 1PPE_r_u.pdbA.pdb 1PPE_l_u.pdbB.pdb > protAB.pdb

### Running the Application:

Before run the test, please export your rosetta bin and database directory, 
as well as executable of mpirun

    $ export MPI_RUN=$YOUR_MPIRUN_EXECUTABLE
    $ export ROSETTA_BIN=$YOUR_ROSETTA_BIN_DIRECTORY
    $ export ROSETTA_DATABASE=$YOUR_ROSETTA_DATABASE_DIRECTORY

Run the command.sh script provided in folder generate_initial_conformation

    $ ./command.sh

### Example Outputs:

Example outputs are found in the folder generate_initial_conformation.  

1. P.pdb - final decoy, will be used as the initial starting conformation for all
the tests 
2. score.sc - score file

Common rules in this work
-------------------------

For rigid-body docking refinement, we have applied rigid-body mover
UnbiasedRigidBodyPerturbNoCenterMover and sidechain movers including
PerturbChiSidechainMover, PerturbRotamerSidechainMover and
JumpRotamerSidechainMover with the general Metropolis-Hastings framework. The
acceptance of a move is decided by the Metropolis Criterion. In order to avoid
the two binding partners diffuse away from each other in Monte-Carlo move,
 a very loose encounter constraint is applied and acts on the distance of the
mass center of the two binding partners. To achieve local search for docking
refinement, the rigid-body space is restricted with respect to the initial
conformation by translation of 20Å and rotation of 90° in
UnbiasedRigidBodyPerturbNoCenterMover. 

For productive simulation, 2,000,000 Monte-Carlo steps need to be run, and 
snapshots are stored every 1,000 steps. At the end of the simulation, only
decoys generated with the reference setting ( hardrep, temperature=0.15) are
collected and analyzed as final results. To save space and computer time,
the example output in this folder are from much shorter simulation. For 
productive simulations, please also change the related parameter "trial" in
dock.xml file.

Details of each protocol please also refer to the rosetta scripts file
dock.xml in each folder.

All the four protocols need to be run with MPI. 

### Example Outputs

Take the example outputs in wte_remc_docking for example:

1. decoys_P_0001_rt.out - silent trajectory file with rotation matrices and 
   translation vectors started with “RT”, and scores
2. decoys_P_0001_traj.out - silent trajectory file with decoys and scores
3. scores.fsc - silent score file of the trajectory
4. decoys.out - final decoy of a trajectory; this file is a relict of using the 
   JD2-framework and can be generally ignored. 
5. trial.stats - acceptance rate for each mover in each replica
6. tempering.stats - exchange rate between replicas
7. we_bias.grid - well-tempered ensemble bias information at each replica, 
   including grid size, grid range, bias energy in each bin, number of 
   conformations dropped into each bin. Only exist when BiasEnergy is applied, 
   for example in wte_remc_docking and wte_h_remc_docking.
8. decoys_P_0001_m_n.out - checkpoint silent decoy files with m indicating 
   replica number and n indicating the checkpoint number, used for restarting 
   the simulation. When BiasEnergy is applied, in the checkpoint silent file, 
   WTE bias energy information is stored as REMARK started with `REMARK 
   BIASENERGY`

### Analysis

    $ ref=5 # the number of reference replica
    $ scripts/silent_data.py decoys_P_0001_rt.out Lrmsd Irms score I_sc temp_level Fnat_n bias | awk -v ref_rep=“$ref” ’$5==ref_rep’ > collected_data
    $ scripts/collect_tempering_stats.py tempering.stats # collect the average exchange rate over the whole simulation
    $ scripts/collect_trial_stats.py trial.stats # collect average acceptance for each mover at each replica over the whole simulation

MC Docking
----------

This protocol using standard monte-carlo sampling protein-protein docking
refinement. FixedTemperatureController with temperature 0.15 in Rosetta is
used with the MetropolisHastings framework. The magnitude of the step size and
sampling weight of the rigid-body and sidechain movers are fixed along the
entire simulation. Rigid-body mover has a much lower sampling weight than the
sidechain movers to serve the purpose of refinement. 

### Executable/Script

Rosetta/main/source/bin/rosetta_scripts.mpi.linuxgccrelease

### Running the Application

Run the command.sh script provided in folder mc_docking: 

    $ ./command.sh -n $N_PROC   # N_PROC should be (2 + nstruct)

REMC Docking
------------

This protocol using parallel tempering sampling protein-protein docking
refinement. Temperatures are drawn from geometric progression withe the lowest
temperature same as used in mc-docking. HamiltonianExchange is used to control
the temperatures of each replica with the Metropolis-Hastings framework.
Exchange is attempted between neighbor temperatures every 1,000 steps. The
magnitude of the step size and the sampling weight all the movers are
modulated according to the temperature in the initialization such that in the
lower levels more frequent sidechain moves and few small rigid-body moves are
applied and in higher levels less frequent sidechain moves and more bigger
rigid-body moves are applied. 

### Executable/Script

Rosetta/main/source/bin/rosetta_scripts.mpi.linuxgccrelease

### Running the Application

Run the command.sh script provided in folder remc_docking: 

    $ ./command.sh -n $N_PROC   # N_PROC should be (2 + nstruct * n_replica)

WTE REMC Docking
----------------

This protocol applied well-tempered ensemble technique with parallel tempering
to sampling protein-protein docking refinement. By increasing the tunable
factor gamma of BiasedEnergy, we can reduce the number of replicas, but
maintain an approximately same exchange rate between replicas. 

### Executable/Script

    Rosetta/main/source/bin/rosetta_scripts.mpi.linuxgccrelease

### Running the Application

Run the command.sh script provided in folder wte_remc_docking: 

    $ ./command.sh -n $N_PROC   # N_PROC should be (2 + nstruct * n_replica)

WTE H REMC Docking
------------------

In this protocol, two dimensional replica exchange, with variable of the first
dimension as temperature, and of the second dimension as the scaling of
softness of repulsive Lennard-Jones potential. In this second scaling
dimension, we have tested with five levels: standard (hard rep), soft50%,
soft55%, soft60% and soft65%. In the dimension with variable of temperature,
we have five temperatures with lowest equal 0.15. In total 25 levels are run
in parallel and exchange between neighbor levels is attempted every 1,000
steps periodically along the two dimensions.

### Executable/Script

    Rosetta/main/source/bin/rosetta_scripts.mpi.linuxgccrelease

## Running the Application

Run the command.sh script provided in folder wte_remc_docking: 

    $ ./command.sh -n $N_PROC   # N_PROC should be (2 + nstruct * n_replica)

FlexPepDock AbInitio Protocol Capture
=====================================

Written by: Barak Raveh, Nir London, Lior Zimmerman, Ora Schueler-Furman  

---

This is a protocol for the de-novo folding and docking of peptides to proteins (a major extension of the previous FlexPepDock refinement protocol)
For this protocol, no initial information of the peptide backbone is requiered, only a placement of an arbitrary peptide strating 
structure of the peptide within the approximate binding pocket.

Setting up the demo
-------------------

Before running the demo, make sure to change the following variables to you 
local environment:

    FlexPepDockAbInitio/prepack_example:                 set PATH_TO_EXE and PATH_TO_DB
    FlexPepDockAbInitio/run_example:                     set PATH_TO_EXE and PATH_TO_DB
    FlexPepDockAbInitio/scripts/frags/make.sh:set        set PATH_TO_EXE and PATH_TO_DB
    FlexPepDockAbInitio/scripts/clustering/cluster.sh:   set PATH_TO_EXE and PATH_TO_DB 
    FlexPepDockAbInitio/scripts/prep_abinitio.sh:        set pathToDemo and pathToVall
    FlexPepDockAbInitio/scripts/frags/make_fragments.pl: set $scratch and $pathToDemo

Running the demo
----------------

These flags are all found in input_files/flags:

IO flags:

    -s start.ppk.pdb                                        # The start structure of the peptide-protein complex
    -native native.pdb                                      # A reference structure for RMSD calculations - MANDATORY!!!
    -out:pdb_gz                                             # silent output flags
    -out:file:silent_struct_type binary
    -out:file:silent decoys.silent
    -scorefile score.sc                                     # name of scorefile

If using multiple processes and no silent file:

    -multiple_processes_writing_to_one_directory

Number of structures to produce (for demo):

    -nstruct 5                                              # number of structures to produce 

Number of structures to produce (for production run):

    -nstruct 50000

FlexPepDock flags:

    -flexPepDocking:lowres_abinitio
    -flexPepDocking:pep_refine                              # Refine after ab-initio
    -flexPepDocking:flexpep_score_only                      # add aditional interesting scores to scorefile

Packing flags:

    -ex1
    -ex2aro
    -use_input_sc
    -unboundrot native.pdb

Fragment picker flags:

    -frag3 frags/frags.3mers.offset
    -frag9 frags/frags.9mers.offset
    -flexPepDocking:frag5 frags/frags.5mers.offset
    -flexPepDocking:frag5_weight 0.25
    -flexPepDocking:frag9_weight 0.1

Example Rosetta Command Line:

    $PATH_TO_EXE/FlexPepDocking.release -database $PATH_TO_DB @flags

Overall protocol execution (demo):

1.  scripts/prep_abinitio.sh 2b1z (preparation step)

    This will create a 'frags' dir with frgments of the peptide for your run. 
    This will also create links to the native.pdb, start.pdb and flags files 
    from the input_files dir.

2.  prepack_example (prepacking step)

    This will create a start.ppk.pdb structure which is the pre-packed complex 
    to start the simulation from. Also, a ppk.score.sc score file for the 
    repacked structure, as well as a prepack.log log of the run.

3.  run_example (docking step)

    This will create the models of the interaction in a silent (compressed) 
    file (5 for the demo, use 50,000 for real life problems) as well as an 
    initial score file for the models. 

4.  scripts/scoring/rescore.sh score.sc (rescoring step)

    This step is NOT needed if you use Rosetta version 3.3 onwards. If you use 
    version 3.2, use this step to create a new score file (named 'newscore.sc') 
    which includes the fpdock-abinitio recommended score for ranking and 
    clustering of the models (reweighted_sc) as the last column. In this case, 
    use 'newscore.sc' instead of 'score.sc' in the next step (clustering).

5.  scripts/clustering/cluster.sh [ntop] 2 score.sc native.pdb decoys.silent reweighted_sc (clustering step)

    This will cluster the [ntop] lowest energy models (use 5 in this demo, 500 
    for real life problems). The script relies on the specified score file 
    'score.sc' (use 'newscore.sc' if using Rosetta version 3.2), and ranks 
    models according to the specified column 'reweighted_sc' (a reweighted 
    version of the Rosetta score), with the specified clustering radius of 2A. 
    The reprasentative models will be shown in the file 
    'clusters_by_reweighted_sc.txt', sorted by the score of their lowest-energy 
    representatives.

Version
-------
Latest version applies to svn revision 45531 (Oct 2011)


References
----------
Raveh B, London N , Zimmerman L & Schueler-Furman O (2011)
Rosetta FlexPepDock ab-initio: Simultaneous Folding, Docking and Refinement of 
Peptides onto their Receptors PLoS ONE 6(4): e18934.
# Spin Labels
The associated page for this protocol capture can be found at :
https://structbio.vanderbilt.edu:8443/display/MeilerLab/ProtocolCapture

A copy of the page is in MeilerLab-ProtocolCapture.pdf

The input, config, and bin directories have some common materials for the
substeps of the protocol. Each subdirectory has its own input, config, and
bin directories with specific materials for that step. Each subdirectory also
has its own README.txt with specific information about that step. Each step is
entirely self contained, but you might enjoy working through them in the order
given below. 

- [[create_mtssl_mutant/|protocol_capture/spin_labels/create_mtssl_mutant/README]]
- [[relax_mtssl_mutant/|protocol_capture/spin_labels/relax_mtssl_mutant/README]]
- [[relax_mtssl_mutant_membrane/|protocol_capture/spin_labels/relax_mtssl_mutant_membrane/README]]
- [[rotamer_conformation_recovery/|protocol_capture/spin_labels/rotamer_conformation_recovery/README]]
- [[epr_distance_distribution_agreement/|protocol_capture/spin_labels/epr_distance_distribution_agreement/README]]
- [[calculate_cone_model_parameters/|protocol_capture/spin_labels/calculate_cone_model_parameters/README]]
# Using NCAAs in Protein-Peptide Interface Design

## Author(s)
- P. Douglas Renfrew (renfrew@nyu.edu)
- Eun Jung Choi
- Brian Kuhlman

## Reference
[[Incorporation of Noncanonical Amino Acids into Rosetta and Use in Computational Protein-Peptide Interface Design|http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0032637]]

## Brief Description
This is four separate protocol captures that describe the four different ways Rosetta was used in the accompanying publication: creating Rosetta ResidueType parameter files for NCAAs, creating backbone dependent rotamer libraries for NCAAs, calculating explicit unfolded state reference energies for NCAAs, and the running of the DougsDockDesignMinimizeInterface protocol using NCAAs. The protocol caputres for creating ResidueType paramater files, rotamer libraries and explicite unfolded state energies describe the process used in the publication but are written from the standpoint of a researcher looking to add an additional NCAA to Rosetta. 

## Running
See README files for individual files for descriptions about each part was run. 

## Version
The code for this protocol capture is completely checked into trunk rXXXXX and database rYYYYY. 
-------------------------------------------------------------------------------------------------
   INSTRUCTIONS FOR CALCULATING EXPLICIT UNFOLDED STATE ENERGIES FOR NONCANONICAL AMINO ACIDS
-------------------------------------------------------------------------------------------------

Calculating the explicit unfolded state energies is the third of three steps toward being able to use a noncanonical amino acid (NCAA) in Rosetta. To add a new NCAA or to better understand how the NCAAs in the related publication were added one should have already completed or understand the steps in HowToMakeResidueTypeParamFiles and HowToMakeRotamerLibraries. 

The explicit unfolded state energies of an amino acid represent the energy of an amino acid in the unfolded state of a protein and is used to replace the reference energies in Rosetta. The UnfoldedStateEnergyCalculator uses a fragment based method to calculate the average unfolded state energies for each ResidueType. The protocols works on a large set of protein structures that are split in to randomly generated fragments. The central residue of each fragment is mutated to the residue of interest. The fragment is repacked. The unweighted energy for each energy method in the scoring function is recorded for the central residue. After the energies for all fragment central residues are collected, a boltzmann-weighted-average average energy is calculated for each term. 

Calculation explicit unfolded state energies for a NCAA requires three steps:
 - Obtaining a set of input pdbs
 - Running the UnfoldedStateEnergyCalculator protocol on the set of pdbs
 - Modifying the unfolded state energies file in the database

-------------------------------------------
   STEP 01 OBTAINING A SET OF INPUT PDBS
-------------------------------------------

Since the UnfoldedStateEnergyCalculator protocol uses fragments from protein structures, we need a set of high quality structures to work with. Through their PISCES server, the Dunbrack laboratory maintains lists of structures in the Protein Data Bank organized based on xray resolution, precent sequence similarity, and r-factors**. These lists are a convenient way to get a set of high quality structures. In this example we will use a list culled on May 20, 2011. It contains 1801 pdb files that have an an xray resolution of at least 1.6 angstroms, less than 20% sequence identity, and r-factors of less than 0.25. To get the pdbs simply use a supplied script to download the pdbs from the Protein Data Bank ftp servers. 

$ cd inputs
$ ../scripts/get_pdbs.bash cullpdb_pc20_res1.6_R0.25_d110520_chains1859

There should be 1801 gzipped pdb files and a text file containing a list of them called cullpdb_pc20_res1.6_R0.25_d110520_chains1859_list_pruned in the inputs directory. Rosetta will sometimes fail to correctly read in particular pdbs files. The cullpdb_pc20_res1.6_R0.25_d110520_chains1859_list_pruned file is a list of the pdbs which have been screened to be read successfully by Rosetta. 

**Citation: G. Wang and R. L. Dunbrack, Jr. PISCES: a protein sequence culling server. Bioinformatics, 19:1589-1591, 2003. 

----------------------------------------------------------------
   STEP 02 RUNNING THE UNFOLDEDSTATEENERGYCALCULATOR PROTOCOL
----------------------------------------------------------------

The UnfoldedStateEnergyCalculator is relatively easy to run. The command line options are described bellow:

frag_size: single integer value, sets the number of residues in each fragment, should be an odd number and has a default of 5 which is what was used in the accompanying publication
residue_name: string value, sets the three letter code of the residue type which the central residue will be mutated to
repack_fragments: boolean value, controls if the fragments will be repacked before scoring and defaults to true
native_sequence: boolean value, controls if the central residue will be mutated before scoring and defaults to false

Additionally it is strongly recommended to add the following flags as they will make Rosetta handle more pdb files and improves runtime by disabling default features that will be negated by the fragmenting and prepacking

ignore_unrecognized_res: causes Rosetta to ignore unrecognized residue types and 
ex1 and ex2 and extrachi_cutoff 0: force rosetta to use additional rotamer during the fragment repacking
mute all and unmute devel.UnfoldedStateEnergyCalculator and unmute protocols.jd2.PDBJobInputer: reduces the size of the log file significantly by turning off unnecessary output
no_optH true: turns off the hydrogen optimization done when the protein is first read in 
detect_disulf false: turns off disulfide detection

Continuing the ornithine example we have used in the two previous protocol captures, to calculate the unfolded state energies one would run the following command.

$ cd outputs
$ PATH/TO/bin/UnfoldedStateEnergyCalculator.macosgccrelease -database PATH/TO/rosetta_database -ignore_unrecognized_res -ex1 -ex2 -extrachi_cutoff 0 -l ../inputs/cullpdb_pc20_res1.6_R0.25_d110520_chains1859_list_pruned -residue_name C40 -mute all -unmute devel.UnfoldedStateEnergyCalculator -unmute protocols.jd2.PDBJobInputer -no_optH true -detect_disulf false >& ufsec_log_c40 &

NOTE: The extension on your executable my be different.

The run will take between 30-60 seconds per pdb file.

The log file contains lots of useful information. It contains the unweighted energies for each of the energy methods for each of the individual fragments. At the end it will print the average unweighted energies for each ResidueType as well as the Boltzmann weighted average unweighted energies. Boltzmann weighted average unweighted energies are used because some backbones just can't tolerate a mutation to a particular ResidueType and there are extremely high repulsive energies for some fragments that skew the average value. Using the Boltzmann weighting removes the higher energy outliers in a more elegant fashion than a hard energy cutoff.

-----------------------------------------------------
   STEP 03 MODIFY THE UNFOLDED STATE ENERGIES FILE
-----------------------------------------------------

Once the UnfoldedStateEnergyCalculator has finished running the Boltzmann weighted average unweighted energies need to be added to the database. The line you want is the "BOLZMANN UNFOLDED ENERGIES". These are the Boltzmann weighted average unfolded energies for each energy method. The file you need to modify is unfolded_state_residue_energies_mm_std.

Using the ornithine line as an example, the line form the log file is... 

BOLZMANN UNFOLDED ENERGIES:  fa_atr:    -2.462 fa_rep:     1.545 fa_sol:     1.166 mm_lj_intra_rep:     1.933 mm_lj_intra_atr:    -1.997 mm_twist:     2.733 pro_close:     0.009 hbond_sr_bb:    -0.006 hbond_lr_bb:     0.000 hbond_bb_sc:    -0.001 hbond_sc:     0.000 dslf_ss_dst:     0.000 dslf_cs_ang:     0.000 dslf_ss_dih:     0.000 dslf_ca_dih:     0.000

We could add the following to the unfolded_state_residue_energies_mm_std file in the database using the command bellow.

$ echo "C40 -2.462 1.545 1.166 1.933 -1.997 2.733 0.009 -0.006 0.000 -0.001  0.000" >> minirosetta_database/scoring/score_functions/unfolded/unfolded_state_residue_energies_mm_std 

The ResidueType can now be used in almost any Rosetta protocol that is compatible with the MM_STD scoring function.------------------------------------------------------------------------
   INSTRUCTIONS FOR CREATING NONCANONICAL AMINO ACID ROTAMER LIBRARIES
------------------------------------------------------------------------

Creating a Noncanonical Amino Acid (NCAA) rotamer library is the second of two steps toward being able to use a NCAA in Rosetta. To add a new NCAA or to better understand how the NCAAs in the related publication were added one should have already completed or understand the steps in HowToMakeResidueTypeParamFiles.  

Rotamer libraries are sets of common side chain conformations that generally correspond to local minima on the side chain conformational energy landscape. Side chain conformations are usually represented as a set of mean angles and a standard deviation to indicate variability. Rotamer libraries are used in for two main purposes in Rosetta: to provide starting points for side chain optimization routines, and the relative frequency is used as a pseudo-energy. Traditionally rotamer libraries are created by collecting statistics from protein structures. Rosetta uses the backbone dependent Drunbrack rotamer libraries. Since there are not enough structures containing NCAAs they must be generated.

Running the MakeRotLib protocol consists of four steps
 - creating and input template and generating the MakeRotLib options files
 - running the MakeRotLib protocol on each option file
 - Assembling the individual rotamer libraries in a single file
 - modify the ResidueType parameter file to be aware of our new rotamer library


--------------------------------
   STEP 01 MAKING INPUT FILES
--------------------------------

Rosetta primarily uses backbone dependent rotamer libraries. Backbone-dependent rotamer libraries list provide side chain conformations sampled from residue positions whose backbone dihedral angles fall in particular bins. In the case of the Drunbrack rotamer libraries used by Rosetta the bins are in 10 degree intervals for for both phi and psi for a total of 1296 (36*36) phi/psi bins. To replicate this for the NCAAs we need to create a set of side chain rotamers for each member of a set of phi/psi bins.

The MakeRotLib protocol takes an option file as input. It requires an options file for each phi/psi bin. The first step in running it is creating these 1296 options files. Continuing from the HowToMakeResidueTypeParamFiles protocol capture we are again using ornithine as an example. Ornithine has 3 sidechain dihedral angles (chi). We want to sample each chi angle from 0 to 360 degrees in 30 degree intervals, and based on the chemistry of the side chain we predict that were will probably be three preferred angles for each chi angle at 60, 180, and 300 degrees for a total of 27 rotamers (3x3x3). We setup our MakeRotLib options file template as shown bellow.

<<<<< C40_rot_lib_options_XXX_YYY.in start >>>>>
AA_NAME C40
PHI_RANGE XXX XXX 0
PSI_RANGE YYY YYY 0
NUM_CHI 3
CHI_RANGE 1 0  330  30
CHI_RANGE 2 0  330  30
CHI_RANGE 3 0  330  30
CENTROID 300 1 300 1 300 1
CENTROID 300 1 300 1 180 2
CENTROID 300 1 300 1  60 3
CENTROID 300 1 180 2 300 1
CENTROID 300 1 180 2 180 2
CENTROID 300 1 180 2  60 3
CENTROID 300 1  60 3 300 1
CENTROID 300 1  60 3 180 2
CENTROID 300 1  60 3  60 3
CENTROID 180 2 300 1 300 1
CENTROID 180 2 300 1 180 2
CENTROID 180 2 300 1  60 3
CENTROID 180 2 180 2 300 1
CENTROID 180 2 180 2 180 2
CENTROID 180 2 180 2  60 3
CENTROID 180 2  60 3 300 1
CENTROID 180 2  60 3 180 2
CENTROID 180 2  60 3  60 3
CENTROID  60 3 300 1 300 1
CENTROID  60 3 300 1 180 2
CENTROID  60 3 300 1  60 3
CENTROID  60 3 180 2 300 1
CENTROID  60 3 180 2 180 2
CENTROID  60 3 180 2  60 3
CENTROID  60 3  60 3 300 1
CENTROID  60 3  60 3 180 2
CENTROID  60 3  60 3  60 3
<<<<< C40_rot_lib_options_XXX_YYY.in end >>>>>

AA_NAME <three letter code for the amno acid> 
PHI_RANGE <phi value for this bin> <phi value for this bin> 0 : The phi range functionality is not functional. Both values need to be the same and the interval set to 0
PSI_RANGE <psi value for this bin> <psi value for this bin> 0 : The psi range functionality is not functional. Both values need to be the same and the interval set to 0
NUM_CHI <number side chain dihedral angles> : This should be the same as in the parameter file.
CHI_RANGE <chi number> <starting value> <ending value> <interval> : The number of CHI_RANGE fields needs to equal the values specified for NUM_CHI.
CENTROID <Rotamer number for chi 1> <starting value> {<rotamer number for chi 2> <starting value>}{etc.} : CENTROIDS specify the starting points for the K-means clustering described in the related publication. A CENTROID field is needs for each potential rotamer. The number of CENTROID fields defines the number of rotamers listed in the resulting rotamer library.

To generate the 1296 input files we use a provided script that simply replace the XXX and YYY with the phi and psi values. The script is run as shown bellow.

$ cd inputs
$ ../scripts/make_inputs C40_rot_lib_options_XXX_YYY.in

The number of chi angles and the CHI_RANGE sampling interval are the primary determinants of the run time as they determine the number of rotamers that will be tested for each phi/psi bin. It is recommended to have at least 500 samples per chi. In the ornithine example we sample in 30 degree intervals for each of the 3 chi angles giving us a total of 1728 (12x12x12) conformations tested for each phi/psi bin. For a residue with a single chi 1 degree bins will suffice. 

---------------------------------------------
   STEP 02 RUNNING THE MAKEROTLIB PROTOCOL
---------------------------------------------

The next step is to run the MakeRotLib protocol on each of the input files we created in step one. This is the most time consuming portion of the process and should probably be done on a cluster. As cluster setups vary, an example for a single MakeRotLib options file is provided. The other 1295 should be run identically.

$ cd outputs
$ PATH/TO/ROSETTA/bin/make_rot_lib.macosgccrelease -database PATH/TO/rosetta_database -rot_lib_options_file ../inputs/C40_rot_lib_options_-60_-40.in >& C40_rot_lib_options_-60_-40.log &

NOTE: The extension on your executable maybe different.

THe only options passed to the executable are the path to the database and the MakeRotLib options file. After the run completes a file called C40_-60_-40_bbdep.rotlib should be in the output directory. This is the backbone dependent rotamer library for a phi of -60 and a psi of -40.

The log file from the rosetta run in includes quite a bit of useful output.

There are three main sections to the log output,  "ROTAMERS", "CENTROIDS" and "FINAL ROTAMERS" sections. Each one shows the following data: phi, psi, omega and epsilon backbone dihedral angles, probability, total energy, torsion energy, intra-residue repulsive, intra-residue attractive, the number of chi angles, the assigned cluster number, the set of input chi angles, the set of minimized chi angles, the standard deviation, and the distance from that point to each of the cluster centroids.

The log file also displays the number of conformations per cluster, the average distance between the cluster center and the members of that cluster. Lack of conformations in a cluster and a large (>30) average cluster centroid distance suggests that that cluster is higher in energy. 

-----------------------------------------------------------------------------
   STEP 03 ASSEMBLING THE INDIVIDUAL ROTAMER LIBRARIES IN TO A SINGLE FILE
-----------------------------------------------------------------------------

After the MakeRotLib protocol has been run on all of the MakeRotLib options files the individual rotamer libraries for each phi psi bin need to be assembled in to a single file. This is accomplished with a provided script as shown bellow. 

$ cd outputs
$ ../scripts/make_final_from_phi_psi.pl C40

The single file rotamer library should be called C40.rotlib. The file should be placed in the ncaa_rotlibs directory in the database. 

$ cp C40.rotlibs PATH/TO/DATABASE/rosetta_database/ncaa_rotlibs/

--------------------------------------------------
   STEP 04 MODIFYING THE RESIDUE TYPE PARAM FILE
--------------------------------------------------

The last step is modifying the residue type parameter file to use the new rotamer library. To do this we need to add the name of the rotamer library, the number chi angles it describes, and how many bins there are for each chi angle to the Residue type parameter file. The ornithine rotamer library is called C40.rotlib and the rotamer library describes 3 chi angles and each of those 3 chi angle has 3 rotamer numbers. So we would use the following commands to add that information to the file we created in HowToMakeResidueTypeParamFiles.

$ echo "NCAA_ROTLIB_PATH C40.rotlib" >> PATH/TO/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/l-ncaa/ornithine.params
$ echo "NCAA_ROTLIB_NUM_ROTAMER_BINS 3 3 3 3" >> PATH/TO/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/l-ncaa/ornithine.params


-------------------------------------------------------------------
   INSTRUCTIONS FOR RUNNING THE DOUGSDOCKDESIGNMINIMIZE PROTOCOL
-------------------------------------------------------------------

The DougsDockDesignMinimize (DDDM) protocol was used in the accompanying manuscript to redesign the protein/peptide interface of Calpain and a fragment of its inhibitory peptide calpastatin. The protocol was written for this specific protein/peptide interaction and modifications to the code will be necessary to run the protocol on a different system. A modified form was used as the example protocol in the advanced section of the Rosetta 3.0 release manual.

Minor modifications to the protocol have been made from the version used to produce the designs in the accompanying publication. The interfaces of the Rosetta libraries have changed since the initial implementation of the protocol and these modifications were necessary to allow the protocol to work with the current release. 

The protocol has two loops referred to as the inner and outer loops. The outer loop controls the number of structures generated, outputting structures only if certain filters are met. The inner loop iterates between two phases: a perturbation phase and a design phase. During the perturbation phase three types of perturbations are used. Two types of perturbations are applied to the peptide: rotational and translational rigid body permutation of the peptide in the binding pocket of the protein, and small and shear perturbations of the backbone dihedral angles in all residues of the peptide. The third type of perturbation is preformed on residues 1-45 of the N-terminus which comprise the peptide binding site and surrounding residues of the protein are perturbed using small and shear perturbations of the backbone dihedral angles. One of the three perturbations is randomly chosen to be used each time the function in called and is followed by a coarse side chain optimization (RotamerTrials). During the design phase of the protocol the side chains positions are optimized using a more intensive side chain optimization routine (PackRotamers) followed by minimization of the backbone and side chain dihedrals angles as well as the jump between the peptide and protein. The number of perturbations per design as well as the magnitude of perturbations are controlled using command line options described bellow.

There are three basic steps to running the protocol as it was used in the accompanying manuscript.
 - generating the input resfiles and folders
 - running the DougsDocDesignMinimize protocol for each mutation at each position 
 - running the analysis scripts

Example output is provided for a small version of the full run.

------------------------------------------
   GENERATING INPUT RESFILES AND FOLDERS
------------------------------------------
The Protein Databank code for the Calpain/Calpastatin structure used in the designs is 1NX1. The files contain two copies of calpain (chains A and B) and two copies of the calpastatin inhibitory peptide (chains C and D). To reduce the size of the system, designs were done on chain A and chain C only because chains A and C have lower B-factors than chains B and D. The calpain/calpastatin interface is distal to the calpain/calpain interface and the resides that make up the calpain/calpain interface were held fix during the protocol. Additionally the calcium atoms and water molecules were removed. Rosetta doesn't work with with water molecules and the calcium ions are not near the protein peptide interface. This input pdb was repacked and both the side chain and backbone dihedrals minimized using a modified version of the fixbb app. The input PDB file is called 1NX1_clean_repack_min_all.pdb

Each of the NCAAs added in the accompanying publication was tried at each position in the peptide. To do this and to keep all of the output pdb organized a script is provided that creates folders and generates resfiles based on templates for each sequence. To generate the resfiles and folders run the following commands...

$ cd run_dir
$ ../scripts/make_folders_resfiles.tcsh

NOTE: To save space, the script has been modified to only produce resfiles and folders for position 610 in the peptide and residue type MPA (4-methyl-phenylalanine). To modify the script to produce folders and resfiles for each position simply uncomment the lines that read "#foreach i ( 601 602 603 604 605 606 607 608 609 610 611 )" and "#foreach j ( ABA APA HLU ..... C92 C93 C94 )" comment out the line that reads "foreach i ( 610 )" and "foreach j ( MPA )".

Using the templates in the run_dir the make_folders_resfiles.tcsh script makes folders and resfiles to preform all of the DDDMI runs.


-------------------------------------------------------
   RUNNING DOUGSDOCKDESIGNMINIMIZEINTERFACE PROTOCOL
-------------------------------------------------------

To run the protocol modifications need to made to files in the database. In the file rosetta_database/chemical/residue_type_sets/fa_standard/residue_types.txt all of the paths to the residue type parameter files under the L-NCAA heading need to be uncommented by removing the "#" from the front of the line. Additionally the rotamer libraries for the NCAA are not provided in the default Rosetta database because they are more than 400MB. The rotamer libraries for the NCAAs added in the accompanying publication are provided as supplemental information. 

NOTE: Turning on all of the additional residue types dramatically increases the number of residue types and the memory footprint of Rosetta. The memory foot print can be reduced by commenting out unnecessary patches in the rosetta_database/chemical/residue_type_sets/fa_standard/patches.txt file. For the DougsDockDesignMinimizeProtocol all but the NtermProteinFull.txt and CtermProteinFull.txt can be safely commented out by placing "#" symbols at the beginning of each line of the patches.txt except for the lines that say "NtermProteinFull.txt" and "CtermProteinFull.txt".

To run the protocol as in the accompanying publication preform the following commands starting at the HowToRunDougsDockDesignMinimizeProtocol directory.
$ cd run_dir
$ cd pos_610_MPA
$ /PATH/TO/ROSETTA/bin/doug_dock_design_min_mod2_cal_cal.macosgccrelease -database /PATH/TO/rosetta_database -s ../../inputs/1NX1_clean_repack_min_all.pdb -resfile ../resfile_pos_603_MPA -nstruct 255 -inner_num 45 -pert_num 25 -ia_ener 100 -use_input_sc -pdb_gz

The above command generates 255 structures and will take approximately 5 minutes per structure depending on your hardware.

NOTE: The extension of your executable maybe different than the above. Also in the publication the nstruct command line option was 255. However to save space for the protocol capture an nstruct of 10 was used.

A script is provided that will preform the above command for each folder created by the make_folders_resfiles.tcsh script.

$ cd run_dir
$ ../scripts/run_script.bash

NOTE: You will need to set the path to your database and executable in the run_script.bash.

There are a number command line options to control the running of the application.

pert_mc_temp: The MC temperature to use for the perturbation phase of the DDDM protocol, defaults to 0.8 kT.
pert_dock_rot_mag: The rotation magnitude for the ridged body perturbation in the perturbation phase of the DDDM protocol, defaults to 0.5. 
pert_dock_trans_mag: The translation magnitude for the ridged body perturbation in the perturbation phase of the DDDM protocol, defaults to 0.25. 
pert_pep_small_temp: The temperature for the internal MC object in the small mover of the peptide perturbations, defaults to 0.8 kT. 
pert_pep_shear_temp: The temperature for the internal MC object in the shear mover of the peptide perturbations, defaults to 0.8 kT. 
pert_ter_small_temp: The temperature for the internal MC object in the small mover of the termini perturbations, defaults to 0.8 kT. 
pert_ter_shear_temp: The temperature for the internal MC object in the shear mover of the termini perturbations, defaults to 0.8 kT. 
pert_pep_small_H: The maximum angle of perturbation for helical secondary structure for the peptide small mover, defaults to 1 degree.
pert_pep_small_L: The maximum angle of perturbation for loop secondary structure for the peptide small mover, defaults to 1 degree.
pert_pep_small_E: The maximum angle of perturbation for strand secondary structure for the peptide small mover, defaults to 1 degree.
pert_pep_shear_H: The maximum angle of perturbation for helical secondary structure for the peptide shear mover, defaults to 1 degree
pert_pep_shear_L: The maximum angle of perturbation for loop secondary structure for the peptide shear mover, defaults to 1 degree.  
pert_pep_shear_E: The maximum angle of perturbation for strand secondary structure for the peptide shear mover, defaults to 1 degree.
pert_ter_small_H: The maximum angle of perturbation for helical secondary structure for the termini small mover, defaults to 1 degree
pert_ter_small_L: The maximum angle of perturbation for loop secondary structure for the termini small mover, defaults to 1 degree.  
pert_ter_small_E: The maximum angle of perturbation for strand secondary structure for the termini small mover, defaults to 1 degree.
pert_ter_shear_H: The maximum angle of perturbation for helical secondary structure for the termini shear mover, defaults to 1 degree
pert_ter_shear_L: The maximum angle of perturbation for loop secondary structure for the termini shear mover, defaults to 1 degree.  
pert_ter_shear_E: The maximum angle of perturbation for strand secondary structure for the termini shear mover, defaults to 1 degree.
pert_pep_num_rep: Number of small and shear iterations for the peptide, defaults to 100. 
pert_ter_num_rep: Number of small and shear iterations for the terminus, defaults to 100. 
pert_num: Number of iterations of perturbation loop per design, defaults to 100.
inner_num: Number of iterations of the inner loop, defaults to 100. 
ia_ener: Upper energy limit for final design/interface analysis checkpoint, defaults to 0.0. 
desn_mc_temp: The temperature to use for the design/minimization phase of the DDDM protocol, defaults to 0.8 kT.

For the most part the defaults should suffice and are what was used in the paper. The "ia_ener" command line option is dependent on the complex that the protocol is run on. Setting it unreasonably low will cause the protocol to run forever and never output any structures. Setting it to 100 above allows all but the most aggressions structures to be output. The nstruct command line option controls the outer loop. 

-----------------------------
   ANALYZING THE RESULTS
-----------------------------

Each of the output pdb files contains information about the protein peptide complex that can used to evaluate the designs.  For example at the end of the 1NX1_clean_repack_min_all_0004.pdb is shown bellow. For each filter the value is calculated for the protein and peptide together (COMPLEX) and separated by 1000 angstroms (SEPARATE) and the difference between the two (DIFF). ENERGY is the Rosetta energy. The ENERGY_COMPLEX is the primary determinant to how good a design is and the ENERGY_DIFF can give an estimate for the binding energy. SASA is the solvent accessible surface area. SASA_DIFF is indicative of sequences that make a more protein-peptide contacts and can be used for screening designs for example placing a very large side chain at a constrained interface position can cause the peptide to be pushed out of the binding pocket which would be reflected in a smaller magnitude SASA_DIFF. HB_ENER is the hydrogen bonding component of the Rosetta energy. Larger HB_ENER_DIFF values indicate that the design is making more or better hydrogen bonds across the protein peptide interface. PACK is the RosettaHoles score and is a measurement of how well the protein is packed. A PACK_COMPLEX that is larger than the PACK_SEPARATE is favorable and suggests that the complex is better packed than the protein alone. Additionally the RosettaHoles score penalizes holes that cannot be occupied by solvent so larger PACK_DIFF score indicate that the designed peptide is capable of filling cavities in the protein that are inaccessible to solvent.

ENERGY_COMPLEX:	   -63.3501
ENERGY_SEPERATE:   -48.599
ENERGY_DIFF:	   -14.7512
SASA_COMPLEX:	   11774.1
SASA_SEPERATE:	   13092.3
SASA_DIFF:	   -1318.15
HB_ENER_COMPLEX:   -146.728
HB_ENER_SEPERATE:  -145.037
HB_ENER_DIFF:	   -1.69085
PACK_COMPLEX:	   0.515277
PACK_SEPERATE:	   0.466995
PACK_DIFF:	   0.0482824

A script is provided that pulls the information out of a set of pdb files and sorts it based on the ENERGY_COMPLEX metric. 

$ cd run_dir
$ cd pos_610_MPA
$ ../../scripts/get_interface_data.tcsh

The script produces a file call out.ALL that contains a single line for each pdb file with the metrics in the above order.------------------------------------------------------------------------------------------------
   CREATING ROSETTA RESIDUETYPE PARAMETER FILES FOR NONCANONICAL AMINO ACIDS
------------------------------------------------------------------------------------------------

Creating a Residue Type Parameter file is the first of three steps to being able to use a noncanonical amino acid (NCAA) in Rosetta. To add a new NCAA or to better understand the how the NCAAs in the accompanying publication were added one should follow the steps in this folder followed by the steps in HowToMakeNCAARotamerLibraries.  

ResidueType Parameter files are the basic input file for the creation the core::chemical::ResidueType in Rosetta. ResidueTypes are used for storing information about the atoms, bonds, and chemical connectivity of both polymer (DNA/RNA, peptide) and ligand residues in Rosetta. Including an idealized conformation in an internal coordinate format. They also contain parameters used by some of the energy methods, and rotamer libraries.

This directory includes modified open babble source, molfile2params_polymer.py script, instructions, inputs and outputs for creating a param file for a alpha-peptide backbone with a noncanonical side chain. The example alpha-peptide is ornithine.

There are five basics steps to generating a ResidueType parameter file as in the accompanying publication outlined bellow:
- Creating an initial structure
- Minimizing that structure
- Converting it to a hybrid molfile 2000/3000 format
- Modifying the molfile with additional information used by the molfile_2_params_polymer.py script
- Running the molfile_2_params_polymer.py script

---------------------------------
   STEP 00 BUILDING OPENBABBLE
---------------------------------
Openbabble is an open source program distributed under the GLPv2 license used to create a specially formated molfile to be used as input for the molfile_2_params_plymer.py script. The copy included here has been modified to produce molfile 2000 format files with the molfile 3000 bond types. 

You can build the modified open babble source using the standard configure, make, make install. I like to install it locally but you can install it anywhere.

Note: The path for the configure script will need to change on your system.

$ cd openbable
$ mkdir install
$ tar xvzf openbabel-2.2.0.tar.gz
$ cd openbabel-2.2.0
$ ./configure --prefix="/PATH/TO/HowToMakeResidueTypeParamFiles/openbabel/install/" 
$ make
$ make install

the executable babble should now be in /PATH/TO/HowToMakeResidueTypeParamFiles/openbabel/install/bin/

------------------------------------------
   STEP 01 GENERATING INITIAL STRUCTURE
------------------------------------------
The first step is to generate an initial structure. I usually do this in PyMOL but any molecular editor will work. PyMOL has rudimentary but very functional editing capabilities. The goal is to get the structure close to the final structure. Don't stress too much because the minimization will clean up most things. Make sure you double check chirality. To make sure the geometry for the atoms near the ends of the molecule is correct you will need to put capping groups on the molecule. I would recommend an acetyl (ACE) on the nitro-terminus of and a n-methyl (NME) on the carboxy terminus. This structure (ACE-X-NME) is usually called a dipeptide despite the fact that it only has one complete residue.

PRO TIPS:
 - Within PyMOL go to the "Mouse" drop down and select "3 button editing".
 - Within PyMOL type "help edit_keys" in the command box will bring up the instructions for the editor. 
 - Within PyMOL type "set valence, 1" in the command box to show single, double, triple bonds.
 - If you are making multiple ResidueTypes, be consistent with atom order. It makes future steps easier. I like to put all the atoms for the capping groups first, followed by backbone heavy atoms, side chain heavy atoms, and then then hydrogens.

EXAMPLE:
See the example pdb files in the folder stage_01_initial_structures and shown bellow...
<<<<< ornithine.pdb start >>>>>
ATOM      1  C   ACE     1       0.984   0.045  -0.578  1.00  0.00           C
ATOM      2  O   ACE     1       1.815  -0.721  -1.083  1.00  0.00           O
ATOM      3  CH3 ACE     1      -0.445  -0.396  -0.349  1.00  0.00           C
ATOM      4 1HH3 ACE     1      -0.468  -1.251   0.313  1.00  0.00           H
ATOM      5 2HH3 ACE     1      -1.013   0.407   0.098  1.00  0.00           H
ATOM      6 3HH3 ACE     1      -0.904  -0.668  -1.288  1.00  0.00           H
ATOM     26  N   NME     3       4.937   2.266   0.684  1.00  0.00           N
ATOM     27  CH3 NME     3       5.341   3.324   1.592  1.00  0.00           C
ATOM     28  H   NME     3       5.534   1.511   0.316  1.00  0.00           H
ATOM     29 1HH3 NME     3       4.689   4.177   1.478  1.00  0.00           H
ATOM     30 2HH3 NME     3       5.286   2.977   2.614  1.00  0.00           H
ATOM     31 3HH3 NME     3       6.355   3.627   1.379  1.00  0.00           H
ATOM      7  N   ALA     2       1.392   1.405  -0.195  1.00  0.00           N
ATOM      8  CA  ALA     2       2.806   1.518  -0.539  1.00  0.00           C
ATOM      9  C   ALA     2       3.512   2.478   0.390  1.00  0.00           C
ATOM     10  O   ALA     2       2.929   3.435   0.912  1.00  0.00           O
ATOM     11  CB  ALA     2       2.895   1.945  -2.014  1.00  0.00           C
ATOM     12  C01 ALA     2       2.214   0.926  -2.946  1.00  0.00           C
ATOM     13  N01 ALA     2       1.680   0.434  -5.295  1.00  0.00           N
ATOM     14  C02 ALA     2       2.334   1.404  -4.405  1.00  0.00           C
ATOM     15  H   ALA     2       0.770   2.167   0.251  1.00  0.00           H
ATOM     16  HA  ALA     2       3.283   0.529  -0.413  1.00  0.00           H
ATOM     17 2HB  ALA     2       2.409   2.923  -2.195  1.00  0.00           H
ATOM     18 3HB  ALA     2       3.943   2.035  -2.355  1.00  0.00           H
ATOM     19  H01 ALA     2       2.698  -0.045  -2.840  1.00  0.00           H
ATOM     20  H02 ALA     2       1.461  -0.405  -4.777  1.00  0.00           H
ATOM     21  H03 ALA     2       1.161   0.837  -2.679  1.00  0.00           H
ATOM     22  H04 ALA     2       2.300   0.206  -6.059  1.00  0.00           H
ATOM     23  H05 ALA     2       3.387   1.492  -4.674  1.00  0.00           H
ATOM     24  H06 ALA     2       1.850   2.375  -4.509  1.00  0.00           H
ATOM     25  H07 ALA     2       0.828   0.835  -5.661  1.00  0.00           H
END
<<<<< ornithine.pdb end >>>>>

-----------------------------------
   STEP 02 MAKING GAUSSIAN INPUT
-----------------------------------
Rosetta assumes that We need to minimize the initial structure we made in PyMOL to get a good set of ideal bond lengths and angles. We will use Gaussian, an is a ab initio quantum mechanics package from Schodinger Inc., to do this but You can use the molecular modeling program of your choice. To do this we will take the coordinates from the pdb file and put them into a gaussian input file. Like the one shown bellow and in stage_02_gaussian_input. A discussion of the complete gaussian input structure is beyond the scope of this document. Documentation for gaussian can be found here http://www.gaussian.com/ . 

In short the input is as follows...
line 1 sets the path for the checkpoint file
line 2 describes the level of theory, options, and convergence criteria. You may need to change the basis set to something smaller if you are using stuff bellow the 4th line of the periodic table. 
lines 3-5 are comments
line 6 is the charge and multiplicity (this is usually 0 1 but ornithine is charged)
line 7-37 are the elemental type and xyz coordinates of the atoms from the pdb file in the same order as the pdb file
line 38 is blank
line 39-40 is the modredundant input that says that we want to keep the torsion formed by atoms 1 13 14 and 15 fixed at 150.00 degrees
line 41 is blank.

Running Gaussian is simple but the minimizations take a long time (a few hours per structure). The command bellow will run gaussian on the input file in the stage_02 folder and put the output in the stage_03 folder.

$ g03 stage_02_gaussian_input/ornithine.com stage_03_gaussian_output/ornithine.log

<<<<< ornithine.com start >>>>>
%Chk=stage02_ace_nme_res_ordered_pdbs/ornithine.chk
# HF/6-31G(d) Opt=ModRedundant SCF=Tight Test 

scan rotamers
 
1  1
C        0.984   0.045  -0.578
O        1.815  -0.721  -1.083
C       -0.445  -0.396  -0.349
H       -0.468  -1.251   0.313
H       -1.013   0.407   0.098
H       -0.904  -0.668  -1.288
N        4.937   2.266   0.684
C        5.341   3.324   1.592
H        5.534   1.511   0.316
H        4.689   4.177   1.478
H        5.286   2.977   2.614
H        6.355   3.627   1.379
N        1.392   1.405  -0.195
C        2.806   1.518  -0.539
C        3.512   2.478   0.390
O        2.929   3.435   0.912
C        2.895   1.945  -2.014
C        2.214   0.926  -2.946
N        1.680   0.434  -5.295
C        2.334   1.404  -4.405
H        0.770   2.167   0.251
H        3.283   0.529  -0.413
H        2.409   2.923  -2.195
H        3.943   2.035  -2.355
H        2.698  -0.045  -2.840
H        1.461  -0.405  -4.777
H        1.161   0.837  -2.679
H        2.300   0.206  -6.059
H        3.387   1.492  -4.674
H        1.850   2.375  -4.509
H        0.828   0.835  -5.661

1 13 14 15  -150.00 F
13 14 15 7  150.00 F

<<<<< ornithine.com end >>>>>

------------------------------------------------
   STEP 03 CONVERT GAUSSIAN OUTPUT TO MOLFILE
------------------------------------------------
The next step is to produce a hybrid molfile v2000/3000 file. The program babble we built in the first step can convert the gaussian output to a molfile that has the v2000 structure but uses the v3000 bond types. In particular the molfile_2_params_polymer.py script needs to have aromatic bonds types (type 4) for aromatic rings instead of alternating single and double bonds (Kekule structure). 

The command bellow will convert the gaussian output to molfile format for the example...

$ openbabel/install/bin/babel -i g03 stage_03_gaussian_output/ornithine.log -o mol stage_04_molfile/ornithine.mol

<<<<< ornithine.mol start >>>>>
 OpenBabel01101117083D

 31 30  0  0  0  0  0  0  0  0999 V2000
    0.4866    2.2296   -0.2354 C   0  0  0  0  0
    1.2081    1.9315   -1.1549 O   0  0  0  0  0
    0.4369    3.6268    0.3352 C   0  0  0  0  0
   -0.3187    4.1933   -0.1996 H   0  0  0  0  0
    0.1852    3.6362    1.3888 H   0  0  0  0  0
    1.3920    4.1075    0.1799 H   0  0  0  0  0
   -2.6920   -1.1572   -0.6650 N   0  0  0  0  0
   -4.0550   -1.5847   -0.3922 C   0  0  0  0  0
   -2.3500   -1.2554   -1.5934 H   0  0  0  0  0
   -4.1238   -1.9440    0.6231 H   0  0  0  0  0
   -4.7601   -0.7727   -0.5220 H   0  0  0  0  0
   -4.3066   -2.3874   -1.0709 H   0  0  0  0  0
   -0.3415    1.3293    0.3431 N   0  0  0  0  0
   -0.6023    0.0286   -0.2270 C   0  0  0  0  0
   -2.0235   -0.3660    0.1877 C   0  0  0  0  0
   -2.4547   -0.0133    1.2526 O   0  0  0  0  0
    0.3558   -1.0672    0.2829 C   0  0  0  0  0
    1.8184   -0.8036   -0.0848 C   0  0  0  0  0
    4.1559   -1.6146    0.0004 N   0  0  0  0  0
    2.7190   -1.9434    0.3678 C   0  0  0  0  0
   -0.9580    1.6065    1.0765 H   0  0  0  0  0
   -0.5197    0.1005   -1.3044 H   0  0  0  0  0
    0.2486   -1.1369    1.3608 H   0  0  0  0  0
    0.0400   -2.0197   -0.1352 H   0  0  0  0  0
    1.9103   -0.6656   -1.1566 H   0  0  0  0  0
    4.2542   -1.4797   -0.9970 H   0  0  0  0  0
    2.1487    0.1188    0.3755 H   0  0  0  0  0
    4.7999   -2.3436    0.2758 H   0  0  0  0  0
    2.4931   -2.8811   -0.1191 H   0  0  0  0  0
    2.7087   -2.0883    1.4385 H   0  0  0  0  0
    4.4568   -0.7571    0.4437 H   0  0  0  0  0
  2  1  2  0  0  0
  1  3  1  0  0  0
  1 13  1  0  0  0
  4  3  1  0  0  0
  6  3  1  0  0  0
  3  5  1  0  0  0
  9  7  1  0  0  0
  7  8  1  0  0  0
  7 15  1  0  0  0
 12  8  1  0  0  0
 11  8  1  0  0  0
  8 10  1  0  0  0
 14 13  1  0  0  0
 13 21  1  0  0  0
 22 14  1  0  0  0
 14 15  1  0  0  0
 14 17  1  0  0  0
 15 16  2  0  0  0
 24 17  1  0  0  0
 18 17  1  0  0  0
 17 23  1  0  0  0
 25 18  1  0  0  0
 18 20  1  0  0  0
 18 27  1  0  0  0
 26 19  1  0  0  0
 19 28  1  0  0  0
 19 20  1  0  0  0
 19 31  1  0  0  0
 29 20  1  0  0  0
 20 30  1  0  0  0
M  END
<<<<< ornithine.mol end >>>>>

-----------------------------------
   STEP 04 MODIFYING THE MOLFILES
-----------------------------------
The molfile2params_polymer.py script requires some additional data to be added to the end of the molfile. This data is specified at the end of the file after the bond information. It is a list of variable names and then a list of values. The variable are described bellow. 

ROOT: Single numerical value. Atom number (according to the order in the molfile). Where the atom tree is rooted for this residue type. Should be the nitrogen of the central residue.
POLY_N_BB, POLY_CA_BB, POLY_C_BB, POLY_CO_BB: Atom number (according to the order in the molfile). The backbone nitrogen, alpha-carbon, carbonyl-carbon and carbonyl-oxygen. These get special rosetta atom types and so are listed here special.
POLY_IGNORE: List of atom numbers (according to the order in the molfile). These are the atoms for the capping groups with the exception of the upper and lower connect atoms. They will not be listed in the atoms and bonds in the params file but are used in determining the atom types.
POLY_UPPER, POLY_LOWER: Atom number (according to the order in the molfile). These are the atoms in the capping groups that connect to the residue. They will not be listed in the atoms and bonds in the params file but are used in determining the atom types and they are listed in the internal coordinate section.
POLY_CHG: Single numerical value. Overall charge on the residue. 
POLY_PROPERTIES: List of alpha-numerical values. These get used by Rosetta at various places in the program. You can say something like "if ( pose.residue(10).type().is_protein() ) { // do something }".
END: The end of the file.

Note: There are 2 spaces between the "M" and the variable name.
Note: If you have to make multiple residue types keeping the atoms in the same order makes assigning all these numbers easier. 

Additional info for peptide ornithine...
M  ROOT 13
M  POLY_N_BB 13
M  POLY_CA_BB 14
M  POLY_C_BB 15
M  POLY_O_BB 16
M  POLY_IGNORE 2 3 4 5 6 8 9 10 11 12
M  POLY_UPPER 7
M  POLY_LOWER 1
M  POLY_CHG 1
M  POLY_PROPERTIES PROTEIN POLAR CHARGED
M  END

---------------------------------------
   STEP 05 RUNNING MOLFILE2PARAMS.PY
---------------------------------------
Finally, the molfile2params_polymer.py script will convert the modified molfile to a params file. The commands bellow work for the example files and produce a ResidueType parameter file for ornithine with the three letter code C40. 

$python scripts/molfile_to_params_polymer.py --clobber --polymer --no-pdb --name C40 -k ornithine.kin stage_05_modified_molfile/ornithine.mol

There may need additional tweaking that needs to happen to the params files to make them work correctly. Compare these to the ones in the database for reference. In particular it is important to double check that the chi angles and the atoms that comprise each of them are specified correctly. Additionally it is important to check that the dependancies in the internal coordinates are correctly setup with respect the the atoms specified for the chi angles. For example, in ornithine the chi1 is defined as

CHI 1  N    CA   CB   CG 

and the internal coordinates for it and it hydrogens are

ICOOR_INTERNAL    CG   -62.608549   66.937720    1.530976   CB    CA    N  
ICOOR_INTERNAL   1HG   121.953676   70.110328    1.084548   CG    CB    CD 
ICOOR_INTERNAL   2HG   116.869954   70.361318    1.082495   CG    CB   1HG 

Since the gamma hydrogens are defined relative to the CG they will move when the chi1 is rotated.

Pro Tip: The molfile2params.py script can produce a kinemage file. This is handy as it lets you check how rosetta will build the atom tree, and the atom type assignments, and other parameters. You can open kinemage files using the KiNG program from the Richardson lab at Duke (http://kinemage.biochem.duke.edu/software/king.php).

----------------------------------------------------------------
   STEP 06 ADDING THE RESIDUE TYPE PARAM FILE TO THE DATABASE
----------------------------------------------------------------

The last step is to add the new param file to the rosetta database and make Rosetta aware that it is there. Simply copy the param file to the appropriate directory for the case of ornithine it is "l-ncaa"

$ cp C40.parms minirosetta_database/chemical/residue_type_sets/fa_standard/residue_types/l-ncaa/ornithine.params

Next add a line to the end of the residue_types.txt file.

$ echo "l-ncaa/ornithine.params" >> minirosetta_database/chemical/residue_type_sets/fa_standard/residue_types.txtsrc/python/apps/curated/
	A directory with published protocols and curated applications.
	Must contain tutorial on the use of the protocol.
	Must describe the scientific benchmark performed to validate the
	protocol.
	Must report the svn revision number used in the benchmark.

src/python/apps/public/
	A directory with protocols people are playing around with and have
	found useful but aren't quite ready to move into the curated section.
	Must contain tutorial on the purpose and use of the protocol.

src/python/apps/pilot/
	An area where people can put anything that they are working on but
	haven't completely finished yet.
# Stepwise Enumerative Assembly
## Author
Rhiju Das, rhiju@stanford.edu

## Protocol Name
Was called: "StepWiseAssembly", but that was too generic. 
We need a catchier name. SEA? Enumerosetta?

## Brief Description

An enumerative ansatz for high resolution RNA and Protein Folding
Rhiju Das, Parin Sripakdeevong, Das Lab, Stanford Biochemistry
Tuesday Aug. 3, 2010

## Abstract

High-resolution structure modeling is severely limited by difficulties in conformational sampling. Current Rosetta approaches use a low-resolution sampling stage, fragments from experimental structures, or a Monte-Carlo-like search -- typically all three. We describe an alternative "ansatz" that builds well-packed models in small steps, enumerating several million conformations for each residue, and covering all possible paths. We present results on non-canonical RNA motifs as well as highly irregular protein loops that have been intractable for prior fragment assembly or analytic loop closure approaches. In all cases, the method either reaches atomic accuracy or exposes flaws in Rosetta’s high-resolution energy function. Blind tests on a tetraloop-receptor motif are underway, as well as extension of the method to more complex systems such as aptamers and knotted cyclotides.


## Running
### Example Rosetta Command Line:

```
  stepwise_protein_test.macosgccrelease  -rebuild -fasta mini_1alc.fasta  -cluster:radius 0.1  -score:weights score12_no_hb_env_dep.wts  -pack_weights pack_no_hb_env_dep.wts -add_peptide_plane -align_pdb mini_1alc_H.pdb -native mini_1alc_H.pdb  -s1 noloop_mini_1alc_H.pdb  -input_res1 `seq 1 11` `seq 20 28` -sample_res 12 -out:file:silent build_first_residue.out  -calc_rms_res `seq 12 19` -fixed_res `seq 1 11` `seq 20 28` -database ~/minirosetta_database/
```

### Example Overall Command Line (if overall protocol is run via a script or other program)
```
  grinder_dagman.py  -loop_start_pdb noloop_mini_1alc_H.pdb  -loop_res  `seq 12 19`   -align_pdb mini_1alc_H.pdb  -fasta mini_1alc.fasta  -nstruct 200  -native mini_1alc_H.pdb   -final_number 50 -denovo 1

[This is available in input/README_SETUP. ]

[The "seq" phrases make use of linux's seq command; for macs you either need to install this, or type out 12 13 14 15 16 17 18 19 for `seq 12 19`.]
```
The result should be:
```
  protein_build.dag 
```
An example is available in output/.

"dag" stands for "directed acyclic graph". This is in the format recognized by CONDOR's "dagman", which was the original queuing paradigm.  However, I found it really slow; lots of ltency in queuing jobs. So I ended up writing my own scripts to run on Stanford's BioX2 cluster, which uses LSF (Load Sharing Facility):
```
  SWA_pseudo_dagman_continuous.py  -j 20 protein_build.dag  

[This is available in input/README_SUB]
```
We have also written scripts to carry out the calculation on condor and torque clusters, and are almost finished with versions that use Amazon's EC2/S3. It usually takes a day or so to rewrite what we have for arbitrary systems. Our next step may be to write a version for Amazon's ElasticMapReduce, but it gets complicated since we actually have a series of map/reduce steps, not just one. 


## Versions

1. The right Rosetta version:

https://svn.rosettacommons.org/source/branches/das_lab/mini
Revision: 36561

[This branched off trunk in the winter of 2009.]

2. Several scripts to generate a loop modeling job are in:
 https://svn.rosettacommons.org/source/workspaces/rhiju/python
 Revision: 36561

[the scripts  include:
  grind_dagman.py
 stepwise_post_process_cluster.py
 stepwise_post_process_combine_and_filter_outfiles.py
 stepwise_pre_process_setup_dirs.py
 extract_lowscore_decoys.py

 and various helper python scripts...
]

 Also note that some of the python scripts look in ~rhiju for other python scripts -- I think if you change the paths at the top of grind_dagman.py to rosetta and to the python directory, you'll be good to go. If someone fixes this to be more generic, please let me know. We have to do it before publication anyway...

3. Finally, the queuing scripts are in:

 https://svn.rosettacommons.org/source/branches/das_lab/SWA_dagman_python 
 Revision: 36561

  I should probably combine the important scripts from #2 with the scripts in #3.

## References to published works using this protocol

This is unpublished (we're waiting for blind predictions to be tested... almost there!). Let us know if you extend the work, so we can trade tips and queuing strategies.

## Other Comments:
1. The run can be sped up if you use a subset of the protein around the loop of interest... 

2. To start extending the method to design, I think you just have to comment out some lines in src/protocols/swa/protein/StepWiseProteinPacker.cc:
    ```
    234		pack_task_->restrict_to_repacking();
    313		pack_task_->restrict_to_repacking();
    ```
    but obviously there will be more steps to calibrate the procedure for design applications.

3. We also have extensive scripts written for RNA building, (stepwise_rna_test.cc), and are planning to refactor and unify RNA/protein into one code flow in the Winter of 2010. 
Membrane Relax
==============

Author: Rebecca F. Alford (rfalford12@gmail.com)  
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)  
Last Updated: January 2015  
Rosetta Revision #58069 

---

High-resolution refinement is key for advancing low resolution structures from x-ray
crystallography to atomic level detail. For membrane proteins, this method can also
reveal an ensemble of possible membrane embeddings: the position and orientation of 
the biomolecule with respect to the membrane bilayer. 

The membrane relax application combines the Rosetta FastRelax algorithm with the
all atom energy function for membrane proteins and a gradient-based technique 
for optimizing the membrane embedding. First, a series of small backbone moves, 
rotamer trials, and minimization are used to refine the protein structure. In addition, 
the membrane position is optimizied by minimizing the "jump" or connecting relating
the MEM residue to the biomolecule. 

Publication describing the method: 
* Alford RF, Koehler Leman J, Weitzner BD, Duran A, Elazar A, Tiley D, Gray JJ 
  (2015) An integrated framework advancing membrane protein modeling and design 
  PLoS ONE (in preparation) 

## Executable/Script ##
The membrane framework relax application is implemented in Rosetta script. This script, 
called membrane_relax.xml is included in the main directory of this protocol capture. 

It can be run with the following executable: 

    Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease

## Generating Inputs ##
Two inputs are required for the membrane relax application: 

1. PDB for the protein structure of interest

2. Span file describing the location of trans-membrane spans

Steps for generating these inputs are found below. A set of example inputs can 
also be found in example_inputs/. Here, metarhodopsin II (PDB ID: 3pxo) is 
used as an example: 

1. PDB File: Generate a PDB file where the membrane protein structure is transformed 
   into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php).

2. Span File: Generate a spanfile from the PDB structure using
   the spanfile_from_pdb application described in the MP_spanfile-from-pdb protocol
   capture in Rosetta/demos/protocol_captures/2014. An example commandline using 
   3pxo is also provided here: 

        Rosetta/main/source/bin/spanfile_from_pdb.linuxgccrelease -database /path/to/db -in:file:s example_inputs/3pxo_tr.pdb

   For this example, this command will produce 1 output file: 
   * 3pxo_tr.span: Spanfile containing predicted trans-membrane spans

## Steps of the protocol ##
Here, we describe the steps required to run the MP_Relax protocol. As an example, all steps 
use the PDB 3pxo: 

1. Required Options: Options (flags) needed to run this application. A file with these flags, 
   relax_flags, is also provided for 3pxo in this demo: 

        flags                                  descriptions
        --------------------------------------------------------------------------------------------------
        -parser:protocol membrane_relax.xml    Use the membrane relax protocol Rosetta script
        -in:file:s                             Input PDB Structure: PDB file for protein structure
        -membrane_new:setup:spanfiles          Spanfile describing trans-membrane spans of the starting structure
        -membrane_new:scoring:hbond            Turn on membrane depth-dependent hydrogen bonding weight
        -relax:fast                            Use the FastRelax mode of Rosetta Relax (uses 5-8 repeat cycles)
        -relax:jump_move true                  Allow the MEM and other jumps to move during refinement
        -nstruct                               Number of structures to generate
        -packing:pack_missing_sidechains 0     Wait to pack until the membrane mode is turned on
        -out:pdb                               Output all PDB structures of refined models
        -out:file:scorefile                    Specify destination for score file

2. Recommended # of Decoys

   - For demo run: 1
   - For production runs: 1000

3. Command line: 

    To run this application, use the following command line: 

        Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease -database /path/to/db @relax_flags

Note on timing: Refinement in Rosetta is a time consuming application. Depending on avaialble
computing power and size of the protein, refinement of an individual decoy can take between 10-15min
for ~200 residues and between 0.5-1.0hrs for proteins > 200 residues. 

## Example Outputs ##
The following outputs will be generated from the relax protocol. A version of these outputs are also
provided in the example_outputs/ directory: 

* 3pxo_tr_0001.pdb      : Output refined model of 3pxo
* relax_scores_3pxo.sc  : Rosetta scores (including membrane scores) for refinement run

## References
1. Tyka MD, Keedy DA, Andre I, DiMaio F, Song Y, et al. (2011) Alternate states of proteins revealed by detailed energy landscape mapping. J Mol Biol. 

2. Barth P, Schonbrun J, Baker D (2007) Toward high-resolution prediction and design of transmembrane helical protein structures. Proc Natl Acad Sci 104: 15682–15687. 

3. Fleishman SJ, Leaver-Fay A, Corn JE, Strauch E-M, Khare SD, et al. (2011) RosettaScripts: A Scripting Language Interface to the Rosetta Macromolecular Modeling Suite. PLoS ONE 6: e20161. 
# Protocol Name:Loop Building for Membrane Proteins with Toplogy Broker

## Author
Elizabeth Dong


## Brief Description
The Topology Broker application in Rosetta can be used to build loops and flexible regions in membrane proteins. The protocol allows for the inclusion of membrane span and lipophilicity information as well as membrane weight set for the scoring function. This protocol rebuilds the intracellular helical region between TM 5 & 6 and C-terminus for Squid Rhodopsin (PDB ID: 2Z73A).


## Associated RosettaCon Talk 

### Title, Authors & Lab , Year, Session and Day of talk
  Elizabeth Dong, Anette Schreiber, Karen Gregory, Kristian Kauffman, Jeff Conn, Jens Meiler
  Meiler Lab, Vanderbilt University, Nashville, TN
  08/05/10, Session 9
  
### Abstract
	  Selective modulators of metabotropic glutamate receptor subtype 5 (mGluR5), a class C G-protein coupled receptor, provide novel treatment strategies for disorders that disrupt cognitive function. Identifying the specific residues on mGluR5 that contact these small molecules would allow for a deeper understanding of the binding interaction and aid in the development of therapeutic compounds. Construction of the mGluR5 model entailed identification of TM segments in the sequence using JUFO9D, modeling the 7 TM helices based on the three mammalian GPCR crystal structures using Rosetta and modeling the loops using the Toplogy Broker application in Rosetta, which allows for the inclusion of membrane span and lipophilicity information as well as the use of a membrane weight set to apply to the scoring function. Residues of mGluR5 critical for the binding of allosteric modulators were determined through Rosetta Ligand docking studies informed by experimental functional data. The experimentally validated models demonstrate the success of Rosetta to model GPCRs.
 
## Running

### Flags 
```
-run:protocol broker  #initiate call to broker
-broker:setup ./input_files/setup_broker.tpb #defines constraints on folding protocol
-frag3 ./input_files/aa2Z73A03_05.200_v1_3 #fragment files for folding of protein
-frag9 ./input_files/aa2Z73A09_05.200_v1_3 #fragment files for folding of protein
#Patches to the scoring function ensure that membrane potentials are used in the folding protocol
#make sure to have these either in the local directory or your database directory under scoring/weights
-stage2_patch score_membrane_s2.wts_patch
-stage3a_patch score_membrane_s3a.wts_patch
-stage3b_patch score_membrane_s3b.wts_patch
-stage4_patch score_membrane_s4.wts_patch
#allows setup of membrane options
-abinitio
	-membrane
#options for membrane scoring functions
-membrane
	-no_interpolate_Mpair
	-Menv_penalties
#tells folding protocol to close loops
-close_loops
-non_ideal_loop_closing
#not sure what this does but Yeifan included it.
-score
	-find_neighbors_3dgrid
-no_prof_info_in_silentout #no time-columns appears in score/silent - files
#input files these options should actually be supplied in the command line
-in
	-file
		-fasta ./input_files/2Z73A.fasta
		-spanfile ./input_files/2Z73A.span #generate from Octopus prediction (http://octopus.cbr.su.se/) using /TopologyBroker_GPCR/scripts/octopus2span.pl
		-lipofile ./input_files/2Z73A.lips4 #generate using /TopologyBroker_GPCR/scripts/run_lips.pl
#-out
#	-file
		#output options the results.silent_binary_out should be supplied on command line
		#-silent results.silent_binary_out
		#-silent_struct_type binary
#number of structures to generate supply on command line
#-nstruct 1
```

### Example Rosetta Command Line
```
r_broker.linuxgccrelease -database ~/minirosetta_database -out:file:residue_type_set centroid -out:file:silent rbroker_run1.out -nstruct 1 @flags.txt
```


## Versions
### svn revision number: 37327

## Other Comments: 
Before running the example, put all score_membrane*.wts_patch in the local directory

To generate *.span file: 
generate from Octopus prediction (http://octopus.cbr.su.se/) using /TopologyBroker_GPCR/scripts/octopus2span.pl

To generate *.lips4 file:
run the script /TopologyBroker_GPCR/scripts/run_lips.pl with the following command line:
```
run_lips.pl <fasta file> <span file> <path to blastpgp> <path to nr database> <path to alignblast.pl script>
```
example: 
```
run_lips.pl 2Z73A.fasta 2Z73A.span /sb/meiler/Linux2/x86/blast/blast-2.2.18/bin/blastpgp /sb/meiler/scripts/sequence_analysis/db/nr /mini/src/apps/public/membrane_abinitio/alignblast.pl 
```

Topology Broker only runs in centroid mode as of now. To extract pdb from *.out, use:
```
~/extract_pdbs.linuxgccrelease -database ~/minirosetta_database/ -in:file:silent rrbroker_run1.out -in:file:residue_type_set centroid -out:file:residue_type_set centroid -out:output
```
# Replica Docking

## Authors
This file was written in Jan 2013 by Zhe Zhang (zhe.zhang@tum.de) and corresponding PI is Oliver Lange (oliver.lange@tum.de)


## General Description
This demo contains all the files neccessary to replicate the results from the PLoS ONE RosettaCon
collection paper "Replica Exchange drastically improves sampling in low resolution docking stage
of RosettaDock" by Zhe Zhang, and Oliver Lange (2012).

All files, including the benchmark (`tar -xf dock_targetlib.tar.gz`) tested in the paper, and commands to setup the runs as well
 as post analysis are provided. For example outputs please refer to the silent-files (also trial.stat
from replica exchange runs) in example\_runs

## How to Run Demo
To run these demos:

1. export the directory where this protocol\_capture is located
    ```
    export PROTOCOL_CAPTURE=rosetta/rosetta_demos/protocol_capture/2012/
    ```

2. export your rosetta bin and database directory if you want to directly repeat the runs under directory example\_runs/ (need to first delete the outputted silent-files) instead of generating separate runs of your own. But this would only work for linux with slurm.
    ```
    export ROSETTA3_BIN=rosetta/rosetta_source/bin
    export ROSETTA3_DB=rosetta/rosetta_database
    ```

3. input pdbs preparation

   	all the targets are from Dockground benchmark3.0 (http://dockground.bioinformatics.ku.edu/UNBOUND/request\_new.php), in which unbound docking partners have been superimposed over its corresponding complex.
    ```
	get_pdb.py 1bvn_u1.pdb A	# this should output a file with only atom records from chain A, 1bvn_u1.pdbA.pdb
	replace_chain.py 1bvn_u2.pdb B > 1bvn_u2_B.pdb	# overwrite its chainID to B, 1bvn_u2_B.pdb
	get_pdb.py 1bvn_u2_B.pdb B	# 1bvn_u2_B.pdbB.pdb

	cat 1bvn_u1.pdbA.pdb 1bvn_u2_B.pdbB.pdb > protAB.pdb	# superimposed native structure, used for rmsd related calculation.
	scripts/initial_randomize.sh	# output P.pdb, both docking partners are randomly reoriented, and then slided into contact.
					# used as the input pdb for docking.
    ```
	write flag "-partners A_B" into file "partners" for later use.
    ```
	scripts/get_disulf_pairs.sh	# get the disulfide residue pairs, later used in refinement
    ```
	All the input files in dock_targetslib (`tar -xf dock_targetlib.tar.gz`)  are prepared in this way.

4. setup target library using automated setup tools available with the CS-Rosetta toolbox(www.csrosetta.org).

	This step assembles target related input files to build the target library (`dock_targetlib` as an example). By default 
    the library is stored in folder cs\_targetlib at home of your workspace. You can also specify a directory using the flag 
    `-target_prefix` as follows. Absolute path is recommended for `-target_prefix`.

	The following commands for setup demo target 1bvn have been wraped up in `scripts/test_target_setup.sh`
    ```
	#RosettaDock:
	#setup target for RosettaDock as in published paper "Protein-Protein Docking with Simutaneous Optimization of Rigid-body 
	#Displacement and Side-chain Conformations", Jeffrey J. Gray et al., J. Mol. Biol. (2003)
	setup_target -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
    -disulf disulf_file -native protAB.pdb -pdb P.pdb -partners partners

	#ReplicaDock:
	#setup target for ReplicaDock as described in Zhang and Lange, PLOS One 2012.
	setup_target -method replica_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
    -native protAB.pdb -pdb P.pdb -partners partners

    #or you can conviently copy the inputs from a previously prepared 'rosetta_dock' setup as follows:
    setup_target -method replica_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
    -transfer_method rosetta_dock
    ```

5. setup run using automated setup tools. This step creates a run-ready directory as specified with flag `-dir`, in which job-scripts, 
    input files, RosettaScripts-xml as well as flag files are contained. For flag `-dir`, absolute path is recommended. For job scripts, 
    you can use different types (e.g. moab) according to your queuing system.

	In example\_runs additional comments have been added subsequently to guide you through the automatically generated input files. 
    All job scripts in example\_runs are for queuing system with slurm.

	The commands to setup run as follows for single-mache/interactive mode are wraped up in `scripts/test_run_setup.sh`. After run 
    `scripts/test_run_setup.sh`, use `source production.interactive.job -n $Np` with $Np specifying the processor numbers to start the run 
    in the run directory, (for example `$PROTOCOL_CAPTURE/test_runs/replica_dock/udock_1bvn/run` after you run `scripts/test_run_setup.sh`)

	For quick test only purpose for ReplicaDock/ReplicaDock-LoT, I have set the total trajectory length to be 5000 MC-steps in 
    RosettaScripts(dock_cen.xml, dock_cen_7.xml) in the automated setup tool box, which can be simply found in 
    `csrosetta3/flag_library/methods/_docking_base/` and modified accordingly for production purpose use.

    1. centroid stage
	    1. ReplicaDock, using temperatures [2.0 3.0 5.0]

	        ReplicaDock is run in MPI-mode using RosettaScript to compile: 
            ```
            ./scons.py -j 48 bin/rosetta_scripts.mpi.linuxgccrelease mode=release extras=mpi
            ```

	        Please note that specific numbers of processors have to be used: calculate number of processes using the formula: nstruct \* n_replica + 2. 
            The extra 2 processes are dedicated to the job distributor and File IO. n_replica is the number of temperature levels (here 3), and nstruct 
            can be any positive integer. ReplicaDock outputs the trajectory in the form of two silent-files: one containing decoy+score information 
            (name: decoys\_<input_pdb>\_nnnn\_traj.out), and the second file is a copy of just the score information (scores\_<input_pdb>\_nnnn\_traj.out). 
            Decoytags are of the form P\_tttt\_rrr\_ssssssss where tttt informs about trajctory number, rrr about the replica, and ssssssss about the 
            snapshot number within the trajectory. The temperature\_levels are switched between different replicas. The current temp\_level or temperature
            of a replica at the moment a decoy was recorded is found in the score-columns temp\_level and temperature.
            At the end of a trajectory the final decoy is written to the file 'decoys.out'; this file is a relict of using the JD2-framework and
            can be generally ignored.
            Additionally, the file 'trial.stat' is produced which gives information about acceptance rates in each temperature level.
            ```
            setup_run -method replica_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
            -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/replica_dock -job interactive -extras mpi -score interchain_cen \
            -nstruct 1 -protocol rep_cen -xml uniform -n_replica 3
            ```

            start running:
            ```
            cd $PROTOCOL_CAPTURE/replica_docking/test_runs/replica_dock/udock_1bvn/run/;
            source production.interactive.job -n 5
            ```

	    2. ReplicaDock-LoT, using temperatures [0.6 0.8 1.0 1.2 1.5 2.0 2.5] + min_score

	        Min_score is used to flatten the score function. A reasonble min_score value is determined as the average score of the first 50 
            snapshots of temperature 1.5 when simulated without min-score. As shown in example_runs/replica_dock_LoT/get_min_score/udock_1bvn/run/, 
            several trajectories are run and an average value of -36.519 can be get by:
            ```
            cat scores_P_000* >scores_traj.fsc
            silent_data.py scores_traj.fsc temperature score | awk '$1==1.5{print}' | median.py

            setup_run -method replica_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
            -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/replica_dock_LoT -job interactive -extras mpi -score interchain_cen \
            -min_score -36.519 -nstruct 1 -protocol rep_cen -xml uniform -n_replica 7
            ```

	    3. RosettaDock's original low-resolution stage (shotgun sampling)
            ```
	        setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
            -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/rosetta_dock -job interactive -extras mpi -protocol centroid -batches 2 \
            -score interchain_cen -nstruct 25
            ```

    2. refinement. To refine the decoys generated in the centroid stage, we don't want to copy all the decoy-files into the run directory of refinement, but only specify the path of the decoys files. For this, we use the flag `-start` with absolute path specified together with flag `-pattern` to only include the files with a certain name pattern.

	    1. refinement of the low-resolution ensembles produced by ReplicaDock
        ```
	    setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
        -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/refine_replica_dock -job interactive -extras mpi -protocol refine \
        -pattern "low_decoys_*out" -prefix refine -score docking -nstruct 1 \
        -start $PROTOCOL_CAPTURE/replica_docking/example_runs/replica_dock/udock_1bvn/run/
        ```

	    2. refinement of the low-resolution ensembles produced by RosettaDock's shotgun sampling approach
        ```
	    setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
        -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/refine_rosetta_dock -job interactive -extras mpi -protocol refine \
        -pattern "low_decoys_*out" -prefix refine -score docking -nstruct 1 \
        -start $PROTOCOL_CAPTURE/replica_docking/example_runs/rosetta_dock/udock_1bvn/run
        ```

	    3. refinement of the low-resolution ensembles produced by ReplicaDock-LoT
        ```
	    setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
        -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/refine_replica_LoT -job interactive -extras mpi -protocol refine \
        -pattern "low_decoys_*out" -prefix refine -score docking -nstruct 1 \
        -start $PROTOCOL_CAPTURE/replica_docking/example_runs/replica_dock_LoT/replica_with_minscore/udock_1bvn/run
        ```

	    4. generate RelaxedNative ensembles with 1000 decoys from the superimposed native structure protAB.pdb, used as reference in final analysis
        ```
	    setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
        -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/relax_native -job interactive -extras mpi -protocol refine \
        -out relax_native.out -score docking -nstruct 1000
        ```

6. post-filter of low resolution ensembles. To save space, some .out files have been deleted and only score files shown.

	1. RosettaDock's shotgun approach sampled ensembles:
	use N to denote the total number of decoys produced from RosettaDock low resolution phase. Then we first exclude the decoys with interchain\_contact>10, then select N\*0.36 decoys by score.

	For an example decoys file in `$PROTOCOL_CAPTURE/replica_docking/example_runs/rosetta_dock/udock_1bvn/run/decoys.out`, we do:
    ```
	# get the tags of the decoys selected. 50 decoys in all in the example file, so select 50*0.36=18 decoys 
	# for analysis and further refinement
	for i in $(ls decoys_000?.out); do echo $i; cat $i | grep SCORE: >> decoys.fsc; done

	scripts/silent_data.py decoys.fsc score interchain_contact description | awk '$2<=10{print}' | sort -n -k 1 | head -n 18 \
    | awk '{print $3}' > tag_low

	# extract the selected decoys from decoys.out
	for i in $(ls decoys_000?.out); do echo $i; scripts/extract_tagged_decoys.py $i tag_low > low_$i; done
    ```

	2. replica_dock: T=[2.0 3.0 5.0]
	exclude decoys with temperature 5.0, then exclude decoys with interchain\_contact > 10
    ```
	# get the tags of the snapshots to refine
	cat scores_P_000*fsc > scores_traj.fsc
	scripts/silent_data.py scores_traj.fsc temperature interchain_contact description | awk '$1<4&&$2<=10{print $3}' > tag_low

	# extract selected decoys from the trajectory silent file
	for i in $(ls decoys_P_000?_traj.out); do echo $i; scripts/extract_tagged_decoys.py $i tag_low > low_$i; done
    ```

	3. replica_dock_LoT: T=[0.6 0.8 1.0 1.2 1.5 2.0 2.5] and min_score
	first exclude decoys with interchain\_contact > 10, then keep the 0.5\*N lowest scoring decoys from the remaining set.
    For the example silent-files of trajectory-snapshotse in `$PROTOCOL_CAPTURE/replica_docking/example_runs/replica_dock_LoT/replica_with_minscore/udock_1bvn/run/`, it is done as follows:
    ```
	# get tag. 14 decoys in all for the
	cat scores_* > scores_traj.fsc
	scripts/silent_data.py scores_traj.fsc interchain_contact score description | awk '$1<=10{print}' | sort -n -k 2 \
    | head -n 7 | awk '{print $3}' > tag_low

	# extract selected decoys from the trajectory silent file
	for i in $(ls decoys_P_000?_traj.out); do echo $i; scripts/extract_tagged_decoys.py $i tag_low > low_$i; done
    ```

7. refinement of low resolution ensembles

	refinement is setup with the automated tool box as shown in part 5.2.1-5.2.4.

8. analysis

	scripts for quick checking and analysis are given, i.e. `scripts/outfile_plot.py` and `scripts/hist.py`. To run these two scripts, extra python library matplotlib.pyplot is required.
    ```
	scripts/outfile_plot.py #generates scatter plots for given a silent_file with specified score terms, e.g. rms vs. score
	scripts/hist.py #generates histogram figure for given silent_files with specified score term
    ```
	
###
###
### This file is part of the CS-Rosetta Toolbox and is made available under
### GNU General Public License
### Copyright (C) 2011-2012 Oliver Lange
### email: oliver.lange@tum.de
### web: www.csrosetta.org
###
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <http://www.gnu.org/licenses/>.
###
###
This module is copied from biopython-1.57 and has been cleaned up in a hack and slash fashion to reduce external dependencies. This means that probably many parts of the API 
are now broken. We only use this for PDB import/export and could probably reduce the number of files here drastically or use a different PDB reader/write entirely.

If you want to use any functionality of biophython we suggest to install the full package and import via Bio.PDB. 


Here is the Biophython License Statement

                Biopython License Agreement

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.

Here is a list of Contributors to the origian BioPython package: 
CONTRIBUTORS
============

This is a list of people who have made contributions to Biopython.
This is certainly not comprehensive, and if you've been overlooked
(sorry!), please mention it on the development mailing list.
People are listed alphabetically by surname.

Cecilia Alsmark <Cecilia.Alsmark at domain ebc.uu.se>
Tiago Antao <tiagoantao at gmail.com>
Sebastian Bassi <sbassi at domain asalup.org>
Bill Barnard <bill at domain barnard-engineering.com>
Yves Bastide <ybastide at domain irisa.fr>
Paul T. Bathen
Yair Benita <Y.Benita at domain pharm.uu.nl>
Peter Bienstman <Peter.Bienstman at domain rug.ac.be>
Jose Blanca
Bob Bussell <rgb2003 at domain med.cornell.edu>
Diego Brouard <diego at domain conysis.com>
James Casbon <j.a.casbon at domain qmul.ac.uk>
Hye-Shik Chang <perky at domain fallin.lv>
Jeffrey Chang <jchang at domain smi.stanford.edu>
Brad Chapman <chapmanb at domain arches.uga.edu>
Peter Cock <p.j.a.cock at googlemail dot com>
Marc Colosimo <mcolosimo at domain mitre.org>
Andres Colubri <andres dot colubri at gmail dot com>
Cymon J Cox <cymon at domain duke.edu>
Gavin E Crooks <gec at domain compbio.berkeley.edu>
Andrew Dalke <dalke at domain acm.org>
Michiel de Hoon <mdehoon at domain c2b2.columbia.edu>
Bart de Koning <bratdaking gmail>
Sjoerd de Vries <sjoerd at domain nmr.chem.uu.nl>
Nathan J. Edwards <nje5 at edu domain georgetown>
Kyle Ellrott
Jeffrey Finkelstein <jeffrey.finkelstein at domain gmail.com>
Iddo Friedberg <idoerg at domain burnham.org>
Bertrand Frottier <bertrand.frottier at domain free.fr>
Phillip Garland <pgarland at gmail>
Walter Gillett < https://github.com/wgillett >
Frederik Gwinner
Jason A. Hackney <jhackney at domain stanford.edu>
Thomas Hamelryck <thamelry at domain binf.ku.dk>
Michael Hoffman <hoffman+biopython at domain ebi.ac.uk>
Thomas Holder
Yu Huang <krocea at domain yahoo.com.cn>

#                 Biopython License Agreement
#
# Permission to use, copy, modify, and distribute this software and its
# documentation with or without modifications and for any purpose and
# without fee is hereby granted, provided that any copyright notices
# appear in all copies and that both those copyright notices and this
# permission notice appear in supporting documentation, and that the
# names of the contributors or copyright holders not be used in
# advertising or publicity pertaining to distribution of the software
# without specific prior permission.
#
# THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
# WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
# CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
# OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
# OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
# OR PERFORMANCE OF THIS SOFTWARE.
Kevin Jacobs <jacobs at bioinformed dot com>
Diana Jaunzeikare
Joanna & Dominik Kasprzak
Frank Kauff <fkauff at domain duke.edu>
Siong Kong
Andreas Kuntzagk <andreas.kuntzagk at domain mdc-berlin.de>
Michal Kurowski <michal at domain genesilico.pl>
Uri Laserson <laserson at Massachusetts Institute of Technology dot edu>
Chris Lasher <chris.lasher at gmail.com>
Gaetan Lehman <gaetan.lehmann at domain jouy.inra.fr>
Katharine Lindner <katel at domain worldpath.net>
Erick Matsen <surname at fhcrc dot org>
Tarjei Mikkelsen <tarjei at domain genome.wi.mit.edu>
Konstantin Okonechnikov <k.okonechnikov at domain gmail.com>
Cheng Soon Ong <chengsoon.ong at tuebingen.mpg.de>
Anne Pajon <ap one two at sanger ac uk>
Claude Paroz <claude at two (as digit) xlibre dot net>
Andrea Pierleoni <andrea at the Italian domain biocomp dot unibo>
Mike Poidinger <Michael.Poidinger at domain eBioinformatics.com>
Leighton Pritchard <lpritc at domain scri.sari.ac.uk>
Wolfgang Schueler <wolfgang at domain proceryon.at>
Peter Slickers <piet at domain clondiag.com>
Thomas Sicheritz-Ponten <thomas at domain cbs.dtu.dk>
Frederic Sohm <fsms at domain users.sourceforge.net>
Joao Rodrigues <anaryin at the domain gmail dot com>
Thomas Rosleff Soerensen <rosleff at domain mpiz-koeln.mpg.de>
Eric Talevich <eric.talevich at domain gmail.com>
Bartosz Telenczuk <bartosz.telenczuk at domain gmail.com>
Carlos Rios Vera <crosvera at domain gmail.com>
Johann Visagie <wjv at domain cityip.co.za>
Dan Vogel <dmv at domain andrew.cmu.edu>
David Weisman <david.weisman at domain acm.org>
Bartek Wilczynski <bartek at domain rezolwenta.eu.org>
David Winter <david dot winter at gmail dot com> 
Hongbo Zhu
Christian Zmasek
Harry Zuzan <iliketobicycle at domain yahoo.ca>


Broker Protocol Captures
========================

This directory contains 3 protocol captures:

* [[domain_insertion|/protocol_capture/broker/domain_insertion/README]]
* [[snugdock|/protocol_capture/broker/snugdock/README]]
* [[ubq|/protocol_capture/broker/ubq/README]]
Domain Insertion Protocol Capture
=================================

This protocol runs ab inito for a domain insertion protein per the Broker 
paper. Because the loop closure for this protocol is particularly difficult,
and because the loop closure algorithm in this example relies on fragments
to perform the closure, the loop closure succeeds rather infrequently. The
problem is that, for many non-native loops, fragment coverage of the
required region of conformational space does not exist.

As a result, loop closure failed for me in >90% of cases. Thus, to test
this protocol, you must run with -n 100 or greater. Alternatively, the
closure is always able to close the loops in the native structure, so if
the sampling protocol starts near a correct (or probably any closable)
structure, so you can run with -s native.pdb to produce a closed structure
nearly every time.

An additional tip for producing output more quickly is to disable (or
turn down) the relaxation step at the end. Because relaxation is done
system-wide and not exclusively for the inserted domain, it can have
a very long runtime.

Assuming Rosetta is in your home directory, you run it as follows:

    $ ~/Rosetta/main/source/bin/rosetta_scripts.default.[platform][compiler]release @flags -nstruct [number of structures]

where platform is "linux" or "mac" and compiler is "gcc" or "clang".
AbInitio Protocol Capture
=========================

This protocol runs abinitio with a given beta strand topology and a helix fixed 
at atomic resolution.

Assuming Rosetta is in your home directory, you run it as follows:

    $ ~/Rosetta/main/source/bin/rosetta_scripts.default.[platform][compiler]release @flags -nstruct [number of structures]

where platform is "linux" or "mac" and compiler is "gcc" or "clang". 10000 
structures should be more than enough to see what's going on.
SnugDock Protocol Capture
=========================

This protocol runs our snugdock-inspired antibody docking protocol. It uses 
fragment insertion instead of normal loop modeling, but is otherwise very 
similar to the orginal.

Assuming Rosetta is in your home directory, you run it as follows:

    $ ~/Rosetta/main/source/bin/rosetta_scripts.default.[platform][compiler]release @flags -nstruct [number of structures]

where platform is "linux" or "mac" and compiler is "gcc" or "clang".
GPCR Model Dock Protocol Capture
================================

Written in Jan 2013 by Elizabeth Nguyen (e dot dong dot nguyen at vanderbilt dot edu)

---

This demo contains all the files necessary to replicate the results from the 
paper:

Challenges of recovering native conformations of ligands docked into 
comparative models of G-protein coupled receptors. Nguyen,E.D.*, Norn, C.*, 
Frimurer, T., and Meiler,J. (2012), PLOS One. [[Supplementary 
Material|/protocol_capture/gpcr_model_dock/SupplementaryMaterial_ProtocolCapture.docx]].

Commands and files are provided to be able to replicate the following steps 
(from Figure 1):

1. Structural alignment of GPCR templates
2. Sequence alignment of the target GPCR to templates
3. Thread target sequence onto template backbone coordinates
4. Rebuild missing density
5. Rebuild ECL 1,2 and 3
6. Evaluate comparative models by clustering by full-receptor RMSD and 
   knowledge-based pocket residue filter
7. Generate ligand conformations in MOE
8. Dock ligand into comparative models
9. Analyze results by clustering binding modes by ligand RMSD


Files included in this protocol capture
---------------------------------------

Input files (rosetta_inputs):

    1u19A_clean.pdb 
    2vt4A_clean.pdb 
    2rh1A_clean.pdb 
    3emlA_clean.pdb 
    3oduA_clean.pdb 
    3pblA_clean.pdb 
    3rzeA_clean.pdb 
    3v2wA_clean.pdb 
    3uonA_clean.pdb 
    4dajA_clean.pdb 
    4dklA_clean.pdb 
    4djhA_clean.pdb 
    4ea3A_clean.pdb 
    4ej4A_clean.pdb
    all_gpcrs.fasta 
    1u19A.fasta
    1u19A.aln
    2rh1A_clean.pdb
    1u19A_on_2rh1A.pdb
    1u19A.jufo_ss
    1u19A.psipred_ss2
    1u19A.span
    1u19A.disulfide
    aa1u19A03_05.200_v1_3
    aa1u19A09_05.200_v1_3
    relax.options
    1u19A_on_2rh1A_relax.pdb
    ccd_initial.options
    1u19A_on_2rh1A.loops
    1u19A_on_2rh1A_initial.pdb
    1u19A_rmsd01.pdb
    1u19A.sdf
    1u19A.params
    1u19A_confs.pdb 
    1u19A_cluster01_01.pdb
    1u19A_cluster01_01_ligand.pdb
    dock.options
    dock.xml 
    1u19A_ligand.cluster.mat

Output files (example_outputs):

    1u19A_10percent_RMSD.txt
    cluster3_1u19A.Centers
    cluster3_1u19A.Rows 
    1u19A_cluster01_01_ligand_011u19A_cluster01_01_ligand_0001.pdb
    all.sdf
    cluster3_1u19A_ligand.Centers
    cluster3_1u19A_ligand.Rows

Scripts and applications not included in the Rosetta 3.4 release (scripts):

    evaluate_score_vs_pocket_rmsd
    jufo9d_span.pl
    rmsd.tcsh

Command lines to run this protocol capture
------------------------------------------

* Prepare GPCR crystal structures from the Protein Data Bank.

        rosetta_tools/protein_tools/scripts/clean_pdb.py 2RH1 A > 2rh1A_clean.pdb

* Perform a structural alignment of GPCRs using crystal structures from the Protein Data Bank. 

        mustang -p . -i 1u19A_clean.pdb 2vt4A_clean.pdb 2rh1A_clean.pdb 3emlA_clean.pdb 3oduA_clean.pdb 3pblA_clean.pdb 3rzeA_clean.pdb 3v2wA_clean.pdb 3uonA_clean.pdb 4dajA_clean.pdb 4dklA_clean.pdb 4djhA_clean.pdb 4ea3A_clean.pdb 4ej4A_clean.pdb -o all_gpcrs -F fasta -D 2.5 -s ON

* Sequence alignment of the target GPCR to templates

  Input target sequence 1u19A.fasta and profile alignment all_gpcrs.fasta to 
  http://mobyle.pasteur.fr/cgi-bin/portal.py#forms::clustalO-profile. We used 
  the default settings for this protocol.

* Thread target sequence onto template backbone coordinates

    rosetta_tools/protein_tools/scripts/thread_pdb_from_alignment.py --template=2rh1A_clean --target=1u19A --chain=A --align_format=clustal 1u19A.aln 2rh1A_clean.pdb 1u19A_on_2rh1A.pdb

* Generate secondary structure prediction, constraint file and fragments for bRh. 

  * Secondary structure- Jufo9D: http://meilerlab.org/index.php/servers/show?s_id=5

  * Secondary structure- PSIPRED: http://bioinf.cs.ucl.ac.uk/psipred/

  * Transmembrane span prediction based on Jufo9D:

            perl scripts/jufo9d_span.pl 1u19A.jufo9d > 1u19A.span 

  * Disulfide bond constraint file: Create file that lists residue number of 
    cysteine residues predicted to disulfide bond according to the alignment 
    with the template.

  * Fragment files: http://www.robetta.org

* Rebuilt missing density

        rosetta_source/bin/loopmodel.linuxgccrelease @ccd_initial.options -database rosetta_database 

* Rebuilt ECL 1,2 and 3 with CCD

        rosetta_source/bin/loopmodel.linuxgccrelease @ccd.options -database rosetta_database 

* Rebuilt ECL 1,2 and 3 with KIC

        rosetta_source/bin/loopmodel.linuxgccrelease @kic.options -database rosetta_database

* Analyze results by clustering top ten percent of comparative models by full receptor RMSD.

        bcl.exe PDBCompare -quality RMSD -atoms CA -pdb_list 1u19A_models.ls -aaclass AACaCb -prefix 1u19A_10percent_
        bcl.exe Cluster -distance_input_file 1u19A_10percent_RMSD.txt -input_format TableLowerTriangle -output_format Rows Centers -output_file cluster3_1u19A -linkage Average -remove_internally_similar_nodes 3

* Analyze results by filtering comparative models with a knowledge-based filter.

        scripts/evaluate_score_vs_pocket_rmsd/01_make_distances.csh
        scripts/evaluate_score_vs_pocket_rmsd/02_filter_models.py

* Create ligand conformations in MOE.

  See MOE operating guide. LowModeMD with the MMFFx94 force field and 
  Generalized Born solvation model was used to generate conformations within 
  the specified energy cutoff. The ligand conformations were then saved as an 
  .sdf file for conversion to .pdb and .params files for Rosetta.

        rosetta_source/src/python/apps/public/molfile_to_params.py -n 1u19A -p 1u19A 1u19A.sdf 
 
* Dock ligand into comparative models.

        rosetta_source/bin/rosettascripts.linuxgccrelease @dock.options -database rosetta_database

* Filter binding modes by energy, clustering and experimental restraints

        /scripts/rmsd.tcsh *.pdb
        bcl.exe ScoreSmallMolecule all.sdf output.sdf -comparison RMSD
        bcl.exe Cluster -distance_input_file 1u19A_ligand.cluster.mat -input_format TableLowerTriangle -output_format Rows Centers -output_file cluster3_1u19A_ligand -linkage Average -remove_internally_similar_nodes 3
Ligand-centric Water Docking
============================

Author: Gordon Lemmon  
Citation:
* Gordon Lemmon, Jens Meiler (2012). Toward ligand docking including explicit 
  interface water molecules. PLoS ONE (submitted).

---

Small molecule docking predicts the interaction of a small molecule ligand with 
a protein at atomic-detail accuracy including position and conformation the 
ligand but also conformational changes of the protein upon ligand binding. 
While successful in the majority of cases, leading docking algorithms including 
RosettaLigand fail in some cases to predict the correct protein/ligand complex 
structure. In this study we show that simultaneous docking of explicit 
interface water molecules greatly improves Rosetta’s ability to distinguish 
correct from incorrect ligand poses. This result holds true for both 
protein-centric water docking wherein waters are located relative to the 
protein binding site and ligand-centric water docking wherein waters move with 
the ligand during docking. Protein-centric docking is used to model 99 HIV-1 
protease/protease inhibitor structures. We find protease inhibitor placement 
improving at a ratio of 9:1 when one critical interface water molecule is 
included in the docking simulation. Ligand-centric docking is applied to 341 
structures from the CSAR benchmark of diverse protein/ligand complexes. Across 
this diverse dataset we see up to 56% recovery of failed docking studies, when 
waters are included in the docking simulation.

Purpose and algorithm
---------------------

The simultaneous docking of ligands and water molecules is now possible within Rosetta
This work is being submitted to the PLoS ONE Special collection from RosettaCon 2012

Waters can be docked in two ways. The first way is refered to in the paper as protein-centric 
water docking. Protein centric waters sample the protein binding pocket independent of the ligand.
Ligand-centric waters are positioned initially around the ligand. As the ligand translates and
rotates about the protein binding site, the waters move with it, maintaining their positions relative
to the protein. After low-resolution docking of the ligand. The waters undergo their own independent
but smaller movements.

Tools and Input Files
---------------------

#### Scripts:

* Use `find_waters_pymol.py` to identify waters within the protein/ligand 
  interface. The script takes 3 arguments, each specifying a PDB file 
  (protein.pdb, ligand.pdb, water.pdb). This script requires that you have 
  installed pymol as a python library

* The script `ligand_properties_from_bcl` requires that you download and 
  install BCL, available from the Meiler Lab website, www.meilerlab.org

#### Required flags:

    -in:file:s <pdb file> # a starting structure upon which docking will be performed. Should contain a protein, a ligand, and one or more waters
    -in:file:extra_res_fa # the .params file for your ligand. This is created by providing a .mol file to the script: rosetta_source/src/python/apps/public/molfile_to_params.py
    -treat_residues_with_this_chain_as_separate_chemical_entities <1 letter chain from PDB> # Useful for giving waters the same chain and having Rosetta treat them separately.

#### Optional flags:

    -ex1, -ex1aro, and -ex2 # expand the rotamer sets that are sampled during packing.
    -in:file:native # allows calculation of comparison metrics between Rosetta models and the correct pose if this is known

#### Example Rosetta Command Line (Use Rosetta3.5 or revision 48472):

The three XML files, standard.xml, protein_centric.xml and ligand_centric.xml, demonstrate how to dock waters using
RosettaLigand. Simply run the following command from the directory where this readme is found:

    <path to rosetta_source>/bin/rosetta_scripts.linuxgccrelease @rosetta_inputs/flags.txt -parser:protocol <choice of xml file> -database <path to rosetta_database>

Expected Outputs
----------------

This demo only produces one structure. Add "-nstruct 1000" to produce the amount of sampling used in this paper.
In our paper we sort the top 1000 by total score, then the top 100 by interface score. The top model by interface
score is used as our most likely prediction of ligand pose

To sort models by total score simply:

    grep -H  total_score *.pdb | sort -nk 2 | head -n 100 | cut -d ':' -f 1 > top100.txt

To sort by interface score:

    grep -H interface_delta_X *.pdb `cat top100.txt` | sort -nk 2 | head -n 1
﻿Anchored Design
===============
 This document describes how to use the AnchoredDesign protocol, both in benchmarking and design mode.  As the protocol's components [AnchorFinder](http://www.rosettacommons.org/docs/latest/anchor-finder.html), [AnchoredPDBCreator](http://www.rosettacommons.org/docs/latest/anchored-pdb-creator.html), and [AnchoredDesign](http://www.rosettacommons.org/docs/latest/anchored-design.html) are reasonably extensively documented elsewhere, this protocol capture is meant to be used alongside that online documentation. Presented at RosettaCon2010 (in poster form) was a description of the protocol itself, plus benchmarking results, plus some early design results.  The accompanying paper ([Lewis SM, Kuhlman BA. Anchored design of protein-protein interfaces. PLoS One. 2011;6(6):e20872. Epub 2011 Jun 17.](http://www.ncbi.nlm.nih.gov/pubmed/21698112) (pubmed link)) describes only benchmarking results, but the tools to do design are described here.  A paper on design results is forthcoming.

Note that this protocol capture is somewhat focused on just the AnchorFinder portion (the least important part of the process), because the other portions are documented elsewhere but AnchorFinder largely is not.  


Contained here:
* Instructions on choosing appropriate benchmarks – AnchorFinder or otherwise
* Instructions on preparing those structures – for benchmarking or via AnchoredPDBCreator
* How to benchmark AnchoredDesign against those structures
* How to use AnchoredDesign to design interfaces
* command lines/option files, and discussion of the options

Not contained here:
* A speck of code – that lives in your Rosetta release.
* Submission scripts for running jobs.  I don't know your architecture.  It's all MPI compatible so it's not hard.

Sort-of contained here:
* The raw_documentation directory includes copies of the doxygen-style manual documentation for AnchorFinder, AnchoredDesign and AnchoredPDBCreator.  These copies are guaranteed NOT to be up to date; look in your copy of the code or the web links above instead.

Overview
--------
The purpose of the AnchoredDesign protocol is to create new protein-protein interactions using information (the anchor) borrowed from a known interaction at the same interface of one partner.  Because this protocol is intended to design protein-protein interactions, the obvious test is to see whether it can recover known structures of such interactions.  This document is massively overconcerned with benchmarking because it accompanies the paper in which AnchoredDesign is benchmarked; unless you are actually trying to replicate the benchmarking you can ignore most of those details and skip to the design tools.

The protocol modifies a loop region around an anchor in designing binders.   Selection of structures for benchmarking therefore requires interface loops with anchor regions.

The ideal anchor has several qualities:
* one or a few contiguous residues – the protocol can only have one anchor
* many interactions across an interface – this is the whole point of an anchor
* not embedded in the middle of a secondary structure element (helix/strand) – the anchor needs to be in a loop because the interface space will be sampled via loop remodeling of the anchor loop

Examples of anchors might be a phosphotyrosine inserting into an SH2 domain, a polyproline sequence binding an SH3 domain, etc.

For design, one would choose an anchor based on one's target.  For benchmarking, you are free to choose anything that has an interface loop with a good anchor.  To select benchmarking structures, I wrote the AnchorFinder protocol.  AnchorFinder has some value in highlighting which residues might make good anchors for a given target (although computational alanine scanning, not covered here, is more likely to be useful).

After finding suitable structures with the help of AnchorFinder, the next step is to pick anchors and loops out of those structures, and in general prepare them for Rosetta's use (removing solvent atoms, etc).  At this point you're ready to run AnchoredDesign.


Compiling AnchorFinder
----------------------
(This section applies to benchmarking only)

Your Rosetta code distribution should include an application called AnchorFinder. If you wish to search large numbers of PDBs for potential anchors (I searched a local copy of the entire PDB structure set), then you will wish to modify the code slightly before running it.  Running any part of Rosetta against huge numbers of unprepared, straight-from-the-PDB structures is challenging because the PDB reader in Rosetta is not robust against nonstandard file formats, etc.

To compile AnchorFinder such that it will be robust, examine the manual documentation on RobustRosetta (also included in this protocol capture).  Briefly, this documentation describes changes that A) make Rosetta slower (thus they aren't on by default) and B) cause it to throw C++ exceptions when it hits errors instead of crashing.  The job distributor catches the errors, skips the bad structures, and continues.  You must recompile after making these changes.

You do not want to use compiled executeables OTHER than AnchorFinder with these changes made – they will significantly slow the code down.  AnchorFinder is quite fast so it's not a problem.

When running AnchorFinder, watch your memory usage.  When I used it, there was a patch in the JobDistributor which deleted starting poses for PDBs that had already been processed.  This patch was rejected by the community and since been replaced by a different patch to do the same thing; AnchorFinder is a run-once sort of thing so it has not been tested against the new method.

Using AnchorFinder
------------------
(This section is minimally relevant if not benchmarking)

At this point you should have a compiled copy of AnchorFinder with the necessary changes to the code.  You can then list your PDBs in one or many -l files (or -s) for use in Rosetta.  The format for Rosetta's -l flag is one path per line:

    A.pdb
    B.pdb
    C.pdb
    ...

Depending on your available architecture, it may be better to split the run up into many -l on separate processors.  I don't know what's best for you.

If you want to do the whole PDB – it's a good idea to skip the largest PDB files ahead of time, particularly ribosome structures.  These take a very long time to process through the PDB reader, and due to heavy nucleic acid content are skipped anyway.  You can also either toss the NMR structures ahead of time or use the -obey_ENDMDL flag to only read the first model.

AnchorFinder will automatically remove nonprotein atoms from the Poses before examination.  It also skips anything that is monomeric, has no protein residues, or smaller than 20 residues after processing.

It will then look through the structures searching for regions with certain command-line-defined characteristics.  These characters are:
* length of windows for consideration - 4 or 5 or 6, etc, contiguous residues at a time?  This flag is -window_size.  I suggest 5 residue windows.
* What fraction of this window should have loop secondary structure as assigned by DSSP?  -loopness controls this.  It takes it as a decimal between 0-1, I suggest 0.6 (which translates to 3/5 residues for a 5 residue window) to 1 (all residues loop).
* How many cross-interface interactions per residue should the window have?  I suggest a minimum of 4.  This translates to 20 (redundancy included) cross-interface interactions for a 5 residue window.  By redundancy, I mean residues 43 and 44 on chain A can both interact with residue 234 on chain B and it will count as two interactions.  Specify this with -nbrs_per_residue.
* What file name should the good interactions be printed to?  I leave it as an exercise to the reader to pick their own file name.  Specified with -bestoutfile; defaults to goodfile.out.

Running AnchorFinder, while not particularly slow, is still something you only want to do once.  The defaults suggested above produce lots of output, which can then be further processed quickly without reloading PDBs.  To expedite this, AnchorFinder produces two levels of output.  All residues have their data printed to a file named (pdbname).data – you can reprocess this to get data for differing window lengths, loopnesses, etc.  Windows passing the loopness and interactions filters are printed to the specified output file.

A suggested options file for AnchorFinder is available with this document.

Interpreting AnchorFinder Results
---------------------------------
(This section is minimally relevant if not benchmarking)

After you've run AnchorFinder, you'll have a fairly large pile of output: pdbname.data for all pdbs, plus goodfile.out for the better windows.

`pdbname.data` looks like this:

Rows are residues, columns are chains, data are neighbors in that chain for each residue

    residue chain   PDBdata DSSP    1       2
    1       1       2 D     L       7       0
    2       1       3 D     L       10      0
    3       1       4 D     L       14      0
    ...

The columns are residue and chain in Rosetta numbering, residue/chain in PDB numbering, DSSP value, and then N columns for the N chains in the protein.  The number in those columns is the number of cross-interface neighbors on that chain for that position.

`goodfile.out` looks like this:

    PDB pdb2vk1 window 45 loopness 5 nbrs 0 28 0 0 start 46 A pymol select pdb2vk1 and chain A and resi 46-50
    PDB pdb2vk1 window 108 loopness 5 nbrs 0 25 0 0 start 109 A pymol select pdb2vk1 and chain A and resi 109-113
    PDB pdb2vk1 window 109 loopness 5 nbrs 0 36 0 0 start 110 A pymol select pdb2vk1 and chain A and resi 110-114
    PDB pdb2vk1 window 110 loopness 5 nbrs 0 46 0 0 start 111 A pymol select pdb2vk1 and chain A and resi 111-115
    PDB pdb2vk1 window 111 loopness 5 nbrs 0 46 0 0 start 112 A pymol select pdb2vk1 and chain A and resi 112-116
    PDB pdb2vk1 window 112 loopness 5 nbrs 0 47 0 0 start 113 A pymol select pdb2vk1 and chain A and resi 113-117

Each line identifies the PDB, the window number, its loopness, its number of neighbors on each chain in the PDB (variable # of columns), the starting residue PDB numbering for the window, and a Pymol selection for the window.

Inputs and outputs for this stage from a convenience sample (PDBs 3cy?) are included with this protocol capture.

At this point, the data is yours to play with.  I searched for windows with large numbers of neighbors on only one chain using sifter.py (included), then sorted for those with the largest number of neighbors (sort -n -k1 `input`).  After that it was all manual filtering to choose structures for the benchmarks.

Choosing Loop and Anchor — Benchmarking
---------------------------------------
(This section applies to benchmarking only)

OK, so you ran AnchorFinder, looked at the results, and/or picked what protein you want to run through AnchoredDesign.  How do you choose a loop/anchor?

If you ran AnchorFinder, look at the AnchorFinder result lines that came up as good:

    92 PDB pdb1zr0 window 526 loopness 5 nbrs 0 0 92 0 start 13 D pymol select pdb1zr0 and chain D and resi 13-17

Load this PDB into pymol (1zr0.pdb) and activate the suggested selection.  You'll see that it is in a surface loop of one partner which sticks an arginine straight into its binding partner – a perfect anchor.  (This is a chosen example; not all AnchorFinder hits are this nice.)

Choosing the anchor is entirely up to human effort; here the arginine 15 is an obvious choice.

For choosing loops, I just traveled up and down the chain in both directions until I hit secondary structure, significant backbone-backbone hbonding, or the protein core.  Here I'd choose a loop of D10 to L17 – more N-terminal than that affects the core, and more C-terminal affects a sheet.

Anchor and loop file specifications are included in the release documentation and the examples here.

Note that for the included example, the PDB has been renumbered from 1.  Scripts to do this are occasionally included with Rosetta distributions and not included here.  It will be convenient to also remove waters, ligands, etc.

If you are doing benchmarking, skip to the [Running AnchorDesign](#Running 
AnchorDesign) section.

Choosing a System and Anchor — Design
-------------------------------------
(This section applies to the design case only)

In the design case, you will be choosing your proteins based on what you want designed.  Your target is forced by what targets:
* have crystal structures, and
* are related to your biological problem.

Choosing an anchor then requires:
* a cocrystal of your target with some partner from which to source the anchor.

You can run this cocrystal through AnchorFinder and let it suggest anchors to you, but for one structure you can just look at it yourself.  Look for loops on the partner that insert into the target, or do computational alanine scanning, or examine the literature for mutations that disrupt the interface.

Choosing a Scaffold, Design Positions, and Loops — Design
---------------------------------------------------------
(This section applies to the design case only)

In the design case, you will be replacing your target's partner with some new scaffold to form a mostly de novo interface.  Your scaffold must meet a few requirements:
* flexible, mutateable surface loops (for AnchoredDesign to modify)
* experimentally tractable (hey, your funeral if it's not)
* whatever other functionality you need for your desired design

The protocol was written with the fibronectin monobody scaffold in mind.

Choosing which loops are flexible is dependent on biological knowledge of the scaffold.  In fibronectin's case, many papers have been published establishing the mutability of the BC and FG loops.

Choosing which positions are designable is similarly dependent on your scaffold.  AnchoredDesign carries the assumption that the non-anchor loop positions are designable, and non-loop positions are not, but nothing in the code enforces that.  Use a resfile (documented with the manual) to specify which positions are designable.  The code will automatically prevent design of the anchor (you can turn that off).  The code will automatically prevent design of positions that are not close to either the interface or a flexible loop (you cannot turn that off), so take care in specifying designable positions on opposite faces of your protein.  Proximity is redetermined at each design opportunity so positions peripheral to the interface may not be designed regularly.

Choosing Loop Lengths and Anchor Position
-----------------------------------------
(This section applies to the design case only)

OK, so you know which scaffold to use, and which anchor, and which target.  You are ready to create your starting structure for AnchoredDesign, in which the anchor will be inserted into the scaffold, and the anchor will be aligned properly to the target, dragging the scaffold and target together.  The protocol used for this is called AnchoredPDBCreator; further details are below.

One important part of conformational space that AnchoredDesign cannot search is the space of loop lengths and anchor positions.  You may want to try, for a loop of length N, all combinations of loops of length N-3 to N+3, or even more for long loops.  As you are designing the loop to form an interface, there is no reason to believe its native length is particularly relevant.  You will have to do this searching at this stage: create starting structures for all loop lengths, run them all through AnchoredDesign, and pick off the best ones later.

Loops can be shortened directly by just deleting residues mid-loop before handing the scaffold to AnchoredPDBCreator – it can insert a 3 residue anchor into a 6 residue window, and close the gap.  Loop lengthening must be done externally.  One way to lengthen loops is to manually modify a PDB to contain enough residues in the loop (copy-and-paste a residue, renumber as necessary), then use the loop_modeling executeable's build_initial mode to close the loop.  Further instructions are included in their own folder in this packet.

A paired space is anchor placement space.  Besides choosing which anchor to use (try several), exactly where it is placed within a loop can vary.  For a loop of length 7, and an anchor of length 2, (assuming a flexible residue on each side), you have the following 4 choices:

    X = scaffold
    - = loop
    A = anchor
    X1234567X
    X-AA----X
    X--AA---X
    X---AA--X
    X----AA-X

Again, this space is not searched by AnchoredDesign and must be searched by trying all the inputs.

Using AnchorPdbCreator
----------------------
(This section applies to the design case only)

AnchoredPDBCreator is the protocol which assembles an anchor, scaffold, and target into a starting structure for AnchoredDesign.  Its code documentation is included in this packet.

Briefly, AnchoredPDBCreator takes as input 4 files:
* The target structure, as a PDB, with the partner removed
* The anchor structure, drawn from the cocrystal with the target, containing ONLY the residues being used as an anchor, as a PDB
* The scaffold structure, as a PDB, with loop residues added/deleted as desired
* A scaffold_loop specification, which declares which residues in the scaffold are flexible and where the anchor insertion should occur.

It is ABSOLUTELY VITAL to recognize that AnchoredPDBCreator does NOT produce interfaces, it only produces starting structures for AnchoredDesign.  It is entirely plausible that its structures will have the target and scaffold totally eclipsed.  This is fine, AnchoredDesign will fix it.

AnchoredPDBCreator's results should be interpreted by analyzing ONLY the closure of the anchored loop.  Use the result with the best loop geometry.  Loop geometry can be measured by examining the LoopAnalyzerMover output tagged to the end of result PDBs:

LoopAnalyzerMover: unweighted bonded terms and angles (in degrees)

    position phi_angle psi_angle omega_angle peptide_bond_C-N_distance rama_score omega_score dunbrack_score peptide_bond_score chainbreak_score
     pos phi_ang psi_ang omega_ang pbnd_dst    rama  omega_sc dbrack pbnd_sc   cbreak
      17  -106.8   175.8     178.2    1.322   0.998    0.0342   7.01   -2.68   0.0182
      18  -82.33   64.67    -178.5    1.329   0.211    0.0217   3.11   -3.42   0.0203
      19  -83.63   149.4     177.2    1.329   -1.07    0.0795      0   -3.43    0.584
      20  -75.25   171.1    -178.7    1.329  -0.264    0.0161  0.348   -3.43   0.0151
      21  -58.53  -42.95     174.6    1.329   -0.58     0.294      0   -3.43      2.7
      22  -76.02   159.9    -179.8    1.326  -0.811  0.000404   0.97   -3.45   0.0424
      23  -72.63   130.1     179.4    1.325   -1.29   0.00372   0.24   -3.46   0.0281
      24  -94.91   116.5     179.8    1.323   -1.21   0.00028  0.721   -3.45   0.0694
      25  -65.42   150.7     179.4    1.335   -1.58     0.004      0   -3.32     1.38
      26  -64.68   147.9     179.1    1.323   -1.45    0.0079   1.61   -3.32    0.211
      27  -56.44  -66.68      -180    1.329    1.34  8.08e-30   7.87   -3.43 2.37e-05
      28  -124.4  -56.48     177.6    1.329    2.08    0.0568  0.608   -3.43   0.0533
      29  -124.1   28.78    -177.7    1.264   0.341    0.0542   2.39    2.65     2.07
      30   81.57  -134.3    -176.4    1.329      20     0.126   5.06    2.65    0.128
      31  -112.9   147.2     172.7    1.318  -0.744     0.538  0.534   -3.35     1.38
    total_rama 15.9674
    total_omega 1.23676
    total_peptide_bond -38.3223
    total_chainbreak 8.70689
    total rama+omega+peptide bond+chainbreak -12.4113

    LAM_total -12.4113

In this particular example, position 29 is clearly problematic: the peptide bond distance is too short, as reported by the pbnd_dst, pbnd_sc, and cbreak columns.

You should be running AnchoredPDBCreator for at least 100 trajectories before choosing a starting structure.

Running AnchorDesign
--------------------
If you are benchmarking, the crystal structure of the complex is the appropriate input for AnchoredDesign.  If you are designing, the best result from AnchoredPDBCreator is your starting structure.

The input files for AnchoredDesign provide an example with 1zr0 for running AnchoredDesign.  It is a heterodimer so you can pretend it was AnchoredPDBCreator sourced if you want.  (You can also look in the AnchoredDesign integration test at test/integration/tests/AnchoredDesign for such an input).

* The anchor file specifies which residues form anchor.
* The PDB file is the pdb.
* loopsfile_extended is for benchmarking – the extension column is true, which tells AnchoredDesign to forget the starting loop conformation before sampling
* loopsfile_native tells AnchoredDesign NOT to forget the starting loop conformation – this is probably what you would use for design, although there is no reason you can't reject the starting loop conformation for design.  (This is cheating for benchmarking, but useful for making relaxed natives)
* options is the command line options file.  The active options are for a simple “does it run” test; parameters for a longer running “real” test are included.
* A resfile is necessary if you wish to design, and is only active in that option set.  

Interpreting AnchorDesign — Benchmarking
----------------------------------------
(This section applies to the benchmarking case only)

If you are duplicating the benchmarking results, you passed the rmsd flag. AnchoredDesign will have output a lot of RMSD values allowing you to determine the performance of the protocol against the structures you chose to benchmark.  The paper describes the score versus RMSD metrics used to determine quality (including the I_sup_bb_RMSD, ch2_CA_RMSD, and loop_CA_sup_RMSD.  The structures themselves don't really matter; you are ensuring that the low-scoring structures have low RMSD.

Interpreting AnchorDesign — Design
----------------------------------
(This section applies to the design case only)

In the design case, the other fields of the AnchoredDesign output come in to play. There are three classes of output:
* scorefunction terms
* LoopAnalyzerMover output,
* InterfaceAnaylzerMover output.

Generally, you should rank your structures according to total_score (the Rosetta scorefunction).  This tells you what Rosetta thinks is best.

Next, you use the LoopAnalyzerMover output (described above) and InterfaceAnalyzerMover output to determine which structures have flaws not caught by total_score.  Toss structures that those filters think have problems.  Pick the ones you think are best, order the DNA, and pray.  When it works great, feel free to send me kudos, citations, or money!

InterfaceAnalyzerMover
----------------------
InterfaceAnalyzerMover output looks like this:

    Residues missing H-bonds:
    Residue 	 Chain 	 Atom 
    38 	 A 	 NE2
    101 	 A 	 OE1
    248 	 A 	 O
    250 	 A 	 O
    344 	 B 	 N
    384 	 B 	 O
    477 	 B 	 O

    pymol-style selection for unstat hbond res 
    select start_5411_unsat, /start_5411//A/38+101+248+250+ + /start_5411//B/344+384+477+

    pymol-style selection for interface res 
    select start_5411_interface, /start_5411//A/31+32+33+34+35+36+37+38+39+40+41+54+56+57+59+60+61+62+64+65+66+92+95+98+99+100+101+102+103+106+194+195+224+225+226+227+228+229+230+247+248+249+250+251+252+253+265+ + /start_5411//B/314+315+316+317+318+319+320+321+322+323+324+337+339+340+342+343+344+345+347+348+349+350+375+378+379+381+382+383+384+385+386+389+473+476+477+478+479+480+481+482+488+506+507+508+509+510+511+512+513+

The first section documents where Rosetta thinks there are unsatisfied hydrogen bonds at the interface.  This code is known to be oversensitive to missing bonds, but it's better than nothing.

The next sections print PyMOL selections for interface residues for easier visualization.

InterfaceAnalyzerMover also includes columns into the scorefile:

    dSASA_int 2396.33
    dG_separated -35.3379
    dG_separated/dSASAx100 -1.47467
    delta_unsatHbonds 7
    packstat 0
    dG_cross -27.6963
    dG_cross/dSASAx100 -1.15578
    AllGly_dG -2.83564
    cen_dG -10.3844
    nres_int 96
    per_residue_energy_int -1.15006
    side1_score -361.478
    side2_score -267.353
    nres_all 520
    side1_normalized -1.25513
    side2_normalized -1.15238
    complex_normalized -1.74616
    hbond_E_fraction 0.368537

Most of these are experimental and not useful (and not part of AnchoredDesign; InterfaceAnalyzerMover has other clients).  The useful ones are dG_separated/dSASAx100, which measures the Rosetta energy of binding per unit area of SASA (scaled by a factor of 100).  This ensures you pick an interface that is energetic for its size, not large but sloppy.
Fragment Picking for CS-Rosetta
===============================

Presenting author: Dominik Gront (dgront at chem dot uw dot edu dot pl)  
Protocol name: fragment picker : CS-Rosetta style  
Brief description: The protocol substitutes CS-Rosetta application

Source code location
--------------------

* Check out the mini SVN: https://svn.rosettacommons.org/source/trunk/mini/
* Fragment picker is located in: https://svn.rosettacommons.org/source/trunk/mini/src/core/fragment/picking
* Applications are in: https://svn.rosettacommons.org/source/trunk/mini/src/apps/pilot/dgront/fragmentpicker

Running the protocol capture
----------------------------

1. Set up the path to minirosetta database
2. Set up the path to vall database
3. Run the picker:
   ```
   picker.linuxgccrelease @cs-rosetta.flags
   ```
# Temperature-Sensitive Mutation Prediciton

## Author
Chris Poultney

## Brief Description
Uses Rosetta relax and score protocols, in combination with software listed
below, to predict which mutations to a structure will result in a protein
with a temperature-sensitive phenotype.

## Related RosettaCon Talk

### Title, Authors & Lab , Year, Session and Day of talk
Prediction of temperature sensitive mutations, Chris Poultney, Bonneau Lab
Design session, 8/4/2010

### Abstract

The function of a gene product is often interrogated by means of
mutagenesis or knockout studies. However, in the case of essential
genes such perturbations would simply result in the uninformative
embryonic lethal phenotype. A solution is to design
temperature-sensitive (ts) mutants, which exhibit a mutant phenotype
only at high (or low) temperatures, providing a simple means to
"switch on" a mutation at any stage. Finding ts mutations typically
relies on generating and screening many thousands of mutations, which
is an expensive and labor-intensive process. Here we describe an
in-silico method that uses the Rosetta relax protocol and machine
learning techniques to predict a highly accurate "top 5" list of ts
mutations given a protein of interest.

## Running
### Flags
scripts/generate-scripts.sh flags (all required):
- -protein: name of input protein (must be PDB file name without .pdb extension)
- -species: species name abbreviation (used for display only, but required)
- -cutoff: % surface area accessibility cutoff for residues to mutate in input structure
- -mini\_bin: path of Rosetta bin directory on computer where runs will be executed
- -mini\_db: path of Rosetta database directory on computer where runs will be executed

scripts/predict-ts.sh flags (required):
- -protein: name of input protein (must be PDB file name without .pdb extension)

### Command Line

These command lines are not called directly, but from scripts generated by the protocol (see "Other Comments").
All files are written to and read from the ts\_mutant directory; the input\_files and output\_files directories are
ignored.
```
~/rosetta-3.0/bin/relax.linuxgccrelease -database ~/rosetta-3.0/rosetta3_database -s YBR109C-F140A.pdb -native YBR109C.pdb -nstruct 50 -relax:fast -out:file:scorefile YBR109C-F140A.sc -out:pdb_gz

~/rosetta-3.0/bin/score.linuxgccrelease -database ~/rosetta-3.0/rosetta3_database -s YBR109C-F140A_????.pdb.gz -in:file:native YBR109C.pdb -in:file:fullatom -out:file:scorefile YBR109C-F140Arescore.sc
```

### Example Overall Command Line 

There are three steps (see "Other Comments" for details): generating scripts to runs protocols,
running the scripts (usually on a cluster), and making ts predictions. All must be executed from
the ts_mutant directory:
```
scripts/generate-scripts.sh -protein YBR109C -species Scer -cutoff 10 -mini_bin ~/rosetta-3.0/rosetta3_source/bin -mini_db ~/rosetta-3.0/rosetta3_database

for a in *.sh; do qsub -d $(pwd) $a; done

scripts/predict-ts.sh -protein YBR109C
```

## Versions

Rosetta 3.0 release: https://svn.rosettacommons.org/source/branches/releases/rosetta-3.0

In order to compile Rosetta 3.0 on recent Linux distributions, the
compilation settings need to be changed to use gcc 4.3. Switch into
the rosetta-3.0/rosetta3_source directory and execute the following:
```
patch -p0 < [path_to_ts_protocol_dir]/patch/r30gcc43.patch
```
Then compile the relax and score executables, replacing the -j
parameter (number of concurrent compilation threads) as appropriate:
```
scons -j 6 bin/relax.linuxgccrelease bin/score.linuxgccrelease mode=release cxx_ver=4.3 extras=static
```

## Version for Other Codes Used

The following must all be installed and available on your PATH:

sed and awk

Probe 2.12 or better: http://kinemage.biochem.duke.edu/software/probe.php

PyMOL 1.2 or better: http://www.pymol.org/ (some Linux distros make this available as a package)

NCBI BLAST+ tools 2.2.22 or better: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
Good linux install instructions are here: http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/unix_setup.html
The "nr" database is also required, and is most easily installed using update_blastdb.pl per above instructions

Weka 3.6 or better: weka.jar is included with protocol capture in "svm" directory
See COPYING-weka for copyright/redistribution details.

LibSVM 2.8.9: libsvm.jar is included with protocol capture in "svm" directory
See COPYRIGHT-libsvm for copyright/redistribution details.


## Other Comments 

1. Running the temperature-sensitive allele prediction protocol

    All the protocol needs to run is a starting structure. The structure
    must be avaliable as a .pdb file, and consist of exactly one
    chain. Protocol scripts take as an argument the name of the protein,
    which must be the file name without the .pdb extension. For example,
    if the input file is YBR109C.pdb, the protein name is
    YBR109C. YBR109C.pdb is provided as part of the protocol for testing
    purposes.
    
    The protocol is split into three steps: creating the Rosetta run
    scripts, performing the runs, and analyzing/predicting. This is so
    that the execution stage can be run on a different computer from the
    generation and prediction stages, as some clusters do not handle
    manipulating many small files well. All scripts live in the scripts/
    directory.

2. Generating script files

    Generating script files for the Rosetta runs is done using
    scripts/generate-scripts.sh. For details, execute
    ```
    scripts/generate-scripts.sh -usage
    ``` 
    This script requires five arguments: the protein name, the species
    abbreviation, the path to the Rosetta executables, and the path to the
    Rosetta database. IMPORTANT: the Rosetta paths given must be the paths
    on the machine where the Rosetta runs will be performed! In other
    words, if you plan to generate the scripts on one computer and execute
    them on another, be sure the paths are valid for the execution
    computer.
    
    For example, to generate scripts for predictions on the provided yeast
    protein YBR109C at all positions in the native structure with
    accessible surface area of 10% or less, using the Rosetta executables
    at ~/rosetta-3.0/rosetta3_source/bin and the Rosetta database at
    ~/rosetta-3.0/rosetta3_database (Bash derivatives only):
    ```
    scripts/generate-scripts.sh -protein YBR109C -species Scer -cutoff 10 -mini_bin ~/rosetta-3.0/rosetta3_source/bin -mini_db ~/rosetta-3.0/rosetta3_database
    ```
    This will generate shell scripts for each Rosetta run to be performed:
    one for each mutation at each position in the starting structure with
    accessible surface area <= 10%, plus one for the starting structure
    itself. The script for the starting structure will be called
    YBR109C-WT.sh, and the scripts for the mutations will be
    YBR109C-aNNNb.sh, where a is the native residue, NNN is the position
    (which may be any number of digits), and b is the mutated residue.

3. Performing Rosetta runs

    Now that scripts have been generated and are ready to be
    run. Executing the scripts is system-dependent. The following command
    will queue all runs on a cluster running TORQUE:
    ```   
    for a in *.sh; do qsub -d $(pwd) $a; done
    ```

4. Analysis and prediction
    
    Each of the Rosetta runs in the previous steps generates a score
    file. These score files are analyzed and used to predict ts mutations
    by the predict-ts script, which generates two ranked lists of
    predictions, one for each of the SVM classifiers. This stage includes
    a PSI-BLAST processing step, which currently takes 5-10 minutes on a
    reasonably new computer. Running the command below will generate the
    ranked lists:
    ```
    scripts/predict-ts.sh -protein YBR109C
    ```
    The ranked lists of predictions are now available as
    YBR109C-svmlin.txt and YBR109C-svmrbf.txt
Membrane Protein-Protein Docking
================================

Author: Julia Koehler Leman (julia dot koehler1982 at gmail dot com)  
Corresponding PI: Jeffrey J. Gray (jgray at jhu dot edu)  
Last Updated: January 2015  
Rosetta Revision #58069

---

To run the full application, both steps for prepacking and actual docking needs
to be carried out. 

Prepacking Step
---------------

This script prepacks a membrane protein structure inside the membrane to prepare
it for membrane protein docking. It pull the partners apart along the axis between
the center-of-masses projected onto the membrane plane, does a side-chain packing
round, and pulls them back together.

#### Executable/Script

    Rosetta/main/source/bin/docking_prepack_protocol.macosclangrelease

#### Generating Inputs

You need:

1. A PDB file transformed into membrane coordinates and cleaned

2. A span file containing the topology information of the membrane protein

PDB:

1. Generate a PDB file where the membrane protein structure (our case 1AFO) is 
   transformed into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php).

2. Clean the PDB file by using clean_pdb.pl in the folder 
   Rosetta/tools/protein_tools/scripts/:
   $ clean_pdb.pl 1AFO_tr.pdb ignorechain

Span file:

1. Prepare the spanfile from the PDB file according to the instructions in the 
   demo spanfile_from_pdb

#### Running the Application

Adjust the paths in the command.sh script provided in this folder and run it:

    $ ./command.sh

For the paper, we used 10 output models and selected the lowest scoring model by
total Rosetta score. More models might be better, even though probably not totally
necessary since this does mainly sidechain optimization.

#### Example Outputs

Outputs are found in the output folder:

1. 1AFO_AB_0001.pdb  - final decoy

   The model should have MEM residues at the end (HETATM lines at residue number 81) 
   representing the virtual membrane residue. This needs to be removed for mpdocking in
   this revision even though for newer revisions you can leave that in.

2. score_ppk_1AFO.sc - score file

   Total scores can be >0 because of clashes. This is not worrysome because you 
   should run membrane protein docking afterwards which should get rid of the clashes.

MP Dock
-------

This script docks two partners of a membrane protein inside the membrane bilayer.
The app inherits from the regular DockingProtocol and uses the membrane framework.
Before running MPdocking, both partners need to be in a single PDB file (can be
accomplished using mpdocking_setup app), transformed into the membrane, and
prepacked. To do this, use the mpdocking_prepack app - demo available. 

A note on the score function: currently this application uses its own 
scorefunction (which is a combination of the membrane protein scorefunction 
[low-resolution from Vladimir Yarov-Yarovoy and high-resolution from Patrick Barth]
and the docking scorefunction), so PLEASE don't supply a scorefunction on the 
commandline (in fact, the app might crash if you do). For updates on this, please 
check the application documentation at

https://www.rosettacommons.org/docs/latest/Application-Documentation.html

Preliminary data gives good sampling (even below 1A RMSD to the native in many 
cases) for LOCAL docking and the scorefunction is decent but not perfect 
(more than half of the structures give funnels) with lowest energies in 
the 5A range. The scorefunction for this will be further optimized soon.

#### Executable/Script

    Rosetta/main/source/bin/mpdocking.macosclangrelease

#### Generating Inputs 

You need:

1.  A PDB file of a native structure: only for calculating meaningful RMSDs:
    here: native.pdb

    * Generate a PDB file where the membrane protein structure (our case 1AFO) 
      is transformed into PDB coordinates (z-axis is membrane normal). This can 
      be done either by downloading the transformed PDB directly from the PDBTM 
      website (http://pdbtm.enzim.hu/) or by downloading a PDB file from the 
      PDB and running it through the PPM server 
      (http://opm.phar.umich.edu/server.php).

    * Clean the PDB file by using clean_pdb.pl in the folder 
      Rosetta/tools/protein_tools/scripts/:

        $ clean_pdb.pl 1AFO_tr.pdb ignorechain

2.  A PDB file transformed into membrane coordinates, cleaned, prepacked:
    here: 1AFO_AB_0001.pdb

    You can use the native from step 1. If not, make sure both partners are in a
    single PDB file, transformed into membrane coordinates, and cleaned as in step 1.
    Run the MP prepacking protocol described in the demo MPdocking_prepack. 
    If things crash (they really shouldn't), remove the MEM residue in the 
    prepacked and in the native structure. 

3.  A span file containing the topology information of the membrane protein
    here: 1AFO_AB.span

    * Prepare the spanfile from the PDB file according to the instructions in 
      the demo spanfile_from_pdb.

Example inputs are found in the input folder.

#### Running the Application

Adjust the paths in the command.sh script provided in this folder and run it:

    $ ./command.sh

For production runs, build at least 1000 models. 

#### Example Outputs

Outputs are found in the output folder:

1. 1AFO_AB_0001_0001.pdb - final decoy with MEM residue

2. score_mpdock_1AFO.sc  - score file

   If you supplied a proper native, the rms column in the score file should be the
   RMSD to the native structure. The rms in column 3 is the ligand RMSD, i.e. the 
   RMSD over the docked partner, disregarding the fixed partner. The Irms in column
   5 is the interface RMSD.

References
----------

* Alford RF, Koehler Leman J, Weitzner BD, Duran A, Elazar A, Tilley D, Gray JJ 
  (2015): An integrated framework advancing membrane protein modeling and 
  design, PLosONE (in preparation)
Fragment Picking for CS Talos L1 Rama
=====================================

Presenting author: Dominik Gront (dgront at chem dot uw dot edu dot pl)  
Protocol name: fragment picker : CS-Rosetta style  
Brief description: The protocol substitutes CS-Rosetta application  

Source code location
--------------------

* Check out the mini SVN: https://svn.rosettacommons.org/source/trunk/mini/
* Fragment picker is located in: https://svn.rosettacommons.org/source/trunk/mini/src/core/fragment/picking
* Applications are in: https://svn.rosettacommons.org/source/trunk/mini/src/apps/pilot/dgront/fragmentpicker

Running the protocol capture
----------------------------

1. Set up the path to minirosetta database
2. Set up the path to vall database
3. Run the picker:
   ```
   picker.linuxgccrelease @cs-rosetta.flags
   ```
Membrane ΔΔG
============

Author: Rebecca F. Alford (rfalford12@gmail.com)  
Author: Julia Koehler Leman (julia.koehler1982@gmail.com)  
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)  
Last Updated: March 2015  

---

Measuring free energy changes upon mutation can inform our understanding of membrane protein stability and variation and is a step toward design. In this application, we predict ddGs by measuring the difference in Rosetta energy for the native and mutated conformation. This application uses a modified version of the all atom energy function for membrane proteins, which includes the fa_elec term and pH energy (see below). The Membrane ddG application is part of the RosettaMP Framework.

Documentation Link:  
* https://www.rosettacommons.org/docs/wiki/Membrane-ddG

Publication describing the method:  
* Alford RF, Koehler Leman J, Weitzner BD, Gray JJ (2015) An integrated 
  framework advancing membrane protein modeling and design PLosCompBio (Under 
  Review) 

## Executable/Script ##
The membrane ddG application is implemented as a python script in PyRosetta. The scripts described here can be found in this protocol capture. Developmental versions can also be found: 

1. In the Rosetta Source Code: 

        /path/to/Rosetta/source/src/python/bindings/app/membrane/predict_ddG.py

2. In the PyRosetta package: 

        /path/to/PyRosetta/app/membrane/predict_ddG.py

## Generating Inputs ##
Three inputs are required for the ddG application:  

1. PDB file for the protein structure transformed into the membrane coordinate frame.
2. Span file describing the location of trans-membrane spans
3. Residue position to mutate to (uses pose numbering)

Steps for generating these inputs are found below. A set of example inputs can 
also be found in inputs/. Here, OmpLA (PDB ID: 1qd6) is used as an example: 

1. PDB File: Generate a PDB file where the membrane protein structure is transformed 
   into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php).

   Once the PDB is downloaded from the PDBTM, clean the PDB and extract chain 
   C using teh clean_pdb.py script in Rosetta tools using the following command: 

        /path/to/Rosetta/tools/protein_tools/scripts/clean_pdb.py 1qd6_tr.pdb C

   The resulting PDB is referred to as 1qd6_tr_C.pdb from here. 

2. Span File: Generate a spanfile from the PDB structure using
   the spanfile_from_pdb application described in the MP_spanfile-from-pdb protocol
   capture in Rosetta/demos/protocol_captures/2015. An example commandline using 
   1qd6 is also provided here: 

        Rosetta/main/source/bin/spanfile_from_pdb.linuxgccrelease -in:file:s inputs/1qd6_tr_C.pdb

   For this example, this command will produce 1 output file: 
     = 1qd6_tr_C.span: Spanfile containing predicted trans-membrane spans

   Note: For this example, 1qd6 should have 12 transmembrane spans. Adjust the spanfile,
   if needed. 

## Steps for Running the protocol ##
Here, we describe the steps required to run the MP_ddG protocol. First, we describe how to 
assemble a simple PyRosetta script using the membrane framework for ddG predictions (predict_ompLA_ddG.py). Next, we describe use of a general ddG prediction application for larger scale use. 

1. Application-Specific Membrane ddG PyRosetta Protocol
   PyRosetta calculations can be adapted to use the Rosetta Membrane Framework
   with only a few additional steps. These include: 

   * Use AddMembraneMover (in protcols.membrane) to initialize the membrane framework
   * Use MembranePositionFromTopologyMover to orient the pose in the membrane based on the transmembrane spans (optional, but recommended)
   * Setup a membrane energy function
   * Proceed with normal Rosetta functions

   Here, we provide an example application-specific ddG calculation script for computing ddGs of mutation in OmpLA for comparison with experimental values in Moon & Fleming, 2011. The script can be run with no arguments by the following command: 

        ./predict_OmpLA_ddG.py 

   Step-by-step instructions on how to setup this script are provided in the predict_OmpLA_ddG.py script (in this protocol capture). 

   A single output file is created by this script: 

   * ompLA_ddG.out: Predicted ddGs for each mutation

2. Large Scale ddG predictions with the RosettaMP Framework
   Here, we describe the steps to run the MPddG protocol, incorporating both
   repacking and pH effects. Additional options are described in the documentation above. Here, we also use OmpLA as an example: 

   Here, you will need to specify the input PDB, spanfile, and residue position to 
   mutate. By default, ddGs to all canonical residues will be computed. A specific 
   ddG of mutation can be computed using the flag --mut <AA>. In this example, we specify a repack radius of 8.0A. This means all residues within 8A of the mutant position are repacked. We also specify the pH at which predictions are carried out. 

   This application can be run using the following command line: 

        ./predict_ddG.py --in_pdb --in_span inputs/1qd6_tr_C.span --res 181 --repack_radius 8.0 --include_pH true --pH_value 4.0

   Two output files (with default names) are created by this script
   * ddG.out: predicted ddGs per mutation
   * scores.sc: Breakdown of ddGs by Rosetta score term (weighted)

## Example Outputs
Outputs from the two scripts above were renamed for clarity

1. The predict_OmpLA_ddG.py script will write a list of mutations and ddG 
   values to an output file. The columns in the file are residu position, 
   mutant amino acid (1-letter code) and predicted ddG

   An example is provided in example_outputs/OmpLA_predicted_ddGds_specific.out

2. The predict_ddG.py script will create two output files: 

   * A list of mutations and their predicted ddGs. The columns in this file are
     pdb name, mutant amino acid, mutant score, native score, ddG. An example is provided in

        example_outputs/OmpLA_predicted_ddGds_general.out

   * A breakdown of each ddG by Rosetta score term (weighted). An example 
     output file is provided in

        example_outputs/OmpLA_ddG_breakdown.sc

For all scripts - running multiple times with the same output path specified will APPEND to the file and not overwrite it

## Version Info ##
Rosetta Revision #58069 
Python Version: 2.7+  
PyRosetta Release Version: March 2015  

PyRosetta can be downloaded from http://www.pyrosetta.org. Follow the instructions
provided in the README to setup your shell environment

## Additional References ##
1. Chaudhury S, Lyskov S, Gray JJ (2010) PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta.

2. Moon CP, Fleming KG (2011) Side-chain hydrophobicity scale derived from transmembrane protein folding into lipid bilayers. Proc Natl Acad Sci. 

3. Kellogg, Elizabeth H., Leaver-Fay A, and Baker D. “Role of Conformational Sampling in Computing Mutation-Induced Changes in Protein Structure and Stability.” Proteins 79, no. 3 (March 2011): 830–38. doi:10.1002/prot.22921.

4. Kilambi, KP, and Gray JJ. “Rapid Calculation of Protein pKa Values Using Rosetta.” Biophysical Journal 103, no. 3 (August 8, 2012): 587–95. doi:10.1016/j.bpj.2012.06.044.

# C++ RosettaSurface

## Authors
written by Michael Pacella (Graylab), mpacella88@gmail.com

## General Description
This demo will describe how to run the C++ version of the RosettaSurface 
algorithm.  The test case shown here will analyze the adsorption of the 
model peptide LK-alpha to the 104 surface of calcite

## Algorithm
Simultaneous optimization of protein rigid-body orientation, backbone and 
side chain conformations on a solid surface.  

## Commands
```
surface_docking<.exe> -database <path/to/database> @input/flags
```

## Input Files
- lk_alpha_calcite.pdb = input pdb with LK alpha positioned above calcite 104
- calcite.surf = file containing surface vectors specific to the calcite surface in the input pdb
- lk_alpha_3mers<9mers> = 3mer and 9mer fragment files for LK alpha
- flags = arguments for rosetta surface

## Pre-processing
1. Ensure that the input pdb is properly formatted with the protein appearing before the surface in 
the input file and belonging to a separate chain

2.  Ensure that the specified surface vectors match the input surface

3.  Ensure that parameter files corresponding to the molecules comprising the surface exist in 
the rosetta database

## Post-processing
(in the directory with output decoys)

1. GetTop.sh Ads 4
2. cd TOP4.Ads
3. PostProcessRS.sh
    
Make sure that the folder contains either only adsorbed state PDBs (and native.pdb) or solution state PDBs. 
Also make sure that all post processing scripts (found in the scripts/post_processing directory)
are present in this directory as well 

These commands will extract the top 4 adsorbed-state decoys for analysis and generate
secondary structure, protein-protein contact maps, and protein-surface contact maps

## Example output
- Ads_SecStruct.png = adsorbed state secondary structure
- Sol_SecStruct.png = solution state secondary struccture
- Ads_ContactMap = adsorbed state contact map
- Sol_ContactMap = solution state contact map
- Surface_ContactMap.png = protein-surface contact map
- Diff_ContactMap.png = difference map between contacts in the ads/sol state
- lk_alpha_docked_to_calcite_0001.pdb = output final decoy
- score.sc = scorefile
- SolState_lk_alpha_docked_to_calcite_0001.pdb = solution state pdb
- Surface_lk_alpha_docked_to_calcite_0001.pdb = adsorbed state pdb

## Limitations
This app requires a single protein positioned above a solid surface whose parameters are 
present in Rosetta.  The protein needs to appear before the surface in the pdb file
and the two need to be separate chains

Heterodimeric Antibody Design using Multistate Design
=====================================================

Author: Andrew Leaver-Fay

---

Multistate design considers the impact that a sequence has on multiple structures
(states) simultaneously to rule one sequence more favorable (fit) for a particular
purpose than another sequence.  In the case of this protocol, I'm attempting to
design a heterodimer starting with a homodimer.  Multistate design needs to favor
the binding interactions for the heterodimeric species (AB) while disfavoring
the binding interactions for the homodimeric species (AA and BB).  This objective
is encoded in an input-file fitness function.

A new multistate protein design protocol has been implemented in Rosetta.  In an outer loop,
a genetic algorithm explores sequence space; in an inner loop, fixed-sequence rotamer repacking
finds a low-energy rotamer assignments for each of arbitrarily many fixed-backbone states.
(Fixed sequence repacking can be distributed across multiple CPUs with MPI).  The energies for
the set of states are used to define a fitness for a sequence; the function that converts from
the state energies to a fitness is definable in an input file and may be arbitrarily complex.
This protocol has been applied toward the design of bispecific antibodies by redesigning the Fc
interface; the protocol sought sequences favoring the AB heterodimer while disfavoring the AA
and BB homodimers. 

Running the protocol
--------------------

#### Essential Flags

    -entity_resfile  <fname>
    -fitness_file    <fname>
    -ms::generations <int>
    -ms::pop_size    <int>
    -ms::fraction_by_recombination <float>

The entity resfile defines the number of entity elements in the design task, and defines the
amino acid and rotamer search space for each entity element.  The entity resfile file format
is simply a resfile that's proceeded by one line containing one integer, the number of elements
in the entity.

The fitness file defines the set of states which are to be optimized, and the fitness function
that determines the fitness for an entity, given the energies of each of the states after
they've been repacked using the sequence encoded in that entity.  There are six available
commands in the fitness file; STATE, STATE_VECTOR, VECTOR_VARIABLE, SCALAR_EXPRESSION,
VECTOR_EXPRESSION, ENTITY_FUNCTION, and FITNESS.  Each command must be on its own line.
The syntax for command X may be found on in the @details section preceeding the
function definition for:

    void DynamicAggregateFunction::process_X_command

in the file

    /mini/src/devel/pack_daemon/DynamicAggregateFunction.cc

#### Important Flags

    -ms::generations <int>

I have found that the number of generations should be ~ 15 x N where N
is the number of entity elements. In this example, I have 7 residues on each side of the
interface being designed, so I have 14 entity elements, and therefore run for 210 generations.

    -ms::pop_size <int>

I have found that, using ms::generations <15\*N>, the population size of 100 is good.

    -ms::fraction_by_recombination <float>

The fraction of crossover events; the fraction mutated by 
point mutations is 1 - this number.  I have found that high point mutation rates are preferable, and
typically set the fraction by recombination to 5% or lower.

#### Less commonly used flags

    -ms::numresults <int>

The number of results to output, sorted by increasing (worsening) fitness.
Typically there are many sequences near the top sequence that are not very different in sequence space
nor in fitness.  The default is to output only the pdbs relating to the entity with the best fitness.

All of the flags that control the initialization of a packer-task may be used, but their use is
discouraged in favor of specifying behavior for residues in either the entity-resfile or the secondary
resfiles.

#### Example Rosetta Command Line:

    path/to/mini/bin/mpi_msd.linuxgccrelease -entity_resfile input_files/entity.resfile -fitness_file input_files/fitness.daf  -ms::pop_size 100 -ms::generations 210  -ms::numresults 1 -no_his_his_pairE -ms::fraction_by_recombination 0.02 -database /path/to/minirosetta_database

#### Example Overall Command Line (if overall protocol is run via a script or other program):

One MSD run is insufficient for multiple reasons: 1) one design run is never enough, 2) your fitness function
probably has one or more parameters that need sweeping, and should be swept through, 3) if you're modeling negative
states, then you probably need to iteratively generate those negative states.  You probably also
need scripts to control the submission of the MSD jobs to an MPI cluster.  That will vary from cluster to cluster.

In the case of the heterodimer design, I also used two more rosetta applications:

    docking_protocol.linuxgccrelease

to re-dock the homodimeric species following their output from the multistate design protocol
with the flags

    -s <pdbname>
    -database /path/to/minirosetta_database
    -nstruct 20
    -docking:docking_local_refine 1
    -no_his_his_pairE


and

    InterfaceAnalyzer.linuxgccrelease

with the flags

    -s <pdbname>
    -database /path/to/minirosetta_database
    -jd2::no_output
    -jumpnum 1
    -overwrite
    -is_compute_hbond_unsat true
    -is_compute_packstat true
    -mute protocols.toolbox
    -no_his_his_pairE

to measure the interface energy, the buried surface area (dSASA), the "binding energy density" (interface energy / dSASA),
and the number of buried unsatisfied hydrogen bonds.

Versions
--------

* Committed to trunk at revision 36337
* InterfaceAnalyzer was used with version "32988:33373M"
* docking_protocol was used with version "34393:34892M"

Other Comments
--------------
Every state which contributed in some way to the fitness for a particular entity is output
at the conclusion of the MSD run. The state is output with the rotamer assignment computed
when the entity's fitness was first evaluated (states are not repacked prior to being output,
so if you think there is something fishy in a rotamer packing, you can look at the rotamer
assignment). For example, if you used the "vmin" function to select the state 
from a state vector with the lowest energy, then only the state with the lowest 
energy is output.

Output pdbs are named with "msd\_output\_", the rank of the source entity 
(1..numresults), the name for the state variable or the state-vector variable 
(this name comes from the fitness file) and a ".pdb".

Membrane Symmetric Protein-Protein Docking
==========================================

Author: Rebecca F. Alford (rfalford12@gmail.com)  
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)  
Last Updated: January 2015  
Rosetta Revision #58069 

---

This application assembles and docks symmetric protein complexes in the membrane
bilayer. The symmetric native complex is first refined using the MP_Relax
application. The lowest scoring native conformation (by total Rosetta score) is
then used as input to the membrane symmetric docking application, which searches
for possible conformations by reassembling and docking subunits together. 

This application combines the membrane framework, symmetry machinery, and standard
symmetric docking algorithm in Rosetta. Currently, docking of Cyclic (C) symmetries is
supported. 

Reference for this protocol capture:
* Alford RF, Koehler Leman J, Weitzner BD, Duran A, Elazar A, Tilley D, Gray JJ 
  (2015) An integrated framework advancing membrane protein modeling and design 
  PLosCompBio (in preparation) 

## Executable/Script ##

    Rosetta/main/source/bin/membrane_symdocking.linuxgccrelease

## Generating Inputs ##

Three initial input files are required for this protocol: 

1. PDB file for the native symmetric complex (all subunits)
2. Span file describing trans-membrane spans of the full complex
3. Span file describing trans-membrane spans of the asymmetric unit

Steps for generating these inputs are provided below. These inputs can also be found 
in the example_inputs/ directory

1. PDB File: Generate a PDB file where the membrane protein structure is transformed 
   into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php).

   We use the potassium channel KcsA as an example here: 
   * Download the PDB 1bl8 transformed into membrane coordinates from the PDBTM
   * Clean the PDB using 

            Rosetta/tools/protein_tools/scripts/clean_pdb.py 1bl8_tr.pdb ignorechain

     (the resulting PDB is referred as 1bl8_tr.pdb from here)

2. Full (symmetric) and asymmetric Span File: Generate a spanfile from the PDB structure using
   the spanfile_from_pdb application described in the MP_spanfile-from-pdb protocol
   capture in Rosetta/demos/protocol_captures/2015. 

        Rosetta/main/source/bin/spanfile_from_pdb.linuxgccrelease -database /path/to/Rosetta/main/database -in:file:s 1bl8_tr.pdb

   For this example, this command will produce 5 output files: 
   * `1bl8_tr.span`: Predicted trans-membrane spans for the full symmetric complex
   * `1bl8_tr<A-D>.span`: Predicted trans-membrane spans for each chain in the complex

## Steps of the protocol ##

Here, we describe the steps required to run the MP_SymDock protocol. As an example, all steps 
use a C4 Symmetric Potassium Channel (PDB ID: 1bl8) 

1. Initial Refinement: Using the native symmetric complex and full spanfile, generate 
   10 refined models using the MP_Relax protocol. This protocol described in the 
   a protocol capture in Rosetta/demos/protocol_capture/2015/MP_relax. Run the following
   commandline with the given flags file: 

        Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease -parser:protocol membrane_relax.xml @relax_flags

   The following output files will be generated: 
   * `1bl8_tr_<0001-0010>.pdb`    : 10 refined models of 1bl8
   * `relax_scores_1bl8.sc`       : Rosetta scores for each resulting models

   Examples of these outputs can be found in example_refined_models/

   Note on timing: Depending on protein size, refinement is a time consuming step. 
   Each decoy will take 0.5-1.0hrs depending on avaialble processing power. 

2. Input model selection: Use the score file from the refinement step to select the lowest scoring
   refined model by total Rosetta score. In this example, the model 1bl8_tr_0009.pdb has the lowest
   total Rosetta score and will therefore be used as input to the next step. From this point forward, 
   this model will be referred to as 1bl8_refined.pdb and is located in example_symmetry_files/. 

3. Generate inputs for symmetry: To prepare the structure for assembly and docking in the protocol, 
   a set of asymmetric inputs must be generated. These inputs describe the asymmetric unit, which will 
   later be used to re-assemble the complex based on a generated symmetry definition. A version of these
   generated input files are provided in the example_symmetry_files/ directory

   First, create the asymmetric input structure and symmetry definition file from the refined symmetric
   complex using the make_symmdef_file.pl script. An example commandline is provided below: 

        Rosetta/main/source/src/apps/public/symmetry/make_symmdef_file.pl -p 1bl8_refined.pdb -a A -i B:4 > 1bl8.c4.symm

   In this command-line: 
   * `-p` specifies the input PDB file
   * `-a` specifies the chain or chains to use as the asymmetric unit, 
   * `-i` specifies how to organize the remaining chains. In this example, 
     "B:4" means use chain B as the next subunit and arrange subunits as a C4 
     tetramer (4 subunits around the Z axis))

   This command will generate various files, only two of which are key here: 
   * `1bl8_refined_INPUT.pdb`  : PDB file containing the asymmetric subunit
   * `1bl8.c4.symm`            : Symmetry definition file describing the 
     arrangement of subunits in the complex. This file specifically describes 
     needed translations and rotations to regenerate and assemble this complex 
     from the input file. 

   Note: To generate a correct symmetry, make_symmdef_file.pl requires all chains be of equal length. Subunits
   should also be close to <0.5Å rmsd to one another. Any asymmetry may result in an incorrect symmetry definiton. 
   To check, you can visualize the example_symmetry_inputs/1bl8_refined_symm.pdb to ensure this initial setup is 
   correct. 

   Next, you will need the 1bl8_trA.span file containing trans-membrane spans for only the asymmetric unit. 
   If your asymmetric unit contains multiple chains, you may need to assemble this file yourself from the full set
   of spans. 

4. Running the symmetric docking application: Using the asymmetric unit PDB, symmetry definition file, 
   and asymmetric unit span file as inputs, you are ready to run the membrane symmetric docking application. 
   Flags, recommended settings, and commandlines are described below: 

   1. Required Options: Options (flags) needed to run this application are described below. A file with these flags,   
      symdock_flags, is also provided for the 1bl8 example. 

            flags                                     descriptions
            --------------------------------------------------------------------------------------------------
            -in:file:s <pdbfile>                      Input PDB Structure: Asymmetric input structure
            -in:file:native <pdbfile>                 Structure of native symmetric complex for RMSD calculations
            -membrane_new:setup:spanfiles <spanfile>  Spanfile describing spanning topology of asymmetric unit
            -membrane_new:scoring:hbond               Turn on membrane depth-dependent hydrogen bonding weight
            -symmetry:symmetry_definition             Symmetry definition file
            -symmetry:initialize_rigid_body_dofs      Locally sample rigid body conformations during intial complex assembly
                                                      (before docking algorithm)
            -nstruct                                  Number of structures to generate
            -packing:pack_missing_sidechains 0        Wait to pack until the membrane mode is turned on
            -docking:dock_lowres_filter 5.0 10.0      Lower van der Waals scoring criteria during centroid stage
                                                      to allow wider range of rigid body sampling;
                                                      Both numbers are highest allowed scores for a low-resolution
                                                      model to proceed to the high-resolution stage:
                                                      5.0 = highest allowed van der Waals score
                                                      10.0 = highest allowed interchain contact score

  2. Recommended # of Decoys
     * For demo run: 1
     * For production runs: 1000

  3. Command Line

            Rosetta/main/source/bin/membrane_symdocking.linuxgccrelease -database /path/to/Rosetta/main/database @symdock_flags 

## Example Outputs ##
The following outputs will be generated from the symmetric docking protocol. A 
version of these outputs are also provided in the example_outputs/ directory: 

1. `1bl8_refined_INPUT_0001.pdb`: Symmetrically docked output model from the 
   protocol

2. `score.sc`: Scorefile output by Rosetta containing memrbane and symmetry 
   scores for this model

## Additional References ##
1. DiMaio F, Leaver-Fay A, Bradley P, Baker D, André I (2011) Modeling Symmetric Macromolecular 
  Structures in Rosetta3. PLoS ONE 6: e20450. 

2. Barth P, Schonbrun J, Baker D (2007) Toward high-resolution prediction and design of 
  transmembrane helical protein structures. Proc Natl Acad Sci 104: 15682–15687. 

# PyRosetta RosettaSurface

## Authors
Written by Emily Koo and Michael Pacella (Graylab), mpacella88@gmail.com

## General Description
This demo will describe how to run the PyRosetta version of the RosettaSurface 
algorithm.  The test case shown here will analyze the adsorption of the 
biomineralization protein osteocalcin to the 100 surface of the mineral
hydroxyapatite

## Algorithm
Simultaneous optimization of protein rigid-body orientation, backbone and 
side chain conformations on a solid surface.  

## Commands
1. `cd rosetta_inputs/Osteocalcin demo`
2. `surface_docking.py@flags`

## Input Files
- 1Q8H100.pdb = input pdb file, extended osteocalcin positioned above HAp 100
- 1Q8H_disulf = disulfide bond specification file
- 1Q8H_native.pdb = native crystal structure of osteocalcin
- aa1Q8H_03_05.200_v1_3 = 3mer fragments for osteocalcin
- aa1Q8H_09_05.200_v1_3 = 9mer fragments for osteocalcin
- flags = arguments for rosetta surface


## PyRosetta RosettaSurface Setup and Pre-processing

Before running the scripts for the first time, make sure that all the following details are addressed.
 
1. When PyRosetta is installed, the install directory may be different from user to user, so adding 
the scripts directories to the PATH environment variable in .bashrc file will allow the scripts to be 
run from any directory, independent from the actual location of the scripts. The following 
instructions show how the directories are added:

    1. Create/open .bashrc in home directory
        ```
        vi ~/.bashrc
        ```
        
    2. If file exists, skip to 3. Else, add the following lines to the file:
        ```
        # .bashrc

        # Source global definitions
        if [ -f /etc/bashrc ]; then
            . /etc/bashrc
        fi

        export PATH=$PATH
        ```
        
    3. Add all directories containing scripts to end of PATH statement, delimited by a colon
        ``` 
        export PATH=$PATH:/path/to/scripts:/path/to/scripts/pdb_objects/
        ``` 
        where /path/to/scripts and /path/to/more/scripts should be modified to the correct directories.

    4. Make sure that the SetPyRosettaEnvironment.sh is sourced every session. Add 
        ```
        source /path/to/PyRosetta/SetPyRosettaEnvironment.sh
        ```
                
    5. Save and close file, then source it.
        ```
        :wq
        source ~/.bashrc
        ```

2. Make sure the python and bash scripts are given the permission to be executed. Run the following 
command in the directory containing the scripts
    ```
        chmod +x *.py *.sh ./pdb_objects/*
    ```
        
3. Gnuplot is the plotting program used to generate the plots, so the program path has to be modified 
in each Plot\* file for the plots to be generated automatically. The version used is 4.2.

4. To generate plots, simply run:
    ```
        PostProcessRS.sh
    ``` 
    Make sure that the folder contains either only adsorbed state PDBs (and native.pdb) or solution state PDBs. 
    
## Post-processing
(in the directory with output decoys)

1. GetTop.sh Ads 4
2. cd TOP4.Ads
3. PostProcessRS.sh

These commands will extract the top 4 adsorbed-state decoys for analysis and generate
secondary structure, protein-protein contact maps, and protein-surface contact maps

## Example output:
- Ads\_SecStruct.png = adsorbed state secondary structure
- Sol\_SecStruct.png = solution state secondary struccture
- Ads\_ContactMap = adsorbed state contact map
- Sol\_ContactMap = solution state contact map
- Surface\_ContactMap.png = protein-surface contact map
- Diff\_ContactMap.png = difference map between contacts in the ads/sol state

## Limitations
This app requires a single protein positioned above a solid surface whose parameters are 
present in Rosetta.  The protein needs to appear before the surface in the pdb file
and the two need to be separate chains

Visualizing Membranes in PyMOL
==============================

Author: Rebecca F. Alford (rfalford12@gmail.com)  
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)  
Last Updated: January 2015  
Rosetta Revision #58069 
PyMOL Version: 1.7.4  

---

Visualizing the position and orientation of a biomolecule in the membrane bilayer
(embedding) is an important part of analyzing membrane models from Rosetta. This visualization 
allows the user to distingish conformations with poor embeddings from conformations with 
native-like/reasonable embeddings. 

This application displays a set of two parallel planes in PyMOL, representing the position 
and geometry of the membrane bilayer. During a simulation, Rosetta extracts the membrane center, 
normal and thickness from the MEM residue, calculates the position of the planes, and sends
this information to PyMOL in real time. PyMOL then uses this information to draw two CGO plane objects
representing the membrane. 

This tool is part of a standalone application and can also be used in combination with any JD2-supported
Rosetta application. 

Publication describing the method: 

* Alford RF, Koehler Leman J, Weitzner BD, Duran A, Elazar A, Tiley D, Gray JJ 
  (2015) An integrated framework advancing membrane protein modeling and design 
  PLoS ONE (in preparation) 

Executable/Script
-----------------

    Rosetta/main/source/bin/view_membrane_protein.linuxgccrelease

"or"

Pass the `-show_simulation_in_pymol 0` flag with any Rosetta Membrane Framework application

Generating Inputs
-----------------

Two inputs are required for using the standalone visualization app: 

1. A PDB to view 
2. Span file describing the location of trans-membrane spans

Steps for generating these inputs are found below. A set of example inputs can 
also be found in example_inputs/. Here, 1c3w is used as an example: 

1. PDB File: If an output model is not already available from Rosetta, 
   generate a PDB file where the membrane protein structure is transformed 
   into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php).

2. Span File: Generate a spanfile from the PDB structure using
   the spanfile_from_pdb application described in the MP_spanfile-from-pdb protocol
   capture in Rosetta/demos/protocol_captures/2014. An example commandline using 
   1c3w is also provided here: 

        Rosetta/main/source/bin/spanfile_from_pdb.linuxgccrelease -database /path/to/db -in:file:s example_inputs/1c3w_tr.pdb

   For this example, this command will produce 1 output files: 
     = 1c3w_tr.span: Spanfile containing predicted trans-membrane spans

Steps of the protocol
---------------------

Here, we describe the steps required to run the MP_PyMOLViewer protocol. As an example, all steps 
use the PDB 1c3w: 

1.  Required Options: Options (flags) needed to run this application. A file with these flags, pymol_flags, 
    is also provided for 1c3w in this demo: 

        flags                                  descriptions
        --------------------------------------------------------------------------------------------------
        -in:file:s <pdbfile>                   Input PDB Structure: PDB file for protein structure
        -membrane_new:setup:spanfiles          Spanfile describing spanning topology of starting structure 
                                               for full symmetric structure
        -show_simulation_in_pymol 0			       Use the PyMOL viewer to visualize membrane planes for structures
        -keep_pymol_simulation_history 1       Keep pymol frames for making movies/replaying simulations (optional)

2.  Startup the PyMOL PyRosetta Session: 

    1. Open a new session of PyMOL
    2. Run the PyMOLPyrosettaServer.py script using the following command line in the pymol window:  

            run /path/to/Rosetta/main/source/src/python/bindings/PyMOLPyRosettaServer.py

    Once run, a message should appear in the PyMOL terminal window indicating the server was 
    initialized successfully. 

3.  Run Rosetta application:  

    From the regular terminal, run the standalone application or other membrane framework apps
    from the command line: 

        Rosetta/main/source/bin/view_membrane_protein.linuxgccrelease -database /path/to/db @pymol_flags

    Within ~10-20 seconds, 2 parallel planes (PyMOL object entitled membrane_planes) will appear 
    in the PyMOL session. 

    Note: In earlier versions of PyMOL, small lines will be drawn to the MEM residue as this is a HETATM
    in the PDB. These can be turned off manually by hiding the line view for the MEM residue. 

Example Outputs
---------------
The following outputs will be generated from the standalone pymol viewer application. A version of these outputs 
are also provided in the example_outputs/ directory: 

1. `1c3w.pse`          : Example pymol session file including membrane planes objects
2. `1c3w.png`          : Example image of PDB 1c3w (bacteriorhodopsin) embedded in the membrane (from session file)
3. `1c3w_tr_0001.pdb`  : Output PDB file containing membrane position information in MEM residue

Note: Rosetta will also output a score file (score.sc) which is not needed for this analysis. 

Additional References
---------------------

1. Baugh EH, Lyskov S, Weitzner BD, Gray JJ (2011) Real-Time PyMOL Visualization for Rosetta and PyRosetta. PLoS ONE 6: e21931.

2. DeLano W (n.d.) The PyMOL Manual: Compiled Graphics Objects (CGOs) and Molscript Ribbons. Available: http://pymol.sourceforge.net/newman/user/toc.html.

# Rhiju Iron Chef Recipe

## Author
Rhiju Das, rhiju@stanford.edu

## Protocol Name
Rhiju's live demo on hacking code [creating a minimized helix for arbitrary sequence]

## Brief Description

  How not to code [renamed by Charlie & Rich to: Iron Chef Rosetta 1, Live Demo with Rotating Knives]
  Rhiju Das
  Tuesday Aug. 3, 2010

  A 20-minute attempt to demystify writing C++ protocols in mini. 
  

## Running
### Example Rosetta Command Line
```
~/src/mini/bin/protein_helix_assemble.macosgccrelease  -database ~/minirosetta_database/ -seq MRGSHHHHHHGMASIEGRGSLRDLQYALQEKIEELRQRDALIDELELELDQKDELIQMLQNELDKYRSVIRP
```

Above is a CASP target that is mostly helical (obviously the His-tag in the front is not, but hey I only had 20 minutes).


## Versions

The right Rosetta version:

https://svn.rosettacommons.org/source/branches/das_lab/mini
Revision: 36561


## Other Comments: 
The above is the Das lab branch, but I think the only thing I used from it is a function clear_conformation_viewers(), which you can comment out.

I think you can compile this in trunk if you just copy the file:
```
src/apps/pilot/rhiju/protein_helix_assemble.cc
```
to your pilot directory, and then make sure scons knows about it, in the file:
```
src/pilot_apps.src.settings.all
```
Capturing Conformational States in Proteins Using Pseudocontact Shifts (PCS)
============================================================================

Author: Kala Bharath Pilla (kalabharath at gmail dot com)  
Protocol Validated by: : JKLeman (julia dot koehler1982 at gmail dot com)  
PIs: Gottfied Otting, Thomas Huber (t dot huber at anu dot edu dot au)  
Rosetta revision: 56658 from Oct 2014  
Publication describing the method:

* Pilla KB, Koehler Leman J, Otting G, Huber T (2015) Capturing Conformational 
  States in Proteins Using Sparse Paramagnetic NMR Data. PLoS Comp Bio (In 
  preparation)

---

Description
-----------

The method described here is applicable to two-domain proteins, fusion proteins and proteins that exhibit distinct conformational states.where  one-of-two domains/part of the structure is always constant and part of the structure which is variable due to induced conformational change. Experimentally Measured Pseudocontact Shifts (Bertini et al) are used as restraints in Rosetta's ab initio executables minirosetta/AbRelax protocol to direct the conformational sampling that fits the experimentally observed PCS data (Schmitz et al).

This protocol particularly uses PCS data from multiple metal centers in a manner analogous to GPS-satellites to accurately position the 3D corrdinates of the atoms, a.k.a GPS-Rosetta (Yagi et al).

Example System
--------------

To demonstrate this protocol the two domain (NS2B, NS3pro) Dengue viral protease is taken as an example system. The viral protease exhibits distinct conformational states (Open and closed) as found in crystal structures (PDB IDs: 2FOM, 3U1I). NS3pro remains unchanged, but NS2B adopts closed conformation when bound to a ligand and (PDBID 3U1I) and an open conformation in ligand free form as observed in crystal structure (PDBID: 2FOM)

The Protocol
------------

Different conformational states for NS2B can be preferentially sampled only by changing input PCS data without changing the computational procedure. The protocol is demonstrated for each of the conformational states in seperate folders, namely 'open' & 'closed'.


Requirements
------------

1. Pseudocontact shifts for protein backbone atoms from more than 2 metal centers (in '.npc' Numbat format)
2. Starting model structure in pdb format.
3. Rosetta's 9 and 3 mer fragment files.

References
----------

1. Bertini I, Luchinat C, Parigi G (2002) Magnetic susceptibility in paramagnetic NMR. Prog Nucl Magn Reson Spectrosc 40: 249–273. doi:10.1016/S0079-6565(02)00002-X.
2. Schmitz C, Vernon R, Otting G, Baker D, Huber T (2012) Protein structure determination from pseudocontact shifts using ROSETTA. J Mol Biol 416: 668–677. doi:10.1016/j.jmb.2011.12.056
3. Yagi H, Pilla KB, Maleckis A, Graham B, Huber T, et al. (2013) Three-dimensional protein fold determination from backbone amide pseudocontact shifts generated by lanthanide tags at multiple sites. Structure 21: 883–890. doi:10.1016/j.str.2013.04.001.
4. Schmitz, C., Stanton-Cook, M.J., Su, X-C., Otting, G. and Huber, T. Numbat: an interactive software tool for fitting deltaX-tensors to molecular coordinates using pseudocontact shifts. J. Biomol. NMR (2008) 41 179-189
5. de la Cruz L, Nguyen THD, Ozawa K, Shin J, Graham B, et al. (2011) Binding of low molecular weight inhibitors promotes large conformational changes in the dengue virus NS2B-NS3 protease: fold analysis by pseudocontact shifts. J Am Chem Soc 133: 19205–19215. doi:10.1021/ja208435s.
# RosettaLigand Transform

## Conformer Generation

Conformers can be generated with a number of tools, including MOE and OMEGA.  In this case, the Conformer Generation tool included as part of the BioChemical Library (BCL) suite was used.  The following command was used to generate conformers:
```
    bcl.exe molecule:ConformerGeneration -conformers /home/kothiwsk/fragment_search csd_libraries_2012_RF0.1/rotamer_libraries/Nov2013/pdb_refinedsupplemented_lib.sdf.bz2 -ensemble rosetta_inputs/ligands/all_ligands.sdf -conformation_comparer DihedralBins -temperature 1.0 -max_iterations 1000 -top_models 100 -bin_size 30.0 -scheduler PThread 8 -add_h -conformers_single_file conformers
```

You can use any conformer generation tool you have available to you for this step. Your generated conformers should be output to a single SDF file.  Every conformer must have 3D coordinates and hydrogens added.  Conformers of the same ligand should have the same name in the SDF file.  For convenience, an example conformer file is provided at `rosetta_inputs/ligands/all_ligands.sdf`.

## Params File Generation

Params files contain the parameterization information for a ligand.  Every ligand or Residue in a protein structure input into Rosetta must have a corresponding params file.  Rosetta is distributed with a script called molfile_to_params.py which generates these files. However, this script is generally cumbersome for the generation of more than a small handful of ligands. The following a protocol for generating params files for large numbers of ligands:

All the scripts needed for this process are in the tools directory in the Rosetta distribution. each of the scripts below would normally be preceded by Rosetta/tools/hts_tools, but this directory prefix has been omitted for brevity.

1. *Split ligand files*
    
    The conformers for all ligands are initially stored in a single SDF file, but molfile_to_params.py expects 1 SDF file per ligand.  sdf_split_organize.py accomplishes this task. It takes as input a single sdf file, and will split that file into multiple files, each file containing all the conformers for one ligand. Different ligands must have different names in the sdf records, and all conformers for one ligand must have the same name.  Output filenames are based on the sha1 hash of the input filename, and are placed in a directory hashed structure. Thus, a ligand with the name "Written by BCL::WriteToMDL,CHEMBL29197" will be placed in the path /41/412d1d751ff3d83acf0734a2c870faaa77c28c6c.mol.

    The script will also output a list file in the following format:

        ligand_id,filename
        string,string
        ligand_1,path/to/ligand1
        ligand_2,path/to/ligand2
        
    The list file is a mapping of protein names to sdf file paths.

    Many filesystems perform poorly if large numbers of files are stored in the same directory.  The hashed directory structure is a method for splitting the generated ligand files across 256 roughly evenly sized subdirectories, improving filesystem performance.

    The script is run as follows:

        sdf_split_organize.py rosetta_inputs/ligands/conformers.sdf split_conformers/ ligand_names.csv
    
    Be sure the split_conformers/ directory exists before running the script.  Examples of the output of this script are in example_outputs/ligand_prep/

1. *Create Projet Database*

    The ligand preparation pipeline uses an sqlite3 database for organization during the pipeline.  The database keeps track of ligand metadata and the locations of ligand files. The project database is created using the following command:

        setup_screening_project.py ligand_names.csv ligand_db.db3
    
    An example of the project database is in example_outputs/ligand_prep

1. *Append binding information to project database*

    The next step is to create a binding data file. The binding data file should be in the following format:

        ligand_id,tag,value
        string,string,float
        ligand_1,foo,1.5
        ligand_2,bar,-3.7

    The columns are defined as follows:

    * ligand_id -- ligand_id is the name of the ligand, which must match the ligand_id in the file_list.csv file created by sdf_split_organize.py.
    * tag -- The name of the protein the ligand should be docked into. If a ligand should be docked into multiple proteins, it should have multiple entries in the binding data file.  Note that this protocol makes a distinction between protein name, and file name. If you have 4 protein files: foo\_0001.pdb, foo\_0002.pdb, bar\_0001.pdb, and bar\_0002.pdb, then you have two proteins with the names foo and bar. The scripts expect that the protein PDB files begin with the protein name.
    * value -- The activity of the ligand. If you are doing a benchmarking study and know the activity of your ligand, you should enter it here. If you are not doing a benchmarking study, or if ligand activity is not relevant to your study, value can be set to 1.0 (or anything else). This field is currently only used in a few specific Rosetta protocols that are in the experimental stages, and is typically ignored, so it is safe to set arbitrarily in almost every case.
    
    An example input file is provided. you can insert it into the project database with the following command:

        add_activity_tags_to_database.py ligand_db.db3 rosetta_inputs/ligand_activities.csv

1. *Generate Params Files*

    The next step is to generate params files. make_params.py is a script which wraps around molfile_to_params.py and generates params files in an automated fashion.  Params files will be given random names that do not conflict with existing Rosetta residue names (no ligands will be named ALA, for example).  This script routinely results in warnings from molfile_to_params.py, these warnings are not cause for concern. Occasionally, molfile_to_params.py is unable to properly process an sdf file, if this happens, the ligand will be skipped.  In order to run make_params.py you need to specify the path to a copy of molfile_to_params.py, as well as the path to the Rosetta database.

    make_params.py should be run like this:
    
        make_params.py -j 2 --database Rosetta/main/database --path_to_params Rosetta/main/source/src/python/apps/public/molfile_to_params.py ligand_db.db3 params/

    In the command line above, the -j option indicates the number of CPU cores which should be used when generating params files. If you are using a multiple core machine, setting -j equal to the number of available cpu cores. Be sure that the params/ directory exists before running the script.

    The script will create a directory params/ containing all params files, pdb files and conformer files.
    
    An example of the output params/ directory is found in example_outputs/ligand_prep

1. *Create job files*
    
    Because of the memory usage limitations of Rosetta, it is necessary to split the screen up into multiple jobs.  The optimal size of each job will depend on the following factors:

    * The amount of memory available per CPU
    * The number of CPUs being used
    * The number of atoms in each ligand
    * The number of conformers of each ligand
    * The number of protein residues involved in the binding site.

    Because of the number of factors that affect RosettaLigand memory usage, it is usually necessary to determine the optimal job size manually. Jobs should be small enough to fit into available memory.

    To make this process easier, the make_evenly_grouped_jobs.py script will attempt to group your protein-ligand docking problem into a set of jobs that are sized as evenly possible. The script is run like this:

        make_evenly_grouped_jobs.py --create_native_commands rosetta_inputs/proteins --n_chunks 1 --max_per_job 1000 params rosetta_inputs/proteins job

    If the script was run as written above, it would use param files from the directory param_dir/, and structure files from the directory structure_dir/. It would attempt to split the available protein-ligand docking jobs into 10 evenly grouped job files (--n_chunks). The script will attempt to keep all the docking jobs involving one protein system in one job file. However, if the number of jobs in a group exceeds 1000, the jobs involving that protein system will be split across multiple files (\texttt{--max_per_job}). The script will output the 10 job files with the given prefix, so in the command above, you would get files with names like "output_prefix_01.js". The script will output to the screen the total number of jobs in each file. All the numbers should be relatively similar.  If a job file at the beginning of the list is much larger than the others, it is a sign that you should reduce the value passed to --max_per_job. If the sizes of all jobs are larger than you want, increase --n_chunks.
    
    Additionally, the script will take the default ligand positions from the ligand pdb files, and the protein files from the rosetta_inputs/proteins directory, and designate these as the "native" pose of the protein-ligand complex. This feature will allow Rosetta to compute ligand RMSDs automatically, and was used in the benchmarking studies described in the manuscript. 
    
    An example job file produced using this script is found in example_outputs/ligand_prep

## Docking

After following the procedure above to prepare your ligands, you are ready to dock the ligands.  The screening job file produced in the previous step contains the paths to the input proteins and ligands and the paths to the necessary params files.  In this example, the ligand pdbs are already positioned in the ligand binding site. 

RosettaLigand protocols are built in the RosettaScripts framework, a modular architecture for creating RosettaLigand protocols. The rosetta_inputs/xml directory contains all of the rosetta protocols were tested in the manuscript, and any of these xml files can be used with the docking commands described below.  See the comments in the XML files for details.

The Rosetta ligand docking command should be run as follows:

    rosetta_scripts.default.linuxgccrelease @rosetta_inputs/flags.txt -in:file:screening_job_file rosetta_inputs/job_01.js -parser:protocol rosetta_inputs/tr_repack.xml -out:file:silent results.out

rosetta_inputs/flags.txt contains flags that are always the same regardless of the input file.

This command will dock every protein-ligand binding pair and place the output in the specified silent file.  In the benchmarking case described in the manual, 2000 models were made for each protein-ligand binding pair.  However, in a practical application 200 models would be appropriate.

## Analysis

### Practical analysis

If this protocol is being used for an application project in which the correct ligand binding position is not known, the lowest scoring model for each protein-ligand binding pair should be selected.  From that point, we recommend filtering by protein-ligand interface score (interface_delta_X), as well as the packstat score\citep{Sheffler:2009bd} which can be computed through the InterfaceAnalyzer mover.  The cutoffs for these filtering steps should depend on the range of scores present, and the number of compounds it is possible to test.

After filtering, the selected compounds should be visually inspected.  If a crystal structure exists with a known binding pose, the predicted binding poses of the unknown compounds should be compared.  Additionally, the overall binding poses of the filtered compounds should be inspected to assess whether or not they make chemical sense.  While this is a qualitative process, human intuition has proven a valuable aid in the drug design process\citep{Voet:2014de}. 

### Benchmarking analysis

Statistical analysis of the benchmarking study provided in this paper was performed using Python. analysis.ipynb is an ipython Notebook (http://ipython.org/notebook.html) containing the code necessary to reproduce these figures, as well as comments and description of that code.  See the iPython documentation for installation and usage instructions. 

# Limitations and Caveats
#RosettaScripts

Directories contain scripts from PLoS paper.
Flxbb design with RosettaScripts
Cycles of fixbb and relax are performed.

Rosetta-3.2 is a stable release of Rosetta
rosetta-41240 is a revision under development.
#Run this example by typing:
$ROSETTA_BIN/rosetta_scripts.linuxgccrelease @flags -database $ROSETTA_DATABASE -parser:protocol <rosetta_xml_file>

# Replace ROSETTA_BIN and ROSETTA_DATABASE with your rosetta directories.
# Replace linuxgccrelease with the extension you are using.
# Replace <rosetta_xml_file> with either "rosetta-3.2.xml" or "rosetta-41240.xml"
# The XML file should be matched with the revision of Rosetta you are using.
# Other versions may require modifications to the XML file 
To run:

1. Change executable and database paths in commandline

2. ./commandline enzdes_test.pdb enzdes.xml

3. output files are score.sc and enzdes_test_0001.pdb
the run.sh script calls teh rosetta executable and the flag file. You need to specify your database in here.
the flag file contains the arguments rosetta needs to know about.
Foldit Standalone Protocol Capture
==================================

Authors: Seth Cooper, David Kim, Firas Khatib, Adrien Treuille (PI), Janos 
Barbero, Jeehyung Lee, Michael Beenen, Andrew Leaver-Fay, Matthew Smith 
(Speaker), David Baker (PI), Zoran Popovic (PI)

Title of Talk: Standalone Foldit: A New and Powerful Interface to Rosetta

Date of Talk: Wednesday, August 4, 2010

Session: Interfaces to Rosetta

---

This protocol describes how to obtain the binaries and code for Foldit 
Standalone.  It also includes the example files used to demonstrate the 
software as well as the scripts and wrappers used to make the scripts.

Protein design with Rosetta typically involves transitioning between
command line tools (Rosetta) to generate designs, and graphical tools
(PyMol) to evaluate designs.  By tightening the feedback loop between
design and evaluation, we may be able to improve our design
methodology.  Standalone Foldit allows the protein designer to make
mutations; minimize ligands, sidechains, and the protein back bone;
and remodel loops, all with a sophisticated scripting interface,
advanced undo functionality, and the Rosetta energy function.

Additionally, Foldit can be used to evaluate and modify protein, DNA,
and RNA structures and any mixture of these.  Standalone Foldit is
under active development and has been used to design active enzymes.
This new interface to Rosetta is a great balance between usability,
interactivity, and power.

Running the protocol
--------------------

FolditStandalone is not a command line protocol, it is a graphical interface to 
Rosetta.

Versions
--------

Standalone Foldit demo was done with svn 35710, newer revisions should work 
similarly.

References
----------

Cooper S, Khatib F, Treuille A, Barbero J, Lee J, Beenen M, Leaver-Fay A, Baker 
D, Popovic Z, Foldit players. Predicting protein structures with a multiplayer 
online game. Nature 2010;466(7307):756-760.

Other comments
--------------

To use Lua script for reverting back to native sequence:

Place the attached scripts somewhere in your path (e.g. your scripts 
directory), and make all four executable by typing, in a terminal:

    chmod u+x lua_to_orig revert_to_orig.lua pos_id_lua.py pdb_fasta.pl

To generate the reversion script:

    lua_to_orig scaffold.pdb > /absolute/lua_script_location/lua_script_name

Then, at the lua prompt in Foldit, type:

    dofile("/absolute/lua_script_location/lua_script_name")
# Four Small Puzzles

## Author
Rhiju Das, rhiju@stanford.edu

## Brief Description

  Follow up to talk:

  An enumerative ansatz for high resolution RNA and Protein Folding
  Rhiju Das, Parin Sripakdeevong, Das Lab, Stanford Biochemistry
  Tuesday Aug. 3, 2010
  
  associated with paper: "Four Small Puzzles that Rosetta Doesn't Solve",
  submitted in Jan. 2011 to RosettaCon 2010 Special Collection of 
  PLoS One

## Abstract

A complete macromolecule modeling package must be able to solve the simplest structure prediction problems. Despite recent successes in high resolution structure modeling and design, the Rosetta software suite fares poorly on deceptively small protein and RNA puzzles, some as small as four residues. To illustrate these problems, this manuscript presents extensive Rosetta results for four well-defined test cases: the 20-residue mini-protein Trp cage, an even smaller disulfide-stabilized conotoxin, the reactive loop of a serine protease inhibitor, and a UUCG RNA tetraloop. In contrast to previous Rosetta studies, several lines of evidence indicate that conformational sampling is not the major bottleneck in modeling these small systems. Instead, approximations and omissions in the Rosetta all-atom energy function currently preclude discriminating experimentally observed conformations from de novo models at atomic resolution. These molecular “puzzles” should serve as useful model systems for developers wishing to make foundational improvements to this powerful modeling suite.


## Running
# Example Rosetta Command Line

 A full set of READMEs is available in input_files. These are directly copied from my run directories, so they better work.

## Example Overall Command Line (if overall protocol is run via a script or other program)

 A few runs involve 'StepWise Assembly' and require a master python script to setup the job, and another master python script to queue up the resulting computation, which is described as a directed acyclic graph. A full set of README_SETUPs and README_SUBs is available in input_files. These are directly copied from my run directories, so they better work too.  The runs work on an LSF queuing system, as is available on Stanford's BioX2 cluster;  they should also run fine via Condor's DAGMAN, though that queuing system has a lot of latency. [We have also written scripts to carry out the calculation on condor and torque clusters, and are almost finished with versions that use Amazon's EC2/S3. It usually takes a day or so to rewrite what we have for arbitrary systems. Our next step may be to write a version for Amazon's ElasticMapReduce, but it gets complicated since we actually have a series of map/reduce steps, not just one. ]


## Versions

1. The right Rosetta version:
```
https://svn.rosettacommons.org/source/branches/das_lab/mini
https://svn.rosettacommons.org/source/branches/das_lab/minirosetta_database
Revision: 36561

[This branched off trunk in the winter of 2009.]
```
2. Several scripts to generate a loop modeling job are in:
```
 https://svn.rosettacommons.org/source/workspaces/rhiju/python
 Revision: 40199

[the scripts  include:
  grind_dagman.py
 stepwise_post_process_cluster.py
 stepwise_post_process_combine_and_filter_outfiles.py
 stepwise_pre_process_setup_dirs.py
 extract_lowscore_decoys.py

 and various helper python scripts...
]
```

3. Finally, the queuing scripts are in:
```
 https://svn.rosettacommons.org/source/branches/das_lab/SWA_dagman_python 
 Revision: 40199
```
## References 

The protein and RNA de novo modeling make use of published protocols, e.g.:

i. Das R, Qian B, Raman S, Vernon R, Thompson J, et al. (2007) Structure prediction for CASP7 targets using extensive all-atom refinement with Rosetta@home. Proteins 69: 118-128.

ii. Mandell DJ, Coutsias EA, Kortemme T (2009) Sub-angstrom accuracy in protein loop reconstruction by robotics-inspired conformational sampling. Nat Methods 6: 551-552.

iii. Das R, Karanicolas J, Baker D (2010) Atomic accuracy in predicting and designing noncanonical RNA structure. Nat Methods 7: 291-294.

There are two manuscripts in preparation on StepWise Assembly for proteins and RNA.

## Other Comments
AbinitioRelax.macosgccrelease -database ~/minirosetta_database  -fasta 2jof.fasta -native 2jof.pdb -frag3  aat000_03_05.200_v1_3.txt  -frag9 aat000_09_05.200_v1_3.txt -out:file:silent 2jof_abrelax.out -out:file:silent_struct_type binary -abinitio:relax -nstruct 1 -ex1 -ex2 -extrachi_cutoff 0

relax.macosgccrelease -s idealize_2jof.pdb -out:file:silent 2jof_nativerelax.out -out:file:silent_struct_type binary -database ~/minirosetta_database  -frag3  aat000_03_05.200_v1_3.txt  -frag9 aat000_09_05.200_v1_3.txt -native 2jof.pdb -nstruct 1  -ex1 -ex2 -extrachi_cutoff 0


# to minimize entire structure:
stepwise_protein_test.linuxgccrelease  -s1 2ci2.pdb  -global_optimize -fixed_res `seq 1 62` -input_res1 `seq 1 62` -fixed_res `seq 1 62` -database ~/minirosetta_database -out:file:silent 2ci2_min.out -fasta 2ci2.fasta

# to carry out KIC loop modeling
loopmodel.macosgccrelease -database ~/minirosetta_database -loops:remodel perturb_kic -loops:refine refine_kic -loops:input_pdb 2ci2_min.pdb -in:file:native 2ci2.pdb -loops:loop_file 2ci2_35_45.loop -loops:max_kic_build_attempts 10000 -in:file:fullatom -out:file:fullatom -out:prefix 2ci2 -out:pdb -ex1 -ex2 -extrachi_cutoff 0 -out:nstruct 1 -out:file:silent_struct_type binary  -out:file:silent 2ci2_kic_loop35_45_from_2ci2_min.out  > kic.log
rna_denovo.macosgccrelease -database  ~/minirosetta_database/ -fasta gcuucggc.fasta -nstruct 1 -out:file:silent gcuucggc.out -minimize_rna -cycles 5000 -mute all -native gcuucggc_RNA.pdb  > farfar.log


rna_denovo.macosgccrelease -database  ~/minirosetta_database/ -fasta gcuucggc.fasta -nstruct 1 -out:file:silent gcuucggc_NATIVE.out -minimize_rna -cycles 5000 -mute all -native gcuucggc_RNA.pdb  -vall_torsions 1f7y.torsions > farfar_native.log

AbinitioRelax.macosgccrelease -database ~/minirosetta_database -fasta 1not_.fasta -frag3  aa1not_03_05.200_v1_3  -frag9 aa1not_09_05.200_v1_3  -out:file:silent 1not_abrelax_CST_increase_cycles.out -out:file:silent_struct_type binary  -nstruct 1  -cst_file 1not_native_disulf_CEN.cst  -abinitio:relax  -cst_fa_file 1not_native_disulf.cst -native 1not.pdb -increase_cycles 10  -score:weights score12.wts  -ex1 -ex2 -extrachi_cutoff 0 > abrelax.log

relax.macosgccrelease -database ~/minirosetta_database -s idealize_1not.pdb -fasta 1not_.fasta -frag3  aa1not_03_05.200_v1_3  -frag9 aa1not_09_05.200_v1_3  -out:file:silent 1not_native_relax.out -out:file:silent_struct_type binary  -nstruct 1    -abinitio:relax  -cst_fa_file 1not_native_disulf.cst -native 1not.pdb -increase_cycles 10 -score:weights score12.wts  -ex1 -ex2 -extrachi_cutoff 0 > native_relax.log


Peptide Backbone and Sequence Design
====================================

Author: Deanne Sammond
RosettaCon Talk:
* Computational design of a new protein-protein interface between Gi1 and a 
  redesigned RGS14 GoLoco, Deanne Sammond, Dustin Bosch, Glenn Butterfoss, 
  Mischa Machius, David Siderovski, Brian Kuhlman, Kuhlman lab, 2010, Session 3 
  on Wednesday August 4th

---

Our project is the computational design of a new high-affinity protein-protein 
interface.  Our model system is an x-ray crystal structure of  Gi1 bound to the 
GoLoco domain from the RGS14 protein.  RGS14 GoLoco spans two domains of Gi1, 
with the C-terminal random coil region binding to the all-helical domain of 
Gi1.  We removed this C-terminal portion of GoLoco, replacing the random coil 
with a de novo designed alpha helix.  The redesigned GoLoco binds to Gi1 with a 
dissociation constant of 810nM, the correct binding of the newly designed 
GoLoco was confirmed using disruptive mutations at the Gi1:GoLoco interface, 
and the correctness of the computational design was assessed with by x-ray 
crystallography.

This protocol builds (or extends) a backbone for a peptide bound to a target 
protein, then designs a low-energy sequence. 

Running the protocol
--------------------

#### Important flags:

    -ex1, -ex2, -exOH, -extrachi_cutoff 1 all seem to be very important for the sequence design run.

#### Example Rosetta Command Line:

    rosetta.mactel aa input_pdb _ -s g000.pdb -loops
    rosetta.mactel -design -l list_of_pdbs -tail -begin 342 -end 351 -chain_ -series bb -protein g000 -resfile g000_resfile -ex1 -ex2 -extrachi_cutoff 1 -exOH -no_his_his_pairE -tight_hb -try_both_his_tautomers 

### How to generate tail designs:

This protocol uses 2 separate rosetta runs — one is centroid mode to build 
backbone coordinates and the other is a design run to find a low-energy 
sequence — and 2 additional scripts.  Step-by-step instructions are below:

1. Generate fragments - you can do this using the [[Robetta 
   server|http://robetta.bakerlab.org/fragmentsubmit.jsp]]. I named my fragment 
   files `aag00003_04.200_v1_3` and `aag00009_04.200_v1_3`.

2. Make starting structure using createTemplate.pl and the g000.zones file.  I 
   also included a fasta file (see g000\_.fasta) because I felt that setting the 
   sequence improved the quality of the centroid models more effectively than 
   using constraints.  For example:

        createTemplate.pl -zonesfile g000.zones -fastafile g000_.fasta -parentpdb gpep_nat.pdb -outpdb g000.pdb

   The resulting file should look something like g000.pdb.  The side-chains are 
   removed from all sequence positions, and the region that will be redesigned 
   is removed.  NOTE: The input file has the sequence positions re-numbered in 
   "Rosetta numbers", so the original pdb starts with sequence position 30, 
   contains a "TER" between chain A and B, and chain B starts with sequence 
   position 496.  But my modified pdb (gpep_nat.pdb) starts with position 1, 
   with no "TER" between chains A and B, and the numbering is sequential 
   through B.

   input: g000.zones, g000_.fasta, gpep_nat.pdb  
   output: g000.pdb

3. Making centroid models

   * Copy fragments to working directory

   * Make a loop file to specity what residues can move.  See g000.loops as an 
     example.  In the loop file, the 1st numer is the # of positions to be 
     built, second number is sequence position to start design, with the final 
     number being the sequence position where the design will end. 

   * Make constraint file to make the designed region fall into the desired 
     location OR direct the redesigned peptide toward the desired orientation 
     using the starting sequence.  I did the latter.  In other words, part 2) 
     above can generate a pdb file with all alanines in the region that will 
     be designed OR the sequence can be specified with a fasta file (see 
     above) so that big hydrophobics fall in buried regions and hydrophilics 
     fall in solvent exposed regions, etc.  An example of a constraint file is 
     g000_.cst.

   * Copy over starting structure g000.pdb from step 2.

   * Run command line:

            rosetta.mactel aa input_pdb _ -s g000.pdb -loops 

   input: g000.loops, g000.cst (I didn't use constraints)  
   output example: aag000_0001.pdb

4. Merge centroid designs with fullatom starting structure (in this case 
   gpep1.pdb) using merge_pdb.csh like this:

        merge_pdb.csh gpep1_nat.pdb [list of pdbfiles].

   You will need to edit merge_pdb.csh if you want to change which residues are 
   being merged.  For example, if you start the design at sequence position 
   342, like we do here, check your gpep_nat.pdb file (original all-atom pdb 
   file) for the line # for the last atom in sequence position 341 and put this 
   # in after "head -", then check your centroid files for the first atom at 
   position 342 and put this # after "tail -".  This step is so that you don't 
   have to repack all of the gpep1.pdb positions during the fullatom 
   simulations.  (NOTE: When building with centroid mode, we don't use a TER in 
   between chains A and B.  The TER needs to be added back in.  Another merge 
   file can be used for this.)

   input: gpep_nat.pdb, list_of_pdbs (example of pdbs in list - aag000_0001.pdb)  
   output example: aag000_0001.m.pdb~ and with the TER added, aag000_0001.m.pdb

5. Making Fullatom models:

        rosetta.mactel -design -l list_of_pdbs (the *.m.pdb merged files from step 4 above WITH a TER added between chains A and B) -tail -begin 342 -end 351 -chain_ -series bb -protein g000 -resfile tail.resfile -ex1 -ex2 -extrachi_cutoff 1 -exOH -no_his_his_pairE -tight_hb -try_both_his_tautomers -linmem_ig 10 -output_hbond_info -decoystats -group_uns 

   * `-linmem_ig 10` is optional.  I used it because I was running on a 
     BlueGene and each node had very limited memory.  `-output_hbond_info`, 
     `-decoystats` and `-group_uns` are also optional.  I used those so the 
     output pdbs could be used with Ron Jacak's h-bond pymol plugin.

   * Interesting "feature" that seems to have appeared in this version - this 
     call to Rosetta is looking for fragment files named aag000l03_05.200_v1_3 
     and aag000l09_05.200_v1_3, whereas the fragment files I used when building 
     centroid models (3) were named aag00003_04.200_v1_3 and 
     aag00009_04.200_v1_3.  I just renamed the fragment files so I could get 
     this done quickly.

   input: g000_resfile  
   output example: aag000_0001.m_0001.pdb

See example of resfile (g000_resfile).  In this file the sequence positions of 
the designed region are allowed to vary, and any neighboring sequence positions 
on the target protein are allowed to relax.

Rosetta Version
---------------
SVN Revision: 29304  
https://svn.rosettacommons.org/source/trunk/rosetta++

Next Generation KIC
===================

Setting up the demo
-------------------

Before running the demo, make sure to change the following variables in your 
local environment:

    PATH_TO_EXE                                         # path to directory with Rosetta binaries
    ROSETTA_BINARY                                      # extension of Rosetta binaries, e.g. linuxgccrelease
    PATH_TO_DB                                          # path to Rosetta database

Running the demo
----------------

IO flags:

    -s 1a8d_MinPacked.pdb                               # The starting structure -- must have residues for the segment to be remodeled, but these don't need to have meaningful coordinates
    -loops:loop_file 1a8d.loop                          # definition of the loop to be remodeled, file format description at http://www.rosettacommons.org/manuals/archive/rosetta3.4_user_guide/d1/d49/loopmodeling.html
    -out:pdb_gz                                         # compress output structures
    -in:file:native 1a8d_MinPacked.pdb                  # native or reference structure for RMSD calculation -- if this flag is not set, an RMSD of 0 will be reported

Number of structures to produce (for demo):

    -nstruct 1                                          # number of structures to produce 

Number of structures to produce (for production runs):

    -nstruct 500                                        # or 500 independent simulations, depending on the cluster setup.

General kinematic closure loop modeling flags:

    -in:file:fullatom
    -loops:remodel perturb_kic
    -loops:refine refine_kic
    -run:test_cycles                                    # fast execution, uncomment for testing purposes only
    -loops:outer_cycles 5                               # for production runs
    -kic_bump_overlap_factor 0.36                       # reduces the threshold for permitted clashes in initial loop closures (default value is 0.4)
    -legacy_kic false                                   # remove a slight bias in the initial KIC implementation towards sampling the C-termial part of the remodeled segment more frequently
    -kic_min_after_repack true                          # minimize rotamers after repacking (happens every -loops:repack_period iterations after loop closure, default 20)
    -corrections:score:use_bicubic_interpolation false  # do not use bicubic interpolation, it has been observed to adversely affect the fraction of sub-Angstrom conformations on the 12-residue benchmark set

Next-generation KIC flags as described in Stein & Kortemme, PLoS ONE, 2013:

    -loops:kic_rama2b
    -loops:kic_omega_sampling
    -allow_omega_move true
    -loops:ramp_fa_rep
    -loops:ramp_rama

Packing flags:

    -ex1
    -ex2
    -extrachi_cutoff 0

Example Rosetta Command Line:

    $PATH_TO_EXE/loopmodel.$ROSETTA_BINARY -database $PATH_TO_DB @flags

Overall protocol execution (demo):

1.  `scripts/pre_min_pack.py <list_of_native_structures> <output_keyword>` (prepacking step)

    This will create a pre-min-packed structure to start the simulation from, 
    as well as a log of the run. The input is a list of starting structures to 
    be processed, as well as a keyword for the output structure and log names. 
    Note that, while standard repacking usually takes a few minutes, with the 
    -min_pack option it can take an hour or more, depending on the size of the 
    structures. This option was developed by Andrew Leaver-Fay and will be 
    discussed in a separate publication (Leaver-Fay et al., Meth Enzym, 2013).


2.  `scripts/submit_NGK.py <input_list> <output_keyword>` (loop remodeling step)

    This will generate structures in which the selected segment is remodeled (1 
    for the demo, use at least 500 for real life problems, more for longer 
    segments). The input is a list of (pre-packed) structures and a keyword for 
    output naming. The Rosetta syntax for .loop files is explained at 
    http://www.rosettacommons.org/manuals/archive/rosetta3.4_user_guide/d1/d49/loopmodeling.html

    Note that if a native or reference structure is provided, the backbone RMSD 
    of the remodeled segment (after superimposition of the fixed parts of the 
    structure) will be reported by loopmodel.release and can thus be parsed 
    from the standard output, if that is redirected into a file. The script 
    will perform redirection into a .log file which later is gzipped. 

    If submitted via qsub on an SGE cluster system, the script will distribute 
    the number of simulations specified with the -t flag across all structures 
    in the list, i.e., to generate 500 structures each for a list of 10 input 
    structures, use -t 1-5000.

    An NGK trajectory (generating one model) typically take 15-30min for 
    12-residue lops.

3.  `scripts/parse_RMSDs.pl <directory with log files>` (parsing step)

    This script will extract the scores and RMSDs for each final model from the 
    respective .log files, and optionally generate Rosetta-energy-vs-RMSD 
    plots. Note that relevant RMSDs will only be reported if a native or 
    reference structure is provided to the NGK run via the -in:file:native 
    flag. If no reference structure is available, users should consider 
    clustering of their results to identify the most commonly sampled 
    conformations (see Figure S1 of Stein & Kortemme, PLoS ONE, 2013).

Authors
-------
* Amelie Stein
* Tanja Kortemme

Rosetta Version
---------------

SVN revision 51851 (Dec 2012)

References
----------

* Stein A, Kortemme T (2013) Improvements to robotics-inspired conformational 
  sampling in Rosetta. PLoS ONE (submitted)

Membrane Relax with Menv_smooth Term
====================================

Author: Vladimir Yarov-Yarovoy

---

Modified version of membrane relax protocol that uses centroid based 
Menv_smooth term (membrane specific residue environment) derived from updated 
database of membrane protein structures.

Running the protocol capture
----------------------------

Standard relax protocol flags:

    -relax:fast
    -relax:constrain_relax_to_start_coords
    -relax:default_repeats 2
    -in:file:fullatom
    -in:file:native input/BRD4.pdb
    -s input/BRD4.pdb
    -nstruct 1
    -out:file:silent_struct_type binary
    -out:file:silent Menv_smooth-BRD4_frlx.out
    -mute core.util.prof
    -mute core.io.database
    -run:no_prof_info_in_silentout
    -out:file:scorefile Menv_smooth_score.sc

Transmembrane regions prediction input:

    -in:file:spanfile input/BRD4.span

Menv_smooth-membrane-relax protocol specific weights:

    -score:weights membrane_highres_Menv_smooth.wts

Membrane normal and membrane center search options:

    -membrane:normal_cycles 40
    -membrane:normal_mag 15
    -membrane:center_mag 2

Example Rosetta Command Line:

    ./bin/relax.linuxiccrelease \
        -relax:fast \
        -relax:constrain_relax_to_start_coords \
        -relax:default_repeats 2 \
        -in:file:fullatom \
        -in:file:native input/BRD4.pdb \
        -in:file:spanfile input/BRD4.span \
        -s input/BRD4.pdb \
        -score:weights membrane_highres_Menv_smooth.wts \
        -membrane:normal_cycles 40 \
        -membrane:normal_mag 15 \
        -membrane:center_mag 2 \
        -nstruct 1 \
        -out:file:silent_struct_type binary \
        -out:file:silent Menv_smooth-BRD4_frlx.out \
        -mute core.util.prof \
        -mute core.io.database \
        -run:no_prof_info_in_silentout \
        -out:file:scorefile Menv_smooth_score.sc \

Version
-------
If checked into trunk, SVN revision number: 39837

Other Comments
--------------

Membrane relax protocol specific input file (BRD4.span in the above example) 
generated using OCTOPUS transmembrane regions prediction server 
(http://octopus.cbr.su.se/) and 
Menv_smooth-membrane-relax/script/octopus2span.pl script as follows:

    Menv_smooth-membrane-relax/script/octopus2span.pl BRD4.octopus

where BRD4.octopus is OCTOPUS topology file gererated by the server.

Sample OCTOPUS topology file:

    ##############################################################################
    OCTOPUS result file
    Generated from http://octopus.cbr.su.se/ at 2008-09-18 21:06:32
    Total request time: 6.69 seconds.
    ##############################################################################


    Sequence name: BRD4
    Sequence length: 123 aa.
    Sequence:
    PIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFV
    WWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGA
    GIV

    OCTOPUS predicted topology:
    oooooMMMMMMMMMMMMMMMMMMMMMiiiiMMMMMMMMMMMMMMMMMMMMMooooooMMM
    MMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMM
    ooo

Sample span file generated by Menv_smooth-membrane-relax/script/octopus2span.pl:

    TM region prediction for BRD4 predicted using OCTOPUS
    4 123
    antiparallel
    n2c
       6    26     6    26
      31    51    31    51
      58    78    58    78
      97   117    97   117

* 1st line is comment line.
* 2nd line shows number of predicted transmembrane helices (4 in the example above) and total number of residues (123 in the example above). 
* 3rd line shows predicted topology of transmembrane helices in the membrane (currently only antiparallel topology is implemented). 
* 4th line and all lines below show start and end residue numbers of each of the predicted transmembrane helices (current format repeats these numbers twice).

Fragment Picking with Psi Jufo SAM L1 Quota
===========================================

Presenting author: Dominik Gront (dgront at chem dot uw dot edu dot pl)  
Protocol name: fragment picker : CS-Rosetta style  
Brief description: The protocol substitutes nnmake  

Source code location
--------------------

* Check out the mini SVN: https://svn.rosettacommons.org/source/trunk/mini/
* Fragment picker is located in: https://svn.rosettacommons.org/source/trunk/mini/src/core/fragment/picking
* Applications are in: https://svn.rosettacommons.org/source/trunk/mini/src/apps/pilot/dgront/fragmentpicker

Running the protocol capture
----------------------------

1. Set up the path to minirosetta database
2. Set up the path to vall database
3. Run the picker:
   ```
   picker.linuxgccrelease @2jsvX-ss_sim_l1-Q.flags
   ```
RASREC Evolutionary Restraints
==============================

This protocol capture has been written by Tatjana Braun in January 2015 and 
demonstrates on an exemplary structure how the results presented in the PLos 
ONE RosettaCon collection paper "Combining evolutionary information and an 
iterative sampling strategy for accurate protein structure prediction" by 
Tatjana Braun, Julia Koehler Leman and Oliver Lange (2015) have been generated.


Directory structure
-------------------

The following files and folders are provided in this protocol capture

* `scripts/`: Python scripts for generating distance restraints, 
  analyzing the final results and setting up a RASREC run
* `rosetta_flags/`: All necessary flagfiles for setting up the initial and 
  refinement RASREC run
* `inputs/`: All necessary input files, including EVFold prediction, fasta 
  file, fragment files and native structure
* `outputs/`: Exemplary outputs for all steps carried out in this tutorial

Detailed Instructions
---------------------

The following section is written in style of a tutorial and reflects exactly 
what has been done for each target in the manuscript. This tutorial will 
predict the structure of 1wvn with all necessary input files provided . All 
commands listed in this tutorial are executed from the protocol_capture folder.

For being able to directly copy all commands from this tutorial store the path 
of this tutorial in  in the following variable:

    PROTOCOL=~/path/to/tutorial
    cd $PROTOCOL

Contact Prediction and Restraint File Generation
------------------------------------------------

The contact predictions used in the manuscript have been generated with the 
EVFold webserver (available at  
http://evfold.org/evfold-web/newprediction.do) using standard parameters. 
The results can be downloaded in form of a compressed folder, which is 
subdivided into several subdirectories. The all-by-all residue pairing 
scores are stored in {jobname}_{scoringmethod}.txt in the ev_couplings 
folder. In case of the standard scoring method (PLM), the file is named 
{jobname}_PLM.txt.

An exemplary prediction for 1wvn is provided in 
inputs/ev_couplings/1wvn_PLM.txt.

From this score file, the N top-ranked residue pairing scores having a minimum 
distance of 5 residues are extracted and translated into Rosetta specific 
distance restraints. This can be done with the following two steps:

1.  Generate NxN contact map

        $PROTOCOL/scripts/create_evfold_contactmap.py -i inputs/ev_couplings/1wvn_PLM.txt -f inputs/1wvn.fasta -o 1wvn_PLM.cmp

    This command generates a file called 1wvn_PLM.cmp in you protocol 
    folder. It should be identical to the one provided in 
    $PROTOCOL/outputs/1_contactmap

2.  Translate top scoring pairs into Rosetta specific distance restraints

        $PROTOCOL/scripts/extract_top_cm_restraints.py 1wvn_PLM.cmp -r_fasta inputs/1wvn.fasta -r_num_perc 1.0

    This command generates the restraint file called 
    "1wvn_PLM_ub8_lb1.5_w1_m5_CB_70_SIGMOID.cst". An already precomputed 
    restraints file is stored in $PROTOCOL/outputs/2_restraints

Structure Prediction with the RASREC protocol
---------------------------------------------

The Rosetta software package version 3.6 or higher has to be obtained from 
www.rosettacommons.org. Rosetta applications are denoted with the extension 
<.ext>, which should be replaced with the system and compiler dependent 
extension. For instance, for gcc compiled Rosetta on a Linux system use 
.linuxgccrelease.

Note: The RASREC protocol requires MPI with a minimum of 4 computes cores 
(higher numbers are highly recommended). Rosetta can be compiled in MPI mode 
with the following commands:

    cd /path/to/Rosetta/main/source
    ./scons.py -j <number of processors> bin mode=release extras=mpi

RASREC requires substantial computer resources. The required time depends on several factors including size, fold complexity, number and information content of restraints. For instance, target 1r9h requires 26 hours on 96 compute cores (2.6 GHz AMD Opteron Processors).
The run time can be reduced by decreasing the standard pool size of 500, however this is not recommended as this will directly affect the final prediction accuracy.

### Fragment Selection

We have run the fragment picker for all our targets with the following command:

    make_fragments.pl -nohoms

The flag -nohoms leads to exclusion of fragments from homologous proteins. This 
flag should be omitted when not used for benchmarking.

Alternatively, fragments can be generated using the webserver Robetta 
(available at http://www.robetta.org/).

For the tutorial, this step can be omitted as fragments are already provided in 
$PROTOCOL/inputs/fragments.

### Starting a RASREC run

All runs in our manuscript have been set up with the CSRosetta Toolbox 
(available at http://www.csrosetta.org/). In case, the Toolbox is not 
available, a folder containing all necessary flag files is provided. Both ways 
(with and without CSRosetta toolbox) to set up a RASREC run will be described 
below.

* Using CSRosetta Toolbox

  The following path assembles target related input files and stores the setup 
  in ~/cs_targetlib.

        setup_target -method rasrec -fasta inputs/1wvn.fasta -frags inputs/fragments/* -native inputs/1wvnA_ref.pdb -native_restrict  inputs/1wvnA_core.rigid -target 1wvn_standard -restraints 1wvn_PLM_ub8_lb1.5_w1_m5_CB_70_SIGMOID.cst.

  A run ready directory as specified with -dir containing all input files and 
  flags is created with the following command:

        setup_run -method rasrec -target 1wvn_standard -pool_size 500 -dir RASREC_runs -extras mpi

  After using this command, the final run ready directory will be 
  RASREC_runs/1wvn/. Along with all necessary input files, run-scripts for 
  different cluster settings have been generated. 

  In case a local machine is used, RASREC can be started with the following commands:

        cd $PROTOCOL/RASREC_runs/1wvn_standard/run
        mkdir logs
        mpirun -np <CORES> minirosetta.mpi.<ext> -out:level 300 -mute all_high_mpi_rank_filebuf -out:mpi_tracer_to_file logs/log -out:file:silent decoys.out @flags_denovo @flags_rasrec @flags_iterative -run:archive -out:nstruct <CORES-1>

* Without the use of the CSRosetta Toolbox

  In case no CSRosetta toolbox is installed on your system, you can setup the 
  RASREC run manually. All flag files are provided in the rosetta_flags. These 
  flags are ready-to-run for this tutorial. All necessary modifications for 
  different targets will be explained below.

  The following commands will create a run folder containing all necessary flags and files for a RASREC run

        cd $PROTOCOL
        mkdir -p RASREC_runs/1wvn_standard/run/logs
        cp rosetta_flags/standard/* RASREC_runs/1wvn_standard/run
        cd RASREC_runs/1wvn_standard/run

  The final RASREC run can either be started with one of the run-scripts (for different cluster settings) or by executing the following commands:

        mpirun -np <CORES> /path/to/Rosetta/main/source/bin/minirosetta.mpi.<ext> -out:level 300 -mute all_high_mpi_rank_filebuf -out:mpi_tracer_to_file logs/log -out:file:silent decoys.out @flags_denovo @flags_rasrec @flags_iterative -run:archive -out:nstruct <CORES-1>

  For running RASREC with user specific files, the following flags have to be adapted:

  flags_denovo:

        -frag3 <3mer fragment file>
        -frag9 <9mer fragment file>
        -in:file:fasta <input.fasta>

  flags_rasrec:

        -broker:setup <broker setup file>
        -in:file:native <native.pdb>
        #specifies residues used for RMSD calculation/ can be omitted if all resdidues should be used
        -evaluation:rmsd NATIVE _full <resiudes.rigid> 

  setup_init.tpb:

        file <restraints.cst>

  A folder containing template files and a README showing which flags need to 
  be changed is provided in $PROTOCOL/rosetta_flags/template

### Successful Termination and Analysis

A RASREC run is finished, once the fullatom_pool\_ has been generated. The 
following Error message occurs after a successful RASREC run and can be 
ignored:

    ERROR: quick exit from job-distributor due to flag jd2::mpi_nowait_for_remaining_jobs --- this is not an error

The final RASREC models are stored in fullatom_pool/decoys.out

#### Extract top Scoring Ensemble

Using the CSRosetta Toolbox, the top scoring ensemble can be extracted from the 
pool of final models with the following commands:

    extract_decoys decoys.out -score 30 > low_30.out 
    pack_pdbs -silent low_30.out > low_30.pdb

Without the toolbox, the PDB ensemble can be extracted as follows:

    $PROTOCOL/scripts/silent_data.py decoys.out score description | sort -nk 1 | head -n 30 | awk '{print $2}' > pdb_list.txt    
    score_jd2.default.<ext> -in:file:silent decoys.out -in:file:tags $(cat pdb_list.txt) -rescore:skip -out:file:silent low_30.out
    extract_pdbs.default.<ext> -in:file:silent decoys.out -in:file:tags $(cat pdb_list.txt)
    for i in $(echo batch*.pdb); do echo "MODEL $i"; cat "$i"; echo "ENDMDL"; done >> low_30.pdb

#### Analyze Energy and Rmsd of low-energy decoys

The average energy and rmsd of 30 lowest-energy models can be analyzed with the 
following command

    $PROTOCOL/scripts/silent_data.py low_30.out score rms_full | sort -nk 1 | $PROTOCOL/scripts/median.py

which yields

      median       Q1       Q3       hi       lo
    -280.880 -281.639 -278.708 -277.760 -284.937    #score
       3.091    2.810    3.261    3.621    2.795    #rms_full

#### Analyze convergence of ensemble

The following command shows the residues being converged within 2A

    ensemble_analysis.default.<ext> -in:file:silent low_30.out -wRMSD 2 -rigid:out core_2.rigid -rigid:cutoff 2 -calc:rmsd

The output is as follows:

    main: computed RMSD on    #Converged residues listed below
    main: RIGID 4 38 0 0 0    # 4-38 
    main: RIGID 42 64 0 0 0    # 42-64
    main: RIGID 66 98 0 0 0  # 66 -98
    main: RIGID 107 142 0 0 0 # 107-142
    main: RIGID 164 166 0 0 0 # 164-166
    main:
    main: number of atoms from 189 for mean RMSD: 130
    main: fraction of residues converged: 0.69    
    main: mean RMSD to average structure: 1.94
    main: mean pairwise RMSD: 2.88
    main: mean pairwise RMSD * superposed_fraction_of_atoms^-1: 4.19

### Refinement with RASREC

If the convergence of the initial RASREC run is not sufficient enough (< 90%), 
a second RASREC run can be carried out. This run will reuse restraints from 
both predicted contact map and the previous result    

#### Repick Restraints

The following command generates two restraint files given a model ensemble and 
a contact map. 

    cd $PROTOCOL
    $PROTOCOL/scripts/repick_restraints_final.py -c 1wvn_PLM.cmp -s RASREC_runs/1wvn_standard/run/fullatom_pool/low_30.pdb -p 1 -o restraints    

As output, the following files are generated:

* `restraints_converged_distances.cst`: converged  distances in low_30.pdb 
  translated to strict bounded potentials around the average distance

* `restraints_filtered_contactmaps.cst`: additional restraints from the 
  contactmap that do not completely disagree with the previous results. Here a 
  more widely bounded potential is used.

#### Setup RASREC run

The flags and patches used for the refinement RASREC run are identical to the 
ones listed in Section 2.3). The two RASREC runs only differ in the restraints 
used for structural guiding. The restraint files are added to a RASREC run in 
the broker file.

* Using the CSRosetta Toolbox

  Setup prediction

        setup_target -method rasrec -fasta inputs/1wvnA.fasta -frags inputs/fragments/* -native inputs/*ref.pdb -native_restrict inputs/*_core.rigid -target 1wvn_rerun -restraints *_converged_distances.cst *filtered_contactmaps.cst

  Generate run folder

        setup_run -method rasrec -target 1wvn_rerun -pool_size 500 -dir RASREC_runs -extras mpi

  The final run folder for the refinement run can be found 
  RASREC_runs/1wvn_rerun/

  Please COMBINE_RATIO for the filtered restraints set in the broker file 
  (setup_init) as follows:

        CLAIMER ConstraintClaimer
        file restraints_filtered_contactmaps.cst
        COMBINE_RATIO 2                            # ADD THIS LINE!
        FULLATOM
        CENTROID
        SKIP_REDUNDANT 0
        FILTER_WEIGHT  1.00
        FILTER_NAME restraints_filtered_contactmaps.cst
        END_CLAIMER

  This is done so that potentially wrong restriants do not affect the structure 
  prediction in a wrong way.

  RASREC can either be executed with one of the run-scripts or with the 
  following commands:

        cd RASREC_runs/1wvn_rerun/run
        mkdir logs
        mpirun -np <CORES> minirosetta.mpi.<ext> -out:level 300 -mute all_high_mpi_rank_filebuf -out:mpi_tracer_to_file logs/log -out:file:silent decoys.out @flags_denovo @flags_rasrec @flags_iterative -run:archive -out:nstruct <CORES-1>

* Without the CSRosetta Toolbox

  The folder with input files for the refinement run can be set up as follows:

        mkdir -p RASREC_runs/1wvn_rerun/run/logs
        cp rosetta_flags/rerun/* RASREC_runs/1wvn_rerun/run
        cd RASREC_runs/1wvn_rerun/run
        mpirun -np <CORES> minirosetta.mpi.<ext> -out:level 300 -mute all_high_mpi_rank_filebuf -out:mpi_tracer_to_file logs/log -out:file:silent decoys.out @flags_denovo @flags_rasrec @flags_iterative -run:archive -out:nstruct <CORES-1>

  In case the protocol is applied to different targets, the parameters have to 
  be changed as described in 2.3.2)

Rosetta Version
---------------

SVN Revision 57296


Fragment Picking with Quota
===========================

Presenting author: Dominik Gront (dgront at chem dot uw dot edu dot pl)  
Protocol name: fragment picker : CS-Rosetta style  
Brief description: The protocol substitutes nnmake  

Source code location
--------------------

* Check out the mini SVN: https://svn.rosettacommons.org/source/trunk/mini/
* Fragment picker is located in: https://svn.rosettacommons.org/source/trunk/mini/src/core/fragment/picking
* Applications are in: https://svn.rosettacommons.org/source/trunk/mini/src/apps/pilot/dgront/fragmentpicker

Running the protocol capture
----------------------------

1. Set up the path to minirosetta database
2. Set up the path to vall database
3. Run the picker:
   ```
   picker.linuxgccrelease @quota-protocol.flags
   ```
Spanfile from PDB
=================

Author: Julia Koehler Leman (julia dot koehler1982 at gmail dot com)  
Corresponding PI: Jeffrey J. Gray (jgray at jhu dot edu)  
Last Updated: 12/18/2014  
Rosetta Revision #58069

---

This script generates a spanfile from a PDB file. A span file is a topology file
read into Rosetta and is generally generated from a membrane proteins sequence 
using the OCTOPUS server (http://octopus.cbr.su.se/) which is then converted into
a spanfile using octopus2span.pl. This app generates a spanfile from a PDB file
for easy testing of membrane applications, when the structure is known.

Reference:
* Alford RF, Koehler Leman J, Weitzner BD, Duran A, Elazar A, Tilley D, Gray JJ 
  (2015): An integrated framework advancing membrane protein modeling and 
  design, PLosONE (in preparation)

## Executable/Script

    main/source/bin/spanfile_from_pdb.macosclangrelease

## Generating Inputs 

1. Generate a PDB file where the membrane protein structure (our case 1AFO) is 
   transformed into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php).

2. Clean the PDB file by using clean_pdb.pl in the folder 
   Rosetta/tools/protein_tools/scripts/:

        $ clean_pdb.pl 1AFO_tr.pdb ignorechain

3. An example input is provided in the input folder: 1AFO_AB.pdb

## Running the Application
Run the command.sh script provided in this folder:

    $ ./command.sh

## Example Outputs
This generates three spanfiles in the input (! unfortunately) folder. For this 
demo, these files have been moved into the output folder. 

    1AFO_AB.span   spanfile of full PDB
    1AFO_ABA.span  spanfile of chain A
    1AFO_ABB.span  spanfile of chain B in isolation (i.e. residue numbering starts from 1)

Please check the span file for errors and report errors to 
julia dot koehler1982 at gmail dot com!
