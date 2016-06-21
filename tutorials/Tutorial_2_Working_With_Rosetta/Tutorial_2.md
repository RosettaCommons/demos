# Working With Rosetta
##Required Software
Rosetta relies upon a number of additional programs to function. In order to productively use Rosetta, you must have:

1. A text editor capable of outputting plain text files; most people use [emacs](https://www.gnu.org/software/emacs/manual/html_node/emacs/index.html), [vi](https://www.washington.edu/computing/unix/vi.html), [nano](https://www.nano-editor.org/dist/v2.0/nano.html), or some more graphical way of interfacing with one of them.
2. A molecular viewing tool. Rosetta does not include a way to actually visualize the output files it creates; a tool such as [PyMOL](https://www.pymol.org/) is necessary to examine the output PDBs. Rosetta does include a PyMOLObserver for directly viewing its output in Pymol; instructions for setting up PyMOLObserver are found [here].
3. A terminal. Unix and Apple users have suitable terminals installed by default. Windows users attempting to run Rosetta remotely on a Mac or Unix machine will require a tool like PuTTY to provide the necessary interface.

It may also help to familiarize yourself with the basics of [command line interfaces](https://bash.cyberciti.biz/guide/Main_Page). Particularly useful commands include:

	ls foo			displays a list of files in the directory foo
	cd foo			moves the current working directory to the directory foo
	mv foo bar 		moves all of foo into bar if bar is a directory or into a file called bar if not
	cp foo bar		moves all of foo into bar if bar is a directory or into a file called bar if not
	ln -s foo 		<path_to_bar> establishes a link called a symbolic link or symlink from foo to bar
	top 			lists the currently running processes, together with their process IDs. This is useful for how much computational resource your process (and others are consuming). They can also be used to identify your process ID
	kill foo_pid	stops the process with process ID foo_pid
	nohup			when prepended to a command, runs that command in the background, preventing it from writing output to the terminal
	

#The Three Main Interfaces to Rosetta
Rosetta is not a single program; it is instead built with a large number of executable files built to be run either individually, sequentially, or in concert. These files may be run directly from the command line, through an XML interface called RosettaScripts, or through Python via Pyrosetta.
##Running the Rosetta Executables directly
The executables that comprise Rosetta can be run directly by calling them through the terminal directly. These commands are constructed by first indicating the path to the executable you want to run, followed by the flags, options, and input files that dictate how you want the program to execute. 
As an example, consider the command 

	$> <path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -s 1ubq.pdb -out:suffix _my_first_rosetta_output @general_fixbb_flags

The first part of this command is the path to and name of the executable we want to run; here, we want to run the executbale `fixbb`, located in `<path_to_Rosetta_directory>/main/source/bin`, which will attempt to design residues onto a peptide backbone provided to it. On your platform, the executable may not be named `fixbb.linuxgccrelease`; to find the version of the executable that will run, navigate to `<path_to_Rosetta_directory>/main/source/bin`, type `ls fixbb`, and then hit _Tab_ to display the versions of `fixbb` that have been compiled on your system. This works for any other executable as well.

The second part, `-s 1ubq.pdb`, is an option, or _flag_, that modifies the execution of the program. In this case, the `-s` option indicates that we want to run `fixbb` on a single structure; the following `1ubq.pdb` is simply the filename of that structure. An alternative to the `-s` option is `-l`, which takes as its argument a newline-delimited list of input files and is broadly conceptually equivalent to running multiple `-s` commands in sequence. Other common options include:

	-nstruct	This option dictates the number of output structures Rosetta is to generate; this is important for [Monte Carlo simulations]().
	-database	This option should not be required for a successful Rosetta install, but indicates the path to the Rosetta database (which is now detected automatically.)
	-out::file	This option indicates the name and location to which Rosetta is to write any output files.
	-ignore_unrecognized_residue This option should not be required for a successful Rosetta install, but indicates that residues that Rosetta cannot process should be left out of the execution.

In order to demonstrate how Rosetta exectuables may be run and their execution controlled via flags, create a folder in your main Rosetta directory named `tutorials` or something similar, and there put the files `1ubq.pdb` and `1qys.pdb`. 

Now run 

	$> <path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -s 1ubq.pdb -out:suffix _my_second_rosetta_output @general_fixbb_flags

from within the tutorials directory you just created and observe that the program generates an output pdb file(`1ubq_0001.pdb`), a score file (`score.sc`) and a log file (`output`).
Try running the following commands:

	$> <path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -s 1ubq.pdb -nstruct 10 -out:suffix _my_third_rosetta_output @general_fixbb_flags
	
	$> <path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -s 1ubq.pdb -out:suffix _my_fourth_rosetta_output -ignore_unrecognized_residue @general_fixbb_flags

The first command should produce ten output files; the second should run similarly, albeit not identically, to the command provided at the start of this tutorial.
Also, try making a list of PDBs named `pdblist` including 1ubq and 1qys in the following format:
	
	1ubq.pdb
	1qys.pdb
	
and running

	$> <path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -l pdblist

##Running Rosetta via RosettaScripts
Running RosettaScripts is similar to running an executable directly from the command line: execute `rosetta_scripts` and pass in the relevant options. The key difference is that the executable run is always rosetta_scripts, which requires the option `parser:protocol` indicating the XML script that it is to run. RosettaScripts is covered in more detail [here](). Try creating a file named `fixbb_script.xml` containing the following:
	<ROSETTASCRIPTS>
	script goes here
	</ROSETTASCRIPTS>
saving it, and running it with 

	$> <path_to_Rosetta_directory>/main/source/bin/rosetta_scripts.default.linuxgccrelease -parser:protocol fixbb_script.xml. You should see output similar to the command run at the start of this tutorial.

##Running Rosetta via PyRosetta
PyRosetta functions are run within Python scripts; it requires that rosetta be imported and initiliazed with the commands
	from rosetta import /*
	rosetta.init()
or similar, after which the Pyrosetta functions may be called within the script itself. Detailed PyRosetta tutorials can be found [here](http://www.pyrosetta.org/tutorials).

#Constructing a Rosetta Job

##Considerations of Scale
Rosetta is capable of running processes at many scales. The deterministic scoring functions may be run on a single processor in a matter of seconds; running an ab initio folding job can require a national-scale supercomputer running for thousands of CPU-hours. It may be worthwhile to review the CPU time required to run Rosetta for different tasks [here]().

##Post-Processing
Because Rosetta is stochastic software, it is often necessary to perform statistical analysis on the structures it generates, called _decoys_. This may be done with a statistics package such as R, documentation for which is [here](https://cran.r-project.org/doc/manuals/r-release/R-intro.pdf).

# Where to Find Things In Rosetta
###Binaries
The Rosetta executables are all symlinked into `<path_to_Rosetta_directory>/main/source/bin` and may be run from there; the actual executables are compiled into `<path_to_Rosetta_directory>/main/source/build`, but should be run from `source/bin` to ensure robustness to updates.
###Source Code
The Rosetta source code exists in `<path_to_Rosetta_directory>/main/source/src`, but does not contain content required for non-developer users of Rosetta.
###Weights files
Rosetta weights files (.wts files), used to parameterize the scorefunction weights, are in `<path_to_Rosetta_directory>//main/database/scoring/weights` .
##Params files
Parameter files are in the corresponding subdirectories of `<path_to_Rosetta_directory>//main/database/chemical/residue_type_sets` .
###Tools
The tools directory, `path_to_Rosetta_directory>/tools`, contains a number of useful tools for manipulating Rosetta inputs and outputs.
###Scripts
RosettaScripts and useful python scripts may be found in `<path_to_Rosetta_directory>/main/source/scripts`.
