# Working With Rosetta

KEYWORDS: CORE_CONCEPTS GENERAL

Written by Frank Teets
Last Modified Jun 21 2016

##Required Software

To productively use Rosetta, a number of additional programs are need:

1. A text editor capable of outputting plain text files; most people use [Emacs](https://www.gnu.org/software/emacs/manual/html_node/emacs/index.html), [Vim](https://www.washington.edu/computing/unix/vi.html), [Nano](https://www.nano-editor.org/dist/v2.0/nano.html), or similar. Note that *word processors* like Microsoft Word and Libre/Open Office are unsatisfactory, as the additional formating they introduce will confuse Rosetta. 
2. A molecular viewing tool. Rosetta does not include a way to actually visualize the output files it creates; a tool such as [PyMOL](https://www.pymol.org/) or [Chimera](https://www.cgl.ucsf.edu/chimera/) is necessary to examine the output PDBs. Rosetta does include a PyMOLObserver for directly viewing its output in PyMOL; instructions for setting up PyMOLObserver are found [here](https://www.rosettacommons.org/docs/latest/rosetta_basics/graphics-and-guis).
3. A terminal. Unix and Apple users have suitable terminals installed by default. Windows users attempting to run Rosetta remotely on a Mac or Unix machine will require a tool like PuTTY to provide the necessary interface.

It may also help to familiarize yourself with the basics of [command line interfaces](https://bash.cyberciti.biz/guide/Main_Page). Particularly useful commands include:

	ls foo			displays a list of files in the directory foo
	cd foo			moves the current working directory to the directory foo
	mv foo bar 		moves all of foo into bar if bar is a directory or into a file called bar if not
	cp foo bar		moves all of foo into bar if bar is a directory or into a file called bar if not
	ln -s <path_to_bar> foo establishes a link called a symbolic link or symlink from foo to bar
	top 			lists the currently running processes, together with their process IDs. This is useful for how much computational resource your process (and others are consuming). They can also be used to identify your process ID
	kill foo_pid	        stops the process with process ID foo_pid
	nohup			when prepended to a command, runs that command in the background, preventing it from writing output to the terminal
	

##The Three Main Interfaces to Rosetta

It is somewhat incorrect to think of Rosetta as a single program. Instead, Rosetta is more a library of functionality to model biomacromolecules. This functionality can be accessed either through a number of executable files built to run either individually, sequentially, or in concert. Alternatively, the functionality may be customized through an XML interface called RosettaScripts or through Python via PyRosetta.

###Running the Rosetta Executables Directly

The executables that comprise Rosetta can be run directly by calling them through the terminal. These commands are constructed by first indicating the path to the executable you want to run, followed by the flags, options, and input files that dictate how you want the program to execute. 

As an example, consider the command 

	$> $ROSETTA3/bin/fixbb.default.linuxgccrelease -s 1ubq.pdb -out:suffix _my_first_rosetta_output @general_fixbb_flags

The first part of this command is the path to and name of the executable we want to run; here, we want to run the executable `fixbb`, located in `$ROSETTA3/bin`, which will attempt to design residues onto a peptide backbone provided to it. On your platform, the executable may not be named `fixbb.linuxgccrelease`; to find the version of the executable that will run, navigate to `$ROSETTA3/bin` and run `ls fixbb*` to display the versions of `fixbb` that have been compiled on your system. This works for any other executable as well.

The second part, `-s 1ubq.pdb`, is an option, or _flag_, that modifies the execution of the program. In this case, the `-s` option indicates that we want to run `fixbb` on one or more structures listed on the command line; the `1ubq.pdb` is the argument to the `-s` option, and indicates the filename of the structure to use. An alternative to the `-s` option is `-l`, which takes as its argument the name of a file containing a newline-delimited list of input files and is broadly conceptually equivalent to running multiple `-s` commands in sequence. Other common options include:

	-nstruct			Takes as an argument the number of output structures Rosetta is to generate per input structure; this is important for Monte Carlo simulations.
	-database			The location of the [Rosetta database](https://www.rosettacommons.org/docs/latest/rosetta_basics/database). On most installations this can be automatically detected, so this option should not be required.
	-ignore_unrecognized_residue	A Boolean (true/false) option which doesn't require an argument, this option tells Rosetta to ignore any residue (like small molecule ligands and crystallization adducts) which it doesn't recognize.

The third part of this command is an options file, or *flags file*, which is simply a newline-delimited list of options. It is functionally equivalent to writing each option individually.

#### Example Run

In order to demonstrate how Rosetta executables may be run and how their execution controlled via options, create a folder named `tutorials` or something similar, and put the files `1ubq.pdb` and `1qys.pdb` there. 

Now run 

	$> $ROSETTA3/bin/fixbb.default.linuxgccrelease -s 1ubq.pdb -out:suffix _my_second_rosetta_output @general_fixbb_flags

from within the tutorials directory you just created and observe that the program generates an output pdb file(`1ubq_0001.pdb`), a score file (`score.sc`) and a log file (`output`).

Try running the following commands:

	$> $ROSETTA3/bin/fixbb.default.linuxgccrelease -s 1ubq.pdb -nstruct 10 -out:suffix _my_third_rosetta_output @general_fixbb_flags
	
	$> $ROSETTA3/bin/fixbb.default.linuxgccrelease -s 1ubq.pdb -out:suffix _my_fourth_rosetta_output -ignore_unrecognized_residue @general_fixbb_flags

The first command should produce ten output files; the second should run similarly, albeit not identically, to the command provided at the start of this tutorial.

Also, try making a file named `pdblist` containing a list of PDBs including 1ubq and 1qys in the following format:
	
	1ubq.pdb
	1qys.pdb
	
and running

	$> $ROSETTA3/bin/fixbb.default.linuxgccrelease -l pdblist

##Running Rosetta via RosettaScripts

Running RosettaScripts is similar to running an executable directly from the command line: execute `rosetta_scripts` and pass in the relevant options. The key difference is that the executable run is always rosetta_scripts, which requires the option `-parser:protocol` indicating the XML script that it is to run. RosettaScripts is covered in more detail [here](). Try creating a file named `fixbb_script.xml` containing the following:

	<ROSETTASCRIPTS>
	script goes here
	</ROSETTASCRIPTS>

saving it, and running it with 

	$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -parser:protocol fixbb_script.xml. You should see output similar to the command run at the start of this tutorial.

##Running Rosetta via PyRosetta

PyRosetta functions are run within Python scripts. Detailed PyRosetta tutorials can be found [here](http://www.pyrosetta.org/tutorials).

##Where to Find Things In Rosetta

####Binaries
The Rosetta executables are all symlinked into Rosetta/source/bin and may be run from there; the actual executables are compiled into Rosetta/source/build, but should be run from source/bin to ensure robustness to updates.

####Source Code
The Rosetta source code exists in Rosetta/source/src, but does not contain content required for non-developer users of Rosetta.

####Weights files
Rosetta weights files (.wts files), used to parameterize the scorefunction weights, are in Rosetta/main/database/scoring/weights .

####Params files
Parameter files are in the corresponding subdirectories (fa_standard or centroid) of /Rosetta/main/database/chemical/residue_type_sets .

####Tools
The tools directory, Rosetta/tools, contains a number of useful tools for manipulating Rosetta inputs and outputs.

####Scripts
Example RosettaScripts XMLs and useful python scripts may be found in Rosetta/main/source/scripts

##Troubleshooting
Rosetta will occasionally throw an informative error message rather than successfully completing a job; these may be indicative of errors in the input or options provided and are not necessarily bugs. Segmentation faults certainly indicate a bug in Rosetta itself, and should be reported to the [Rosetta forums](https://www.rosettacommons.org/forum)
