## Working With Rosetta
#Required Software
	Rosetta relies upon a number of additional programs to function. In order to productively use Rosetta, you must have:
	1. A text editor capable of outputting plain text files; most people use [emacs],[vi],[nano], or some more graphical way of interfacing with one of them.
	2. A molecular viewing tool. Rosetta does not include a way to actually visualize the output files it creates; a tool such as PyMol is necessary to examine the output PDBs.

##The Three Main Interfaces to Rosetta
Rosetta is not a single program; it is instead a large number of individually executable files built to be run either individually or in concert.
#Running the Rosetta Executables directly
	The executables that comprise Rosetta can be run directly by calling them through the command line. These commands are constructed by first indicating the path to the executable you want to run, followed by the flags, options, and input files that dictate how you want the program to execute. 
	As an example, consider the command "$ ./main/source/bin/fixbb.linuxgccrelease -s 1ubq.pdb"
	The first part of this command is the path to and name of the executable we want to run; here, we want to run fixbb, located in main/source/bin. On your platform, the executable may not be named fixbb.linuxgccrelease; to find the version of the executable that will run, navigate to main/source/bin, type "ls fixbb", and then hit Tab to display the versions of fixbb that have been compiled on your system. This works for any other executable as well.
	The second part, -s 1ubq.pdb, is an option, or *flag*, that modifies the execution of the program. In this case, the "-s" option indicates that we want to run fixbb on a single structure; the following "1ubq.pdb" is simply the filename of that structure.
#Running Rosetta via RosettaScripts
	Running RosettaScripts is similar to running an executable directly from the command line: call rosetta_scripts and pass in the relevant options. The key difference is that the executable run is always rosetta_scripts, which requires a flag "parser:protocol" indicating the XML script that it is to run. RosettaScripts is covered in more detail [here].
#Running Rosetta via PyRosetta
	PyRosetta functions are run within Python scripts.
