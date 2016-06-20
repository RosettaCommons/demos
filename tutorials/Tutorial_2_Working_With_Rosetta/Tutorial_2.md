## Working With Rosetta
#Required Software
Rosetta relies upon a number of additional programs to function. In order to productively use Rosetta, you must have:
1. 	A text editor capable of outputting plain text files; most people use [emacs](https://www.gnu.org/software/emacs/manual/html_node/emacs/index.html),[vi](https://www.washington.edu/computing/unix/vi.html),[nano](https://www.nano-editor.org/dist/v2.0/nano.html), or some more graphical way of interfacing with one of them.
2. 	A molecular viewing tool. Rosetta does not include a way to actually visualize the output files it creates; a tool such as [PyMOL](https://www.pymol.org/) is necessary to examine the output PDBs. Rosetta does include a PyMOLObserver for directly viewing its output in Pymol; instructions for setting up PyMOLObserver are found [here].

##The Three Main Interfaces to Rosetta
Rosetta is not a single program; it is instead a large number of individually executable files built to be run either individually, sequentially, or in concert.
#Running the Rosetta Executables directly
The executables that comprise Rosetta can be run directly by calling them through a [command line interface](https://bash.cyberciti.biz/guide/Main_Page). These commands are constructed by first indicating the path to the executable you want to run, followed by the flags, options, and input files that dictate how you want the program to execute. 
As an example, consider the command '$ ./main/source/bin/fixbb.linuxgccrelease -s 1ubq.pdb'
The first part of this command is the path to and name of the executable we want to run; here, we want to run fixbb, located in main/source/bin. On your platform, the executable may not be named fixbb.linuxgccrelease; to find the version of the executable that will run, navigate to main/source/bin, type "ls fixbb", and then hit Tab to display the versions of fixbb that have been compiled on your system. This works for any other executable as well.
The second part, -s 1ubq.pdb, is an option, or *flag*, that modifies the execution of the program. In this case, the "-s" option indicates that we want to run fixbb on a single structure; the following "1ubq.pdb" is simply the filename of that structure. The counterpart to the -s option is -l, which takes as its argument a newline-delimited list of input files and is broadly conceptually equivalent to running multiple -s commands in sequence. Other common options include:
	-nstruct: This option dictates the number of output structures Rosetta is to generate; this is important for [Monte Carlo simulations].
	-database: This option should not be required for a successful Rosetta install, but indicates the path to the Rosetta database (which is now detected automatically.)
	-out::file [...] This option indicates the name and location to which Rosetta is to write any output files.
	-ignore_unrecognized_residue This option should not be required for a successful Rosetta install, but indicates that residues that Rosetta cannot process should be left out of the execution.
Try running the following commands:
	'$ ./main/source/bin/fixbb.linuxgccrelease -s 1ubq.pdb -nstruct 10'
	'$ ./main/source/bin/fixbb.linuxgccrelease -s 1ubq.pdb -ignore_unrecognized_residue'
The first command should produce ten output files; the second should run identically to the command provided at the start of this tutorial.
Also, try making a list of PDBs including 1ubq and 1a2b and running
	'$ ./main/source/bin/fixbb.linuxgccrelease -l pdblist'
#Running Rosetta via RosettaScripts
Running RosettaScripts is similar to running an executable directly from the command line: execute rosetta_scripts and pass in the relevant options. The key difference is that the executable run is always rosetta_scripts, which requires the option "parser:protocol" indicating the XML script that it is to run. RosettaScripts is covered in more detail [here]. Try creating a file named "fixbb_script.xml" containing the following:
	<ROSETTASCRIPTS>
	script goes here
	</ROSETTASCRIPTS>
saving it, and running it with '$./main/source/bin/rosetta_scripts.linuxgccrelease -parser:protocol fixbb_script.xml'. You should see output similar to the command run at the start of this tutorial.
#Running Rosetta via PyRosetta
PyRosetta functions are run within Python scripts; it requires that rosetta be imported with the command
	from rosetta import /*
or similar, after which the Pyrosetta functions may be called within the script itself. PyRosetta tutorials may be called [here].

##Constructing a Rosetta Job

#Considerations of Scale
Rosetta is capable of running processes at many scales. The deterministic scoring functions may be run on a single processor in a matter of seconds; running an ab initio folding job can require a national-scale supercomputer running for thousands of CPU-hours. It may be worthwhile to review the CPU time required to run Rosetta for different tasks [here].

#Post-Processing
Because Rosetta is stochastic software, it is often necessary to perform statistical analysis on the structures it generates, called *decoys*. This may be done with a statistics package such as R, documentation for which is [here].
