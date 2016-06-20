iThis demo shows how to use the fixbb application to invoke the Rosetta packer to optimize the side-chain
conformations of an input structure.  Although the packer can also be used to design new amino acid
sequences, in this case, we hold the sequence fixed, and only vary side-chain conformations.

This demo was created on 20 June 2016 by Vikram K. Mulligan, Ph.D. (vmullig@uw.edu), as part of the Rosetta 2016 Documentation XRW.

BASIC RUN:

After compiling Rosetta, you can run this demo with the following command:

<path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb -in:file:fullatom -resfile resfile.txt -nstruct 5 >log.txt 2>err.txt &

(Note that the executable file's suffix will have to be modified depending on your operating system and compiler.  For example, on Macs with the clang compiler, the suffix would be macosclangrelease instead of linuxgccrelease.)

The following files should be produced:

        1l2y_0001.pdb
        1l2y_0002.pdb
        1l2y_0003.pdb
        1l2y_0004.pdb
        1l2y_0005.pdb
        log.txt
	err.txt
        score.sc

UNDERSTANDING THE RUN:

Open the 1l2y_0001.pdb (output) file and the 1l2y.pdb (input) file in PyMol or another molecular viewer.  Note the changes to side-chain rotamers.  Although the optimal rotamer found by Rosetta for certain residues, such as tryptophan 6, is close to that in the input structure, other, more surface-exposed residues, vary more.  This is to be expected.

Open all of the output PDB files (1l2y_0001.pdb through 1l2y_0005.pdb).  Most likely, these will all be identical.  The packer is a stochastic algorithm, but for a system this small, its output usually converges to the global optimum.

VARYING THE INPUTS:

Now let's understand the various options that we're passing to the fixbb application.  In the command that we ran above, the "-in:file:s 1l2y.pdb" option specifies the input PDB file.  The application loads this file and then calls the packer to run on this input structure, optimizing the side-chains while keeping the backbone fixed.  The "-in:file:fullatom" option tells the application that the input file is a full-atom respresentation (as opposed to a centroid-only model).  This option can be omitted in this case, since it is generally assumed for PDB files.  The "-resfile resfile.txt" option specifies a resfile -- a special type of Rosetta input for controlling the packer.  In this case, our resfile tells the packer to keep the amino acid identity fixed at each position, and to vary only the side-chain conformation.  Finally, "-nstruct 5" tells Rosetta to repeat the protocol five times, producing five independent outputs.  In the case of stochastic algorithms, repeated sampling is important to confirm that an algorithm is successfully converging to the global optimum, or to produce many samples near the global optimum in cases in which the optimum can't quite be reached.

The ">log.txt 2>err.txt" part directs the output to the file log.txt, and any error messages to the file err.txt.  If everything executes properly, err.txt should be empty.

Finally, the terminal ampersand ("&") allows the application to run as a background process, so that you can continue working at the command-line while it executes.

Let's run the application again, but this time, append the flags "-ex1 -ex2 -ex3 -ex4", so that the overall command is:

<path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb -in:file:fullatom -resfile resfile.txt -nstruct 5 >log.txt 2>err.txt &

