# Demo for the fixbb application

KEYWORDS: STRUCTURE_PREDICTION GENERAL

This demo shows how to use the fixbb application to invoke the Rosetta packer to optimize the side-chain
conformations of an input structure.  Although the packer can also be used to design new amino acid
sequences, in this case, we hold the sequence fixed, and only vary side-chain conformations.

This demo was created on 20 June 2016 by Vikram K. Mulligan, Ph.D. (vmullig@uw.edu), as part of the Rosetta 2016 Documentation XRW.

##Basic run:

After compiling Rosetta, you can run this demo with the following command: (where `$ROSETTA3`=path-to-Rosetta/main/source)

```
$> $ROSETTA3/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb -in:file:fullatom -resfile resfile.txt -nstruct 1 >log.txt 
```

(Note that the executable file's suffix will have to be modified depending on your operating system and compiler.  For example, on Macs with the clang compiler, the suffix would be macosclangrelease instead of linuxgccrelease.)

The following files should be produced:

```
1l2y_0001.pdb
1l2y_0002.pdb
1l2y_0003.pdb
1l2y_0004.pdb
1l2y_0005.pdb
log.txt
score.sc
```

At any time, the output files can be deleted selectively by running the clean.sh script ("bash clean.sh").

##Understanding the run:

Open the 1l2y_0001.pdb (output) file and the 1l2y.pdb (input) file in PyMol or another molecular viewer.  Note the changes to side-chain rotamers.  Although the optimal rotamer found by Rosetta for certain residues, such as tryptophan 6, is close to that in the input structure, other, more surface-exposed residues, vary more.  This is to be expected.

Open all of the output PDB files (1l2y_0001.pdb through 1l2y_0005.pdb).  Most likely, these will all be identical.  The packer is a stochastic algorithm, but for a system this small, its output usually converges to the global optimum. Note that for the purpose of this demo, we are using `-nstruct 1` but for real applications you may want to increase this number.

##Varying the inputs:

Now let's understand the various options that we're passing to the fixbb application.  In the command that we ran above, the "-in:file:s 1l2y.pdb" option specifies the input PDB file.  The application loads this file and then calls the packer to run on this input structure, optimizing the side-chains while keeping the backbone fixed.  The "-in:file:fullatom" option tells the application that the input file is a full-atom respresentation (as opposed to a centroid-only model).  This option can be omitted in this case, since it is generally assumed for PDB files.  The "-resfile resfile.txt" option specifies a resfile -- a special type of Rosetta input for controlling the packer.  In this case, our resfile tells the packer to keep the amino acid identity fixed at each position, and to vary only the side-chain conformation.  Finally, "-nstruct 5" tells Rosetta to repeat the protocol five times, producing five independent outputs.  In the case of stochastic algorithms, repeated sampling is important to confirm that an algorithm is successfully converging to the global optimum, or to produce many samples near the global optimum in cases in which the optimum can't quite be reached.

The ">log.txt" part directs the output to the file log.txt.

Finally, the terminal ampersand ("&") allows the application toi run as a background process, so that you can continue working at the command-line while it executes.  Omitting this would result in realtime output to the screen.

Let's delete the output and run the application again, but this time, append the flags "-ex1 -ex2 -ex3 -ex4 -extrachi_cutoff 0 -out:prefix run2_ -out:file:scorefile run2.sc", so that the overall command is:

```
$> $ROSETTA3/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb -in:file:fullatom -resfile resfile.txt -nstruct 1 -ex1 -ex2 -ex3 -ex4 -extrachi_cutoff 0 -out:prefix run2_ -out:file:scorefile run2.sc >log2.txt 
```

You'll notice that this time, the app takes much longer to run.  Examining the output, you'll find a few things.  First, the rotamer chosen for the central tryptophan and for many other side-chains is much closer to the original.  Second, you may see some variation in the five output structures.  So what did we change?  Well, the "-ex1", "-ex2", "-ex3", and "-ex4" flags activate extra rotamers for the first, second, third, and fourth side-chain chi angles, respectively.  The "-extrachi_cutoff" flag determines the number of residue neighbours (nearby residues) that a resiude must have in order to receive the additional rotamers.  By including additional rotamers, it becomes more likely that, at any given position, one of the rotamers included will be very close to the optimal side-chain conformation.  The cost, though, is an increase in the computational complexity (reducing the probability of convergence somewhat), and a concamitant increase in the running time.  Similarly, performing the same task on a larger input structure will also result in a longer running time and a lower probability of convergence.  The "-out:prefix run2_" flag allows us to write PDB files with the prefix "run2_" in the filename, so that we don't try to overwrite the previous output.  Similarly , "-out:file:scorefile run2_.sc" directs the final scores to a new file, so that we don't overwrite the scorefile from the first run.

##Varying the resfile to alter packer behaviour

The same behaviour shown above can also be accomplished by altering the resfile.  If compare resfile2.txt to resfile.txt, you'll notice that the first line is different.  The additional text ("EX 1 EX 2 EX 3 EX 4 EX_CUTOFF 0") accomplishes much the same thing as the extra options that we passed at the commandline earlier.  Deleting the output and running the following shoulkd produce the same behaviour:

```
$> $ROSETTA3/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb -in:file:fullatom -resfile resfile2.txt -nstruct 1 -out:prefix run3_ -out:file:scorefile run3.sc >log3.txt 
```
