==SHORTENING LOOPS IN ROSETTA==

To shorten loops in Rosetta, you will edit the PDB of the original loop to delete the undesired residues, build a loop and fragment file for the newly short protein, and run the loop through loop modeling with build_initial active.

==Loop shortening the PDB==
Take your input PDB and delete the residue lines related to the residue you are deleting.  Here we are shortening 1FNA.pdb by deleting the ALA with the residue number 83. For the deletion, open the PDB file with your favourite text editor and remove the ATOM entry lines of residue 83. In addition remove also all non ATOM entry lines, i.e. lines not starting with the term ATOM.

==Generating Rosetta inputs==
1FNA_del.pdb is the primary Rosetta input. 

We also need a fragments file and a loop file.  The fragments file is best created via the Robetta server, http://robetta.bakerlab.org/fragmentqueue.jsp, or using the fragment insertion tutorial.  You'll need a FASTA file of your protein to generate fragments; don't forget to delete the deleted-residue from your FASTA file and also note that their is shift in the amino acid number between the fasta file and the PDB file. You can find the FASTA file of the PDB protein in starting_files.

The loop file format is:
LOOP START STOP CUT SAMPLE EXTENDED
In our case, the loop start and stop are 71 (T76) and 81 (87P), respectively.  NOTE that this file is in ROSETTA numbering, not in PDB numbering.  The cutpoint is 77, the position before the removed residue (A83 was removed; 77 is P82).


==Running Rosetta==

path/to/loopmodel.[platform][compiler][mode] @rosetta_inputs/options -database path/to/rosetta_database

-database: Specify the path to the Rosetta database, required for any Rosetta simulation

@rosetta_inputs/options: file holding all rosetta commandline flags. See section "Option file" below. 

The only executeable we'll need here is loop modeling. Briefly, we will run loop modeling with the build_initial mode active.  This option triggers Rosetta to rebuild the loop in a closed state. To do the bulk of the loop remodelling, we will choose KIC remodelling.  CCD would work instead.

==Option file==
An options file has been provided (rosetta_inputs/options), annotated with a description of what each flag is doing.

#necessary for pretty much all loop modeling runs to read in PDBs properly
-in:file:fullatom

#path to input pdb
-loops:input_pdb rosetta_inputs/1FNA_del.pdb

#path to loops file
-loops:loop_file rosetta_inputs/loop_file

#what sizes are the fragments?  9 and 3 are traditional.  The flag seems to require a third argument, but you can pass no fragments in that size.
-loops:frag_sizes 9 3 1

#paths to the fragments in the same vein as previous - none for 1mer fragments
-loops:frag_files rosetta_inputs/aa1FNA_09_05.200_v1_3 rosetta_inputs/aa1FNA_03_05.200_v1_3 none

#this flag triggers build initial mode, which fixes the broken loop before re-solving it
-loops::build_initial

#These flags use KIC loop modeling to remodel the loops
-loops:remodel perturb_kic
-loops:refine refine_kic

#These flags could be used instead for CCD remodeling
#-loops:remodel perturb_ccd
#-loops:refine refine_ccd

#output directory
-out:path sample_output

#prefix for output.  Would not be necessary if someone would rewrite loop modeling to use jd2.
-out:prefix 1FNA_del_

#This option controls how many output structures you get; larger is better!  1 is used here because the tutorial can only take so long; in production you'd use 10000 or more
-nstruct 1

==Interpreting results==
You will get PDB files with energy scores at the bottom of the file as your results. To find the best structures, sort the PDB by their total score (using script sort_by_score), and manually examine the top 5% (or top 100, or whatever you have time for) in PyMOL or another viewer.  Using a combination of total score and your protein intuition (this is an art, not a science), pick which you think is best.

Scientifically, this tutorial is underpowered - at best you could use these results to determine IF a loop can be closed after residues have been deleted.  To do true loop modeling, you would first close the loop using a quick protocol like this, then run a more rigorous loop modeling protocol (covered in the main loop modeling documentation)
