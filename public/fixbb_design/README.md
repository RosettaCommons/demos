# Demo for the fixbb application with design

KEYWORDS: DESIGN GENERAL

Demo last modified by Vikram K. Mulligan, Ph.D. (vmullig@uw.edu) during the 2016 Documentation eXtreme Rosetta Workshop (XRW).

This is a demo of a very simple design protocol run on a fixed backbone. If you
have never run Rosetta before, then this is a good first demo to run, because it
is very simple and has few options.

* The demo can be run like this:

        <path_to_Rosetta_directory>/main/source/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb >log.txt &

  The above should take on the order of ten seconds to run.  During that time, you will be able to use the command prompt.  The following files should be produced:

        1l2y_0001.pdb
        log.txt
        score.sc

* Open the structure and the input structure in pymol to observe sequence 
  changes from design.

* Systematically list sequence changes in the form of a sequence profile:

        ls 1l2y_0001.pdb > list.txt  # this would typically be many designed structures all in a list
        python $ROSETTA_TOOLS/protein_tools/scripts/SequenceProfile.py -l list.txt -t 1l2y.pdb

  In the above, the ROSETTA_TOOLS environment variable must be set to point to your Rosetta/tools directory.  Alternatively, you may manually type the location of the Rosetta/tools directory.

* To control which residues are allowed at each sequence position you would add 
  a resfile (included in this demo) like so.  (Note that we're also changing appending
  the suffix "_resout" to the output PDB files so as not to overwrite the files produced
  previously.):

        rosetta/rosetta_source/bin/fixbb.default.linuxgccrelease -in:file:s 1l2y.pdb -resfile resfile.txt -out:suffix _resout > log_resout.txt &

* Open up the resfile.txt file to see its format. Briefly, NATRO leaves the 
  natural rotamer (and amino acid). NATAA leaves the amino acid at a position 
  but allows rotamer to change. ALLAA allows full design with any amino acid. 
  PIKAA followed by a list of single-letter-code amino acids restricts design 
  to just those amino acids.  So, for example:

        1 A PIKAA NT

  indicates that residue 1 can be either N or T.
