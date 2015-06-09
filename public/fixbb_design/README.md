This is a very simple design on a fixed backbone demo. If you have never run Rosetta before this is a good first demo to run, because it is very simple and has few options.

Use the files from the integration test, for example copying them to a new working directory:
rosetta/rosetta_tests/integration/tests/fixbb

## Run like this
rosetta/rosetta_source/bin/fixbb.linuxgccrelease @flags_fullatom_dun10 -database ~/rosetta/rosetta_database/ > log.txt &

## output will be:
1l2y_0001.pdb
log.txt
score.sc

## open the structure and the input structure in pymol to observe sequence changes from design

## systematically list sequence changes in the form of a sequence profile
ls 1l2y_0001.pdb > list.txt
(this would typically be many designed structures all in a list)
python rosetta/rosetta/tools/protein_tools/scripts/SequenceProfile.py -l list.txt -t 1l2y.pdb

## to control which residues are allowed at each sequence position you would add a resfile (included here) like this:
rosetta/rosetta_source/bin/fixbb.linuxgccrelease @flags_fullatom_dun10 -database ~/rosetta/rosetta_database/ -resfile resfile.txt -out:suffix _resout > log_resout.txt &

## open up the resfile.txt file to see its format. Briefly, NATRO leaves the natural rotamer (and amino acid). NATAA leaves the amino acid at a position but allows rotamer to change. ALLAA allows full design with any amino acid. PIKAA followed by a list of single-letter-code amino acids restricts design to just those amino acids.
## So e.g.
1 A PIKAA NT
indicates that residue 1 can be either N or T