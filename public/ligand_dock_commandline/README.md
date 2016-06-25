Ligand Dock Demo
================

KEYWORDS: LIGANDS DOCKING

To run:
(where `$ROSETTA3`=path-to-Rosetta/main/source)

    $> $ROSETTA3/bin/ligand_dock.default.linuxgccrelease @flags

Real run has nstruct 500 e.g. and run 10 separate runs for 5000 total runs.

Gives output:
.l
    silent.out

To extract the scores, run this:

    $> $ROSETTA3/src/apps/public/ligand_docking/get_scores.py silent.out > scores.tab

Then you can extract a structure:

    $> $ROSETTA3/bin/extract_atomtree_diffs.default.linuxgccrelease @flagsextract 

Here there's only one structure and one tag for that structure. In a
real case you'd have many more:

    $> $ROSETTA3/src/apps/public/ligand_docking/best_ifaceE.py -n 1 silent.out

For more information, see the ligand dock entry in the manual:  
http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/app_ligand_docking.html

