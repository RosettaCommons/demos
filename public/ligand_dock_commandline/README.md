Ligand Dock Demo
================

To run:

    [path]/rosetta/rosetta_source/bin/ligand_dock.linuxgccrelease -database ~/rosetta/rosetta_database/ @flags

Real run has nstruct 500 e.g. and run 10 separate runs for 5000 total runs.

Gives output:

    silent.out

To extract the scores, run this:

    [path]/rosetta/rosetta_source/src/apps/public/ligand_docking/get_scores.py
    <silent.out > scores.tab

Then you can extract a structure:

    ~/rosetta/rosetta_source/bin/extract_atomtree_diffs.macosgccrelease
    @flagsextract -database ~/rosetta/rosetta_database

Here there's only one structure and one tag for that structure. In a
real case you'd have many more:

    [path]/rosetta/rosetta_source/src/apps/public/ligand_docking/best_ifaceE.py
    -n 1 silent.out

For more information, see the ligand dock entry in the manual:
http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/app_ligand_docking.html

