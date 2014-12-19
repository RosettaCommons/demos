/Users/julialeman/Documents/julia/git_final/Rosetta/main/source/bin/mpdocking.macosclangrelease \
-database /Users/julialeman/Documents/julia/git_final/Rosetta/main/database \
-in:file:s input/1AFO_AB_noMEM.pdb \
-in:file:native input/native.pdb \
-membrane_new:setup:spanfiles input/1AFO_AB.span \
-out:file:scorefile output/score_mpdock_1AFO.sc \
-out:path:pdb output \
-nstruct 1 \
-docking:partners A_B \
-dock_pert 3 8 \
-show_simulation_in_pymol 0 \
>& output/log_mpdock_1AFO.log &
