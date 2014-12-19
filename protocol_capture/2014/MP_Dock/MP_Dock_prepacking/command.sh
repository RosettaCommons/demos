/Users/julialeman/Documents/julia/git_final/Rosetta/main/source/bin/docking_prepack_protocol.macosclangrelease \
-database /Users/julialeman/Documents/julia/git_final/Rosetta/main/database \
-in:file:s input/1AFO_AB.pdb \
-membrane_new:setup:spanfiles input/1AFO_AB.span \
-out:file:scorefile output/score_ppk_1AFO.sc \
-out:path:pdb output \
-nstruct 1 \
-score:weights mpframework_smooth_fa_2014.wts \
-packing:pack_missing_sidechains 0 \
>& output/log_ppk_1AFO.log &
