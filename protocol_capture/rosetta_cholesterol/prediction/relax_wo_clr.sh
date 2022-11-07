Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
-in:file:s 4ib4_A.pdb \
-mp:setup:spanfiles 4ib4_A.span \
-parser:protocol relax_wo_clr.xml \
-relax:jump_move true \
-packing:pack_missing_sidechains false \
-mp:scoring:hbond true \
-relax:constrain_relax_to_start_coords \