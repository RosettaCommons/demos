Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
-in:file:s 4ib4_A_dock.pdb \
-mp:setup:spanfiles 4ib4_A.span \
-in:file:extra_res_fa CLR.params \
-packing:pack_missing_sidechains false \
-parser:protocol dock.xml \
-nstruct 10 \