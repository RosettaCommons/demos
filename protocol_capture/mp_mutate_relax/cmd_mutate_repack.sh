~/Rosetta/main/source/bin/mp_mutate_relax.macosclangrelease \
-database ~/Rosetta/main/database \
-in:file:s input/helix.pdb \
-in:file:native input/helix.pdb \
-mp:setup:spanfiles input/helix.span \
-mp:mutate_relax:mutation F15C \
-mp:mutate_relax::repack_radius 8.0 \
