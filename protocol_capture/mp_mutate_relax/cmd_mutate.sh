~/Rosetta/main/source/bin/mp_mutate_relax.macosclangrelease \
-database ~/Rosetta/main/database \
-in:file:s input/helix.pdb \
-in:file:native input/helix.pdb \
-mp:setup:spanfiles input/helix.span \
-mp:mutate_relax:mutant_file input/mutations.mut \
-mp:mutate_relax:repack_mutation_only true \
