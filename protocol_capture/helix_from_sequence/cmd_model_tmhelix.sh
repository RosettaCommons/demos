~/Rosetta/main/source/bin/helix_from_sequence.macosclangrelease \
-database ~/Rosetta/main/database \
-in:file:fasta input/tm_helix.fasta \
-mp:setup:transform_into_membrane true \
-mp:transform:optimize_embedding true \
-relax:range:angle_max 0.3 \
-relax:range:nmoves nres \
-out:nstruct 3 \
