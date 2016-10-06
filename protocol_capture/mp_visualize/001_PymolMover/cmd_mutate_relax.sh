~/Rosetta/main/source/bin/mp_mutate_relax.macosclangrelease \
-database ~/Rosetta/main/database \
-in:file:s input/1AFO_opm_AB.pdb \
-in:file:native input/1AFO_opm_AB.pdb \
-mp:setup:spanfiles input/1AFO_opm_AB.span \
-mp:mutate_relax:mutation V15C \
-mp:mutate_relax:relax true \
-relax:range:angle_max 0.3 \
-relax:range:nmoves nres \
-mp:transform:optimize_embedding true \
-mp:mutate_relax:nmodels 1 \
-show_simulation_in_pymol 0 \
-keep_pymol_simulation_history true \