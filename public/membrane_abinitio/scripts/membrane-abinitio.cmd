$ROSETTA3/bin/abinitio2.macosgccrelease \
-abinitio:membrane \
-in:file:native ../rosetta_inputs/BRD4.pdb \
-in:file:frag3 ../rosetta_inputs/aaBRD4_03_05.200_v1_3 \
-in:file:frag9 ../rosetta_inputs/aaBRD4_09_05.200_v1_3 \
-in:file:spanfile ../rosetta_inputs/BRD4.span \
-in:file:lipofile ../rosetta_inputs/BRD4.lips4 \
-membrane:no_interpolate_Mpair \
-membrane:Menv_penalties \
-score:find_neighbors_3dgrid \
-increase_cycles 0.05 \
-out:nstruct 1 \
-out:file:silent BRD4_silent.out \
-run:no_prof_info_in_silentout \
-mute core.io.database \
-mute core.scoring.MembranePotential \
