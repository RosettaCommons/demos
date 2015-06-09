#!/bin/bash

cd generation2

for file in *.pdb
do

stem=`echo $file | sed 's/\.pdb//'`

$ROSETTA3_SRC/bin/loopmodel.default.linuxgccrelease \
	-database $ROSETTA3_DB \
	-in::file::s $stem.pdb \
	-loops::loop_file $stem.loopfile \
	-nstruct 1 \
	-in::file::fullatom \
	-loops::remodel quick_ccd \
	-loops::intermedrelax no \
	-loops::refine no \
	-loops::relax relax \
	-loops::frag_sizes 9 3 1 \
	-loops::frag_files ../aa1q0p_09_05.200_v1_3.gz  ../aa1q0p_03_05.200_v1_3.gz none \
	-loops::extended \
	-loops::random_grow_loops_by 4 \
	-relax::default_repeats 1 \
	-relax::jump_move true \
	-edensity::mapfile ../1q0p.5A.mrc \
	-edensity::realign min \
	-edensity::mapreso 5.0 \
	-edensity::grid_spacing 2.5 \
	-edensity::whole_structure_allatom_wt 0.05 \
	-overwrite

done

cd ..

