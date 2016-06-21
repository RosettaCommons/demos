# to minimize entire structure:
stepwise_protein_test.linuxgccrelease  -s1 2ci2.pdb  -global_optimize -fixed_res `seq 1 62` -input_res1 `seq 1 62` -fixed_res `seq 1 62` -database ~/minirosetta_database -out:file:silent 2ci2_min.out -fasta 2ci2.fasta

# to carry out KIC loop modeling
loopmodel.macosgccrelease -database ~/minirosetta_database -loops:remodel perturb_kic -loops:refine refine_kic -loops:input_pdb 2ci2_min.pdb -in:file:native 2ci2.pdb -loops:loop_file 2ci2_35_45.loop -loops:max_kic_build_attempts 10000 -in:file:fullatom -out:file:fullatom -out:prefix 2ci2 -out:pdb -ex1 -ex2 -extrachi_cutoff 0 -out:nstruct 1 -out:file:silent_struct_type binary  -out:file:silent 2ci2_kic_loop35_45_from_2ci2_min.out  > kic.log
