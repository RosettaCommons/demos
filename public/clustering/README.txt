Cluster demo: Run the following command

The following command is a good place to start:
1) /path/to/rosetta/bin/cluster.linuxgccrelease -database /path/to/database -in:file:s rosetta_inputs/*.pdb -in:file:fullatom -cluster:gdtmm -cluster:radius -1 -cluster:population_weight 0.0 -sort_groups_by_energy 

extra commands would be if the input was in silent file:
-in:file:silent <filename.out> -in:file:silent_struct_type <desired type but binary is best choice>
-limit_cluster_size <n>
-limit_clusters <n>
-cluster:radius <for gdt you want some number between 10-50>
-out:file:silent clustered.out -out:file:silent_struct_type binary
(You can extract the silent file by running the mini app
extract_pdbs.linuxgccrelease)

2)If you have thousands of structures it may be best to first do an energy
cut to filter out high energy structures. The pdb's should be first scored and stored in an out file. Try using the
following app
score_jd2 -in:file:s *.pdb -out:file:silent scored.out
-out:file:silent_struct_type binary -in:file:fullatom

After the pdbs have been scored use the included script to do an energy cut.
make_sub_silent_file_percentile.py scored.out ecut_10.out -1 0.10
where 0.10 is a 10% energy cut.

 
