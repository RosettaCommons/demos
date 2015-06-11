AbinitioRelax.macosgccrelease -database ~/minirosetta_database -fasta 1not_.fasta -frag3  aa1not_03_05.200_v1_3  -frag9 aa1not_09_05.200_v1_3  -out:file:silent 1not_abrelax_CST_increase_cycles.out -out:file:silent_struct_type binary  -nstruct 1  -cst_file 1not_native_disulf_CEN.cst  -abinitio:relax  -cst_fa_file 1not_native_disulf.cst -native 1not.pdb -increase_cycles 10  -score:weights score12.wts  -ex1 -ex2 -extrachi_cutoff 0 > abrelax.log

relax.macosgccrelease -database ~/minirosetta_database -s idealize_1not.pdb -fasta 1not_.fasta -frag3  aa1not_03_05.200_v1_3  -frag9 aa1not_09_05.200_v1_3  -out:file:silent 1not_native_relax.out -out:file:silent_struct_type binary  -nstruct 1    -abinitio:relax  -cst_fa_file 1not_native_disulf.cst -native 1not.pdb -increase_cycles 10 -score:weights score12.wts  -ex1 -ex2 -extrachi_cutoff 0 > native_relax.log


