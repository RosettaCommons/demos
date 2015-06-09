rosetta svn revision 33365

The standard relax protocols can be used. The mtssl parameter file must be provided. The command line below will relax the t4-lysozyme double mutant provided in the input directory.

../bin/relax.linuxgccrelease -database $PATH_TO_ROSETTA_DATABASE -in:file:s input/lysozyme_start_mtssl_mutant.pdb -out:file:fullatom -in:file:extra_res_fa ../input/R1A.params -out:prefix relax_ -run::constant_seed -run::jran 12345 -nstruct 1 >& relax.log 
