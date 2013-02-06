rosetta svn revision 33365

The standard relaxation protocols with membrane flags can be used. The mtssl parameter file needs to be provided. The command line below will relax the double mutant MSBA structure provided in the inputs directory.

bin/relax.linuxgccrelease -database $PATH_TO_ROSETTA_DATABASE -in:file:s input/msba_mtssl_mutant_start_structure.pdb -out:file:fullatom -in:file:extra_res_fa ../input/R1A.params -out:prefix mem_relax_ -run::constant_seed -run::jran 12345 -nstruct 1 -relax:membrane -membrane:normal_cycles 100 -membrane:normal_mag 15 -membrane:center_mag 2 -file:spanfile input/msba.span > & mem_relax.log &
