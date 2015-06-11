rosetta svn revision 33365

The standard fixed backbone design protocol can be used. The mtssl parameter
file must be provided. The residue packing file in the config directory
specifies that only the non-cannonical residue type R1A is allowed at
positions 59 and 159. The command below will generate the 59/159 double mtssl
mutant of the t4-lysozyme protein given in the input directory.

bin/fixbb.linuxgccdebug -database $PATH_TO_ROSETTA_DATABASE -in:file:s input/lysozyme_pseudo_wildtype.pdb -out:file:fullatom -in:file:extra_res_fa ../input/R1A.params  -resfile config/resfile.pack -out:prefix mtssl_mutant_ >& make_mutant.log &
