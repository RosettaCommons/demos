#Spin Labels

KEYWORDS: EXPERIMENTAL_DATA STRUCTURE_PREDICTION

Edited in June 2016 by Parisa Hosseinzadeh (parisah@uw.edu) to enable automated demo testing.

A copy of the page is in MeilerLab-ProtocolCapture.pdf

The input, config, and bin directories have some common materials for the
substeps of the protocol. Each subdirectory has its own input, config, and
bin directories with specific materials for that step. Each subdirectory also
has its own README.txt with specific information about that step. Each step is
entirely self contained, but you might enjoy working through them in the order
given below. 


Step1: create_mtssl_mutant

The standard fixed backbone design protocol can be used. The mtssl parameter
file must be provided. The residue packing file in the config directory
specifies that only the non-cannonical residue type R1A is allowed at
positions 59 and 159. The command below will generate the 59/159 double mtssl
mutant of the t4-lysozyme protein given in the input directory.
```
$> $ROSETTA3/bin/fixbb.default.linuxgccrelease -in:file:s create_mtssl_mutant/input/lysozyme_pseudo_wildtype.pdb -extra_res_fa input/R1A.params -out:file:fullatom -resfile create_mtssl_mutant/config/resfile.pack -out:prefix mtssl_mutant_ >& make_mutant.log 
```
where (`$ROSETTA3`=path-to-Rosetta/main/source)

Step 2: relax_mtssl_mutant

The standard relax protocols can be used. The mtssl parameter file must be provided. The command line below will relax the t4-lysozyme double mutant provided in the input directory.
```
$> $ROSETTA3/bin/relax.default.linuxgccrelease -in:file:s relax_mtssl_mutant/input/lysozyme_start_mtssl_mutant.pdb -extra_res_fa input/R1A.params -out:file:fullatom  -out:prefix relax_ -nstruct 1 >& relax.log 
```

Step3: relax_mtssl_mutant_membrane

The standard relaxation protocols with membrane flags can be used. The mtssl parameter file needs to be provided. The command line below will relax the double mutant MSBA structure provided in the inputs directory.
```
$> $ROSETTA3/bin/relax.default.linuxgccrelease -in:file:s relax_mtssl_mutant_membrane/input/msba_mtssl_mutant_start_structure.pdb -out:file:fullatom -out:prefix mem_relax_ -extra_res_fa input/R1A.params -nstruct 1 -membrane -membrane:normal_cycles 100 -membrane:normal_mag 15 -membrane:center_mag 2 -file:spanfile relax_mtssl_mutant_membrane/input/msba.span >& mem_relax.log 
```
Step 4: rotamer_conformation_recovery

Step 5: epr_distance_distribution_agreement

Step 6: calculate_cone_model_parameters

Are for analysis.
