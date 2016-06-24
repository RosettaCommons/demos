Clustering: Group similar decoys with the Rosetta cluster application.
======================================================================

KEYWORDS: ANALYSIS GENERAL

Clustering is used to group output decoys.  The Rosetta clustering app is only typically good for a total decoy set of 400 or less, as it uses the first 400 structures to create the center of the clusters.   It also cannot handle thousands of decoys efficiently.

The app will eventually be deprecated.  If you have a large number of decoys or wish to use a [better] clustering app, please use the calibur application.

The following demo is intended as a good reference for those wishing to use the current Rosetta clustering application.

The following command is a good place to start:
	
    $> tar -xvzf rosetta_inputs/clustering_rosetta_inputs_demo.tgz -C rosetta_inputs

    $> $ROSETTA3/bin/cluster.default.linuxgccrelease  -in:file:s rosetta_inputs/complex1*.pdb -in:file:fullatom -cluster:gdtmm -cluster:radius -1 -cluster:population_weight 0.0 -sort_groups_by_energy 

where `$ROSETTA3`=path-to-Rosetta/main/source

These flags may be also be useful:

    -in:file:silent <filename.out>
    -in:file:silent_struct_type <desired type but binary is best choice>
    -limit_cluster_size <n>
    -limit_clusters <n>
    -cluster:radius <for gdt you want some number between 10-50>
    -out:file:silent clustered.out -out:file:silent_struct_type binary

You can extract the silent file by running the mini app `extract_pdbs.linuxgccrelease`.

If you have thousands of structures it may be best to first do an energy cut to filter out high energy structures.
The pdb's should be first scored and stored in an out file.
Try using the following app

    $> $ROSETTA3/bin/score_jd2.default.linuxgccrelease -in:file:s rosetta_inputs/complex1*.pdb -out:file:silent scored_silent.out -out:file:silent_struct_type binary -in:file:fullatom

After the pdbs have been scored use the included script to do an energy cut (This script can also be found in ```main/source/scripts/public/python```).

    $> extra_scripts/make_sub_silent_file_percentile.py scored_silent.out ecut_10.out -1 0.10

where 0.10 is a 10% energy cut.

 
