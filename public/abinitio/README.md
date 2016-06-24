#AbInitio: Predict structure from sequence 

KEYWORDS: STRUCTURE_PREDICTION GENERAL

Run on linux 
------------
(where `$ROSETTA3`=path-to-Rosetta/main/source)

    $> $ROSETTA3/bin/AbinitioRelax.default.linuxgccrelease @flags

Run on macs
-----------
    Rosetta/main/source/bin/AbinitioRelax.macosgccrelease @flags

Flags
-----
    -in:file:fasta ./input_files/1elwA.fasta
    -in:file:frag3 ./input_files/aa1elwA03_05.200_v1_3
    -in:file:frag9 ./input_files/aa1elwA09_05.200_v1_3
    -in:file:native ./input_files/1elw.pdb

    -abinitio:relax
    -nstruct 1
    -out:pdb

    -use_filters true
    -psipred_ss2 ./input_files/1elwA.psipred_ss2
    -abinitio::increase_cycles 10
    -abinitio::rg_reweight 0.5
    -abinitio::rsd_wt_helix 0.5
    -abinitio::rsd_wt_loop 0.5
    -relax::fast

Analyze Output
--------------
The output_files directory contains example output.

In `score.fsc` get a score and RMS for each model.
Typical analysis makes a scatter plot of these with RMS on the x-axis and score on the y-axis.
Look for a "funnel" to low energies and RMS in a successful ab initio prediction.
A failed prediction will not have low RMS/energy structures.
For the following arguments for full production run (using the minirosetta compile):

    -abinitio::fastrelax
    -abinitio::increase_cycles 10
    -abinitio::rg_reweight 0.5
    -abinitio::rsd_wt_helix 0.5
    -abinitio::rsd_wt_loop 0.5
    -abinitio::use_filters false
    -database minirosetta_database
    -ex1
    -ex2aro
    -frag3 aa0000103_05.200_v1_3
    -frag9 aa0000109_05.200_v1_3
    -in:file:boinc_wu_zip ploop23_control_fold_data.zip
    -in:file:native 00001.pdb
    -mute all
    -mute all
    -out:file:silent default.out
    -relax::default_repeats 15
    -silent_gz

    resultfiles = default.out.gz queue = 3000

And run a relax:

    -database
    -ex1
    -ex2aro
    -frag3 aa0000103_05.200_v1_3
    -frag9 aa0000109_05.200_v1_3
    -in:file:boinc_wu_zip ploop23_control_fold_data.zip 
    -in:file:fullatom
    -in:file:native 00001.pdb
    -in:file:s 00001.pdb
     minirosetta_database
    -mute all
    -out:file:silent default.out
    -relax::default_repeats 15
    -run:protocol relax
    -silent_gz

    resultfiles = default.out.gz
