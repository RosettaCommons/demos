
## Run on linux
../../bin/AbinitioRelax.linuxgccrelease @flags

## Run on macs
~/rosetta/rosetta_source/bin/AbinitioRelax.macosgccrelease @flags

## flags ##
-in:file:fasta ./input_files/1elwA.fasta
-in:file:frag3 ./input_files/aa1elwA03_05.200_v1_3
-in:file:frag9 ./input_files/aa1elwA09_05.200_v1_3
-database /Users/lucas/rosetta/rosetta_database/

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

####### analyze output #####
###### example output in output_files #####

In score.fsc get a score and rms for each model.
Typical analysis makes a scatter plot of these with rms on x-axis and score on y. 
Look for a "funnel" to low energies and rms in a successful ab initio prediction. A failed prediction will not have low rms/energy structures.

############################################################################
### for a full production run ### (using the minirosetta compile)
arguments =  -frag3 aa0000103_05.200_v1_3 -frag9 aa0000109_05.200_v1_3
-abinitio::increase_cycles 10 -mute all -abinitio::fastrelax
-relax::default_repeats 15 -abinitio::rg_reweight 0.5
-abinitio::rsd_wt_helix 0.5 -abinitio::rsd_wt_loop 0.5
-abinitio::use_filters false -ex1 -ex2aro  -in:file:boinc_wu_zip
ploop23_control_fold_data.zip   -database minirosetta_database
-out:file:silent default.out -silent_gz -mute all   -in:file:native
00001.pdb resultfiles = default.out.gz queue = 3000

And run a relax:
 -run:protocol relax -relax::default_repeats 15 -frag3 aa0000103_05.200_v1_3 -frag9 aa0000109_05.200_v1_3 -mute all    -ex1 -ex2aro  -in:file:boinc_wu_zip ploop23_control_fold_data.zip  -database minirosetta_database -out:file:silent default.out -silent_gz -in:file:native 00001.pdb -in:file:fullatom -in:file:s 00001.pdb resultfiles = default.out.gz
