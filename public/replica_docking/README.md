# Replica Docking

KEYWORDS: DOCKING INTERFACES

## Authors
This file was written in Jan 2013 by Zhe Zhang (zhe.zhang@tum.de) and corresponding PI is Oliver Lange (oliver.lange@tum.de)


## General Description
This demo contains all the files neccessary to replicate the results from the PLoS ONE RosettaCon
collection paper "Replica Exchange drastically improves sampling in low resolution docking stage
of RosettaDock" by Zhe Zhang, and Oliver Lange (2012).

All files, including the benchmark (`tar -xf dock_targetlib.tar.gz`) tested in the paper, and commands to setup the runs as well
 as post analysis are provided. For example outputs please refer to the silent-files (also trial.stat
from replica exchange runs) in example\_runs

## How to Run Demo
To run these demos:

1. export the directory where this protocol\_capture is located
    ```
    export PROTOCOL_CAPTURE=rosetta/rosetta_demos/protocol_capture/2012/
    ```

2. export your rosetta bin and database directory if you want to directly repeat the runs under directory example\_runs/ (need to first delete the outputted silent-files) instead of generating separate runs of your own. But this would only work for linux with slurm.
    ```
    export ROSETTA3_BIN=rosetta/rosetta_source/bin
    export ROSETTA3_DB=rosetta/rosetta_database
    ```

3. input pdbs preparation

   	all the targets are from Dockground benchmark3.0 (http://dockground.bioinformatics.ku.edu/UNBOUND/request\_new.php), in which unbound docking partners have been superimposed over its corresponding complex.
    ```
	get_pdb.py 1bvn_u1.pdb A	# this should output a file with only atom records from chain A, 1bvn_u1.pdbA.pdb
	replace_chain.py 1bvn_u2.pdb B > 1bvn_u2_B.pdb	# overwrite its chainID to B, 1bvn_u2_B.pdb
	get_pdb.py 1bvn_u2_B.pdb B	# 1bvn_u2_B.pdbB.pdb

	cat 1bvn_u1.pdbA.pdb 1bvn_u2_B.pdbB.pdb > protAB.pdb	# superimposed native structure, used for rmsd related calculation.
	scripts/initial_randomize.sh	# output P.pdb, both docking partners are randomly reoriented, and then slided into contact.
					# used as the input pdb for docking.
    ```
	write flag "-partners A_B" into file "partners" for later use.
    ```
	scripts/get_disulf_pairs.sh	# get the disulfide residue pairs, later used in refinement
    ```
	All the input files in dock_targetslib (`tar -xf dock_targetlib.tar.gz`)  are prepared in this way.

4. setup target library using automated setup tools available with the CS-Rosetta toolbox(www.csrosetta.org).

	This step assembles target related input files to build the target library (`dock_targetlib` as an example). By default 
    the library is stored in folder cs\_targetlib at home of your workspace. You can also specify a directory using the flag 
    `-target_prefix` as follows. Absolute path is recommended for `-target_prefix`.

	The following commands for setup demo target 1bvn have been wraped up in `scripts/test_target_setup.sh`
    ```
	#RosettaDock:
	#setup target for RosettaDock as in published paper "Protein-Protein Docking with Simutaneous Optimization of Rigid-body 
	#Displacement and Side-chain Conformations", Jeffrey J. Gray et al., J. Mol. Biol. (2003)
	setup_target -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
    -disulf disulf_file -native protAB.pdb -pdb P.pdb -partners partners

	#ReplicaDock:
	#setup target for ReplicaDock as described in Zhang and Lange, PLOS One 2012.
	setup_target -method replica_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
    -native protAB.pdb -pdb P.pdb -partners partners

    #or you can conviently copy the inputs from a previously prepared 'rosetta_dock' setup as follows:
    setup_target -method replica_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
    -transfer_method rosetta_dock
    ```

5. setup run using automated setup tools. This step creates a run-ready directory as specified with flag `-dir`, in which job-scripts, 
    input files, RosettaScripts-xml as well as flag files are contained. For flag `-dir`, absolute path is recommended. For job scripts, 
    you can use different types (e.g. moab) according to your queuing system.

	In example\_runs additional comments have been added subsequently to guide you through the automatically generated input files. 
    All job scripts in example\_runs are for queuing system with slurm.

	The commands to setup run as follows for single-mache/interactive mode are wraped up in `scripts/test_run_setup.sh`. After run 
    `scripts/test_run_setup.sh`, use `source production.interactive.job -n $Np` with $Np specifying the processor numbers to start the run 
    in the run directory, (for example `$PROTOCOL_CAPTURE/test_runs/replica_dock/udock_1bvn/run` after you run `scripts/test_run_setup.sh`)

	For quick test only purpose for ReplicaDock/ReplicaDock-LoT, I have set the total trajectory length to be 5000 MC-steps in 
    RosettaScripts(dock_cen.xml, dock_cen_7.xml) in the automated setup tool box, which can be simply found in 
    `csrosetta3/flag_library/methods/_docking_base/` and modified accordingly for production purpose use.

    1. centroid stage
	    1. ReplicaDock, using temperatures [2.0 3.0 5.0]

	        ReplicaDock is run in MPI-mode using RosettaScript to compile: 
            ```
            ./scons.py -j 48 bin/rosetta_scripts.mpi.linuxgccrelease mode=release extras=mpi
            ```

	        Please note that specific numbers of processors have to be used: calculate number of processes using the formula: nstruct \* n_replica + 2. 
            The extra 2 processes are dedicated to the job distributor and File IO. n_replica is the number of temperature levels (here 3), and nstruct 
            can be any positive integer. ReplicaDock outputs the trajectory in the form of two silent-files: one containing decoy+score information 
            (name: decoys\_<input_pdb>\_nnnn\_traj.out), and the second file is a copy of just the score information (scores\_<input_pdb>\_nnnn\_traj.out). 
            Decoytags are of the form P\_tttt\_rrr\_ssssssss where tttt informs about trajctory number, rrr about the replica, and ssssssss about the 
            snapshot number within the trajectory. The temperature\_levels are switched between different replicas. The current temp\_level or temperature
            of a replica at the moment a decoy was recorded is found in the score-columns temp\_level and temperature.
            At the end of a trajectory the final decoy is written to the file 'decoys.out'; this file is a relict of using the JD2-framework and
            can be generally ignored.
            Additionally, the file 'trial.stat' is produced which gives information about acceptance rates in each temperature level.
            ```
            setup_run -method replica_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
            -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/replica_dock -job interactive -extras mpi -score interchain_cen \
            -nstruct 1 -protocol rep_cen -xml uniform -n_replica 3
            ```

            start running:
            ```
            cd $PROTOCOL_CAPTURE/replica_docking/test_runs/replica_dock/udock_1bvn/run/;
            source production.interactive.job -n 5
            ```

	    2. ReplicaDock-LoT, using temperatures [0.6 0.8 1.0 1.2 1.5 2.0 2.5] + min_score

	        Min_score is used to flatten the score function. A reasonble min_score value is determined as the average score of the first 50 
            snapshots of temperature 1.5 when simulated without min-score. As shown in example_runs/replica_dock_LoT/get_min_score/udock_1bvn/run/, 
            several trajectories are run and an average value of -36.519 can be get by:
            ```
            cat scores_P_000* >scores_traj.fsc
            silent_data.py scores_traj.fsc temperature score | awk '$1==1.5{print}' | median.py

            setup_run -method replica_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
            -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/replica_dock_LoT -job interactive -extras mpi -score interchain_cen \
            -min_score -36.519 -nstruct 1 -protocol rep_cen -xml uniform -n_replica 7
            ```

	    3. RosettaDock's original low-resolution stage (shotgun sampling)
            ```
	        setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
            -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/rosetta_dock -job interactive -extras mpi -protocol centroid -batches 2 \
            -score interchain_cen -nstruct 25
            ```

    2. refinement. To refine the decoys generated in the centroid stage, we don't want to copy all the decoy-files into the run directory of refinement, but only specify the path of the decoys files. For this, we use the flag `-start` with absolute path specified together with flag `-pattern` to only include the files with a certain name pattern.

	    1. refinement of the low-resolution ensembles produced by ReplicaDock
        ```
	    setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
        -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/refine_replica_dock -job interactive -extras mpi -protocol refine \
        -pattern "low_decoys_*out" -prefix refine -score docking -nstruct 1 \
        -start $PROTOCOL_CAPTURE/replica_docking/example_runs/replica_dock/udock_1bvn/run/
        ```

	    2. refinement of the low-resolution ensembles produced by RosettaDock's shotgun sampling approach
        ```
	    setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
        -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/refine_rosetta_dock -job interactive -extras mpi -protocol refine \
        -pattern "low_decoys_*out" -prefix refine -score docking -nstruct 1 \
        -start $PROTOCOL_CAPTURE/replica_docking/example_runs/rosetta_dock/udock_1bvn/run
        ```

	    3. refinement of the low-resolution ensembles produced by ReplicaDock-LoT
        ```
	    setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
        -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/refine_replica_LoT -job interactive -extras mpi -protocol refine \
        -pattern "low_decoys_*out" -prefix refine -score docking -nstruct 1 \
        -start $PROTOCOL_CAPTURE/replica_docking/example_runs/replica_dock_LoT/replica_with_minscore/udock_1bvn/run
        ```

	    4. generate RelaxedNative ensembles with 1000 decoys from the superimposed native structure protAB.pdb, used as reference in final analysis
        ```
	    setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib \
        -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/relax_native -job interactive -extras mpi -protocol refine \
        -out relax_native.out -score docking -nstruct 1000
        ```

6. post-filter of low resolution ensembles. To save space, some .out files have been deleted and only score files shown.

	1. RosettaDock's shotgun approach sampled ensembles:
	use N to denote the total number of decoys produced from RosettaDock low resolution phase. Then we first exclude the decoys with interchain\_contact>10, then select N\*0.36 decoys by score.

	For an example decoys file in `$PROTOCOL_CAPTURE/replica_docking/example_runs/rosetta_dock/udock_1bvn/run/decoys.out`, we do:
    ```
	# get the tags of the decoys selected. 50 decoys in all in the example file, so select 50*0.36=18 decoys 
	# for analysis and further refinement
	for i in $(ls decoys_000?.out); do echo $i; cat $i | grep SCORE: >> decoys.fsc; done

	scripts/silent_data.py decoys.fsc score interchain_contact description | awk '$2<=10{print}' | sort -n -k 1 | head -n 18 \
    | awk '{print $3}' > tag_low

	# extract the selected decoys from decoys.out
	for i in $(ls decoys_000?.out); do echo $i; scripts/extract_tagged_decoys.py $i tag_low > low_$i; done
    ```

	2. replica_dock: T=[2.0 3.0 5.0]
	exclude decoys with temperature 5.0, then exclude decoys with interchain\_contact > 10
    ```
	# get the tags of the snapshots to refine
	cat scores_P_000*fsc > scores_traj.fsc
	scripts/silent_data.py scores_traj.fsc temperature interchain_contact description | awk '$1<4&&$2<=10{print $3}' > tag_low

	# extract selected decoys from the trajectory silent file
	for i in $(ls decoys_P_000?_traj.out); do echo $i; scripts/extract_tagged_decoys.py $i tag_low > low_$i; done
    ```

	3. replica_dock_LoT: T=[0.6 0.8 1.0 1.2 1.5 2.0 2.5] and min_score
	first exclude decoys with interchain\_contact > 10, then keep the 0.5\*N lowest scoring decoys from the remaining set.
    For the example silent-files of trajectory-snapshotse in `$PROTOCOL_CAPTURE/replica_docking/example_runs/replica_dock_LoT/replica_with_minscore/udock_1bvn/run/`, it is done as follows:
    ```
	# get tag. 14 decoys in all for the
	cat scores_* > scores_traj.fsc
	scripts/silent_data.py scores_traj.fsc interchain_contact score description | awk '$1<=10{print}' | sort -n -k 2 \
    | head -n 7 | awk '{print $3}' > tag_low

	# extract selected decoys from the trajectory silent file
	for i in $(ls decoys_P_000?_traj.out); do echo $i; scripts/extract_tagged_decoys.py $i tag_low > low_$i; done
    ```

7. refinement of low resolution ensembles

	refinement is setup with the automated tool box as shown in part 5.2.1-5.2.4.

8. analysis

	scripts for quick checking and analysis are given, i.e. `scripts/outfile_plot.py` and `scripts/hist.py`. To run these two scripts, extra python library matplotlib.pyplot is required.
    ```
	scripts/outfile_plot.py #generates scatter plots for given a silent_file with specified score terms, e.g. rms vs. score
	scripts/hist.py #generates histogram figure for given silent_files with specified score term
    ```
	
