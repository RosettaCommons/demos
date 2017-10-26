//////////////////////////////////////////////////////////////////////////////////////////
// (c) Copyright Rosetta Commons Member Institutions.
// (c) All the files in this directory and sub-directories are part of the
// (c) Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
//  @brief
//  This is a tutorial set for a de novo model building method described in Wang et. al.,
//  entitled "Accurate de Novo Protein Structure Determination from Near-atomic
//  resolution Cryo-EM Maps". In this tutorial, I use one of the benchmark cases from the
//  manuscript, TRPV1, to demonstrate the method.
//
//  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
//  @author Frank DiMaio, dimaio@u.washington.edu
//
//////////////////////////////////////////////////////////////////////////////////////////

This tutorial includes three folders:
@scripts/
    Wrappers for calling rosetta executables, setting up jobs, and results processing.

@input_files/
    Files necessary for running the protocol on an example case, the 3.4A TRPV1 cryo-EM
    map (the exact commands used in the manuscript).

@denovo_model_building/
    The directory where jobs will be run. Calculations are split into several steps, with
    each step run in an individual subfolder:

    denovo_model_building/
        # round1
        Step1_Place_fragments_into_density/
        Step2.1_Calculate_overlap_scores/
        Step2.2_Calculate_nonoverlap_scores/
        Step3_Simulated_annealing_Monte_Carlo_sampling/

        # round2
        Step4_remove_density_using_average_model/
        Step5_round2_Place_fragments_into_density/
        Step6.1_round2_Calculate_overlap_scores/
        Step6.2_round2_Calculate_nonoverlap_scores/
        Step7_round2_Simulated_annealing_Monte_Carlo_sampling/
        StepFinal_RosettaCM_to_complete_the_model/


##########################################################################################
STEP 0: SETUP

Edit the configuration file, "denovo_model_building_scripts.cfg".
    demo_dir = absolute_path_for_this_demo_in_your_work_station


##########################################################################################
STEP 1: Dock fragments into density

In this step, we dock local backbone fragments into the EM density.  Fragments for
a perticular protein may be downloaded from Robetta (link).  For this case, the
necessary fragment files have been included.

NOTE: As this step is particularly time consuming, it is recommended to run on a
moderately-large cluster (64+ compute cores).  Submission scripts for the condor
queuing system are included; it should be relatively straightforward to convert these
scripts for other queuing systems.

Go to: denovo_model_building/Step1_Place_fragments_into_density/condor_jobs

    A.  Rotation/translation search in a map for each fragment is computationally
        expensive. The search is split at the residue level to enable parallelism.

        Go to the directory:
            denovo_model_building/Step1_Place_fragments_into_density/condor_jobs/

		And run the command:
			./setup_placement_condor_jobs.sh

		This will setup a single condor job for each residue in the protein (in this
		example, 307 jobs).

		Then, to submit the jobs, simply run:
			sh submit_placement_condor_jobs.sh

        For running on a different system, this file will need to be updated to point
        to the new sequence, map, and fragment files.

    B.  After jobs finish, run a script to ensure all output has been produced.  Go to
		the folder:
            denovo_model_building/Step1_Place_fragments_into_density/condor_jobs/

		And run the command:
            ./find_out_unfinished_placement_jobs.sh  | grep -v input_files

        This script will print out all residues that have not finished running (or
		died unexpectedly).  If you need to submit one of these jobs, go to the relevent
        folder, remove "running.lock" and run "condor_submit placement_condor_job".
		If interrupted, the code will resume from wherre it left off.

	C . Once fragment placement is complete, we cluster and extract the best-scoring
        candidate placements.  Go to the directory:
            denovo_model_building/Step1_Place_fragments_into_density/condor_jobs/

		And run the command:
            ./setup_cluster_and_extract_condor_jobs.sh

		This will again setup a single condor job for each residue in the protein.  To
        launch this job, run the following command on the condor system:
			sh submit_cluster_and_extract_condor_jobs.sh

		Finally, once this step is complete, the output will appear in the following
        directory:
            denovo_model_building/Step1_Place_fragments_into_density/candidate_fragment_placements/

		The output from this step is 50 candidate placements for each residue. In this
        example, TRPV1, it should contain 50*307=15350 placements.

		NOTE: If the number of extracted fragments does not match, check if placement
        jobs haven't finished or if clustering failed, and rerun steps A or C.


##########################################################################################
STEP 2: Precompute fragment compatibility scores

In this step, we precompute the compatibility scores of all pairs of fragments identified
in the previous step.  This is subdivided into two jobs, the first computes the overlap
scores and the second the nonoverlap scores.

As with the previous step, this will launch a single compute job for each residue in the
protein.  Submission scripts for the condor queuing system are included; it should be
relatively straightforward to convert these scripts for other queuing systems.

    A. To calculate overlap scores, go to the directory:
            denovo_model_building/Step2.1_Calculate_overlap_scores/
       And run the command:
           ./setup_condor_jobs.sh

       This will setup a single condor job for each residue in the protein.  To
       launch this job, run the following command on the condor system:
			sh submit.sh


    B. For nonoverlap scores, the setup is the same.  Go to the folder:
            denovo_model_building/Step2.2_Calculate_nonoverlap_scores

       And run the command:
           ./setup_condor_jobs.sh

       Launch the condor jobs with:
			sh submit.sh


##########################################################################################
STEP 3: Run Monte Carlo sampling to identify maximally consistent subset of fragments

In this step, different combinations of fragments are combinded to identify the
maximally compatible subset of placements.  These are run in many parallel Monte Carlo
trajectories, in order to identify a mutually consistent subset.

    A. Set up score tables.  This step simply combines the results from steps 2A and 2B.
       Go to the directory:
            ./denovo_model_building/Step3_Simulated_annealing_Monte_Carlo_sampling/
       And run the command:
            ./setup_idx_files.sh

       Running this script will produce four output files, which serve as inputs for the
       next step of the protocol:
            @frags.idx1
                an index of fragments from candidate_fragment_placements/

            @all_density.idx1
                the density score for each fragment

            @all_nonoverlap_scores.weighted.idx1
                two-body score terms: closability score and clash score

            @all_overlap_scores.idx1
                two-body score terms: overlap score

    B. Run Monte Carlo sampling.  Go to the directory:
            ./denovo_model_building/Step3_Simulated_annealing_Monte_Carlo_sampling/
       And run the command:
            ./run_samc_sampling.sh [number_of_jobs]

       This step can be run on one machine, and will launch "number_of_job" threads
       on the node on which it is run.  This may be launched on several machines to
       generate additional trajectories (the outputs will not collide).

	C. Assemble a partial model.  This will combine converged regions from the low-
       energy trajectories, producing a high-confidence partial model.  Go to the
       directory:
            ./denovo_model_building/Step3_Simulated_annealing_Monte_Carlo_sampling/
       And run the command:
		    ./assemble_partial_model.sh

       The output from this step will be in the pdb file:
            denovo_model_building/Step3_Simulated_annealing_Monte_Carlo_sampling/average_model/average.pdb


##########################################################################################
STEP 4. Mask out density for parts of the protein already placed

If the model has at least 70% of residues placed in step three, then we can skip to the
final step.  Otherwise, we mask regions of the map corresponding to protein that has been
placed, and rerun steps 1 through 3.

We begin by masking density from protein that has already been placed.  Go to the
directory:
	denovo_model_building/Step4_remove_density_using_average_model/
And run the command:
	./run.sh


##########################################################################################
STEP 5. Rerunning Step 1 with reduced map

Here (steps 5-7), we again dock fragments into density. The major difference is step D,
where the results from step 1 are combined with the latest results.

    A.  Dock fragments.  Go to:
            denovo_model_building/Step5_round2_Place_fragments_into_density/condor_jobs/
		And run the command:
            ./setup_placement_condor_jobs.sh

        This sets up condor jobs.  Submit them by running the following on your condor
        cluster:
            sh submit_placement_condor_jobs.sh


    B.  After all the jobs finish, check to see if everything ran correctly.  Go to:
            denovo_model_building/Step5_round2_Place_fragments_into_density
		And run the command:
            ./find_out_unfinished_placement_jobs.sh  | grep -v input_files

        Anything output here needs to be rerun (See 1B for details).

	C . Cluster and extract the resulting placements.  Go to:
            denovo_model_building/Step5_round2_Place_fragments_into_density/condor_jobs/
		And run the command:
            ./setup_cluster_and_extract_condor_jobs.sh

        Submit to the condor cluster:
			sh submit_cluster_and_extract_condor_jobs.sh

    D. Finally, combine the results of round 1 and and round 2.  Go to:
            denovo_model_building/Step5_round2_Place_fragments_into_density/
       And run the command:
            ./prepare_round2.sh


##########################################################################################
STEP 6. Precompute fragment compatibility scores from combined fragment set

    A. Calculate overlap scores.  Go to:
            denovo_model_building/Step6.1_round2_Calculate_overlap_scores/
	   And run the command:
           ./setup_condor_jobs.sh

       Submit to the condor cluster:
			sh submit.sh

    B. Calculate nonoverlap scores.  Go to:
            denovo_model_building/Step6.2_round2_Calculate_nonoverlap_scores
	   And run the command:
           ./setup_condor_jobs.sh

       Submit to the condor cluster:
			sh submit.sh

##########################################################################################
STEP 7. Run Monte Carlo sampling on the combined fragment set

    A. Set up score tables.  Go to:
            ./denovo_model_building/Step7_round2_Simulated_annealing_Monte_Carlo_sampling
        run
            ./setup_idx_files.sh

        The outputs are the same as in step 3A.

    B. Run SAMC
        under:
            ./denovo_model_building/Step7_round2_Simulated_annealing_Monte_Carlo_sampling
	   And run the command:
            ./run_samc_sampling.sh [number_of_jobs]

	C. Assemble a partial model.  Go to:
            ./denovo_model_building/Step7_round2_Simulated_annealing_Monte_Carlo_sampling

	   And run the command:
            ./assemble_partial_model.sh

      The output from this step will be in the pdb file:
 			denovo_model_building/Step7_round2_Simulated_annealing_Monte_Carlo_sampling/average_model/average.pdb



##########################################################################################
STEP 8. Complete the partial model using RosettaCM

In this final step, we will complete the model using RosettaCM.  While the commands
below only create 10 models, in some cases (particularly if there are large gaps in the
partial model) it may be beneficial to create more.  Thus, it may be useful to submit
multiple jobs.

	A. Set up RosettaCM.  This sets up the input command line for RosettaCM.  Go to:
            ./denovo_model_building/StepFinal_RosettaCM_to_complete_the_model
	   And run the command:
            ./setup_rosetta_cm.sh

		If you want to run more jobs, launch them with the command:
            ../../rosetta/rosetta_scripts.static.linuxgccrelease @rosetta_cm_flags -mute all

		The output from RosettaCM will be in a compressed format, in the file:
        "rosetta_cm.out"

	B. Pick full-length models and extract them.  This script picks the best 10 models
       using energy and fit-to-density, and extracts them as PDBs:
		    ./pick_models_from_rosettacm.sh

       The convergence of these models also provides some measure of model confidence.

