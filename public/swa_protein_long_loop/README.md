# SWA Protein Loop Long

KEYWORDS: LOOPS GENERAL

## Author
Rhiju Das, rhiju@stanford.edu

# Loop Remodeling by Enumeration: The Core Step of Protein 'Stepwise Assembly'

## Brief Description

This demo is an expansion of the demo swa_protein_main/. This demo involves a more detailed workflow for loops longer than 4-5 residues. It consists of buildup of the loop in 1-2 residue segments from either end, followed by chain closure. The workflow is described as a directed acyclic graph (DAG) of Rosetta jobs with well-defined dependencies.

All possible buildup-paths are followed through a dynamic-programming-like recursion, and up to a 1000 models are retained at each intermediate buildup path. The results are high in accuracy and low in energy, but requires significant computational expense (1000s of CPU-hours) and the ability to run a complex DAG on a cluster.

## Abstract

Consistently predicting protein structure at atomic resolution from sequence alone remains an unsolved problem in computational biophysics. Practical challenges involving protein loops arise frequently in ab initio modeling, comparative modeling, and protein design, but even these cases can become intractable as loop lengths exceed 10 residues and if surrounding side-chain conformations are erased. This demo illustrates a novel approach to protein modeling that is more powerful than prior methods that strive for atomic resolution. The central innovation is a ‘stepwise ansatz’ inspired by recent ab initio RNA algorithms, which resolves a conformational sampling bottleneck through residue-by-residue conformer enumeration and dynamic programming.


Reference: R. Das (2013) "Atomic-accuracy prediction of protein loop structures enabled by an RNA-inspired ansatz", under review.
More info: http://arxiv.org/abs/1208.2680

## Running

### Setup
You need to define an environment variable $ROSETTA with your Rosetta directory. Add to your .bashrc or .bash_profile a line like:

```
export ROSETTA='/Users/rhiju/src/rosetta/'   [change to your Rosetta directory]
```
 
You also need your system to know where the python scripts are for generating the DAG and running the jobs:

```
PATH=$PATH:$ROSETTA/rosetta_tools/SWA_protein_python/generate_dag/
PATH=$PATH:$ROSETTA/rosetta_tools/SWA_protein_python/run_dag_on_cluster/
```

### Example Python Command Line to Generate DAG
**Note: needs to be run in a copy of rosetta_inputs/**

This rebuilds residues 3-8 on the knottin scaffold 2it7 which has had the loop and all the protein's sidechains removed:

```
generate_swa_protein_dag.py  -loop_start_pdb noloop_2it7_stripsidechain.pdb  -native 2it7.pdb -fasta 2it7.fasta -cluster_radius 0.25 -final_number 1000   -denovo 1   -disulfide_file 2it7.disulf  -loop_res 3 4 5 6 7 8
```

If you want to do a quick run to test the overall workflow, you can keep the models within RMSD of 1.0 A to the experimental loop, use coarser backbone sampling, and save only 10 structures per run:

```
generate_swa_protein_dag.py  -loop_start_pdb noloop_2it7_stripsidechain.pdb  -native 2it7.pdb -fasta 2it7.fasta -cluster_radius 0.25 -final_number 1000   -denovo 1   -disulfide_file 2it7.disulf  -loop_res 3 4 5 6 7 8 -n_sample 9 -rmsd_screen 1.0 -nstruct 10
``` 

The outputs are:

```
  protein_build.dag
    text file outlining all the jobs, preprocessing and preprocessing script commands, and their dependencies. In condor DAGMAN format.

  CONDOR/
    directory with job definitions, in format similar to condor job definition format.

  REGION_3_2/, REGION_4_2/, ...
    directories that will hold outputs of each rosetta job. The numbers correspond to the N-terminal residue of the loop fragment reaching from the C-terminus endpoint of the loop, and the C-terminal residue of the loop fragment reaching from the N-terminal takeoff of the loop.  REGION_3_2 corresponds to models in which the loop fragments have met at the boundary between residues 2 and 3.
```

## How to Run the DAG

On condor clusters:

```
condor_submit_dag protein_build.dag
```

We have found condor_dagman can be slow due to latency in queuing rosetta jobs via Condor, unfortunately. A further problem is that there is currently no good universal solution to running DAGs on clusters, although packages like Pegasus and newer versions of Hadoop appear promising. 

For our own purposes, we have developed in-house python scripts (available in `rosetta_tools/SWA_protein_python/run_dag_on_cluster/`) to run the jobs by kicking off a master node that can queue jobs to slave nodes. Most recently, we have been using PBS/torque ('qsub') queueing, and you can run the jobs using 100 cores with the command:

```
 qsub README_QSUB 
```

which will queue from a master node:

```
 SWA_pseudo_dagman_continuous.py -j 100 protein_build.dag  > SWA_pseudo_dagman_continuous.out 2> SWA_pseudo_dagman_continuous.err
```

Further development for LSF clusters and MPI queuing is also under way. Please contact rhiju [at] stanford.edu with questions, or suggestions for supporting more general queuing systems.

**Note:** These scripts should work directly out of trunk for any version of Rosetta after March 2013.


