Canonical Sampling for Protein-Protein Docking Refinement
=========================================================

KEYWORDS: DOCKING STRUCTURE_PREDICTION

Author: Zhe Zhang (zhezhang1986 at gmail dot com)  
Corresponding PI: Martin Zacharias (martin.zacharias at ph dot tum dot de)  
Last Updated: 06/01/2015  

Reference: Zhang Z, Schindler C, Lange OF, Zacharias M (2015): Comparison of
Replica-Exchange Approaches for Protein-Protein Docking Refinement in Rosetta

---

The different protocols described in this paper are tested on unbound docking
targets selected from docking benchmark4.0. In docking refinement practice,
approximate interaction site is known. Thus the initial start conformation is
generated thus from the superimposed unbound structure, by first translating
the second binding partner 15Å, then rotating 60°. This gives the initial
Ligand RMSD approximately sqrt( 15*15 + (2*pi*r*60/2/360)*(2*pi*r*60/2/360) ),
which r denoting the radii of the second binding partner. The translation
direction and rotation axis are both drawn from uniformly distributed vectors
on unit sphere.

Generating an initial conformation
----------------------------------

### Executable/Script:

    Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease

### PDB:

1. Download cleaned up pdb files from http://zlab.umassmed.edu/benchmark/, or 
   directly from http://www.rcsb.org and clean with clean_pdb.py in folder 
   Rosetta/tools/protein_tools/scripts/
2. Replace chain id if neccessary using scripts/replace_chain.py
       $ cat 1PPE_r_b.pdb 1PPE_l_b.pdb > native.pdb
       $ scripts/replace_chain.py 1PPE_r_u.pdb A > 1PPE_r_u.pdbA.pdb
       $ scripts/replace_chain.py 1PPE_l_u.pdb B > 1PPE_l_u.pdbB.pdb
       $ cat 1PPE_r_u.pdbA.pdb 1PPE_l_u.pdbB.pdb > protAB.pdb

### Running the Application:

Before run the test, please export your rosetta bin and database directory, 
as well as executable of mpirun

    $ export MPI_RUN=$YOUR_MPIRUN_EXECUTABLE
    $ export ROSETTA_BIN=$YOUR_ROSETTA_BIN_DIRECTORY
    $ export ROSETTA_DATABASE=$YOUR_ROSETTA_DATABASE_DIRECTORY

Run the command.sh script provided in folder generate_initial_conformation

    $ ./command.sh

### Example Outputs:

Example outputs are found in the folder generate_initial_conformation.  

1. P.pdb - final decoy, will be used as the initial starting conformation for all
the tests 
2. score.sc - score file

Common rules in this work
-------------------------

For rigid-body docking refinement, we have applied rigid-body mover
UnbiasedRigidBodyPerturbNoCenterMover and sidechain movers including
PerturbChiSidechainMover, PerturbRotamerSidechainMover and
JumpRotamerSidechainMover with the general Metropolis-Hastings framework. The
acceptance of a move is decided by the Metropolis Criterion. In order to avoid
the two binding partners diffuse away from each other in Monte-Carlo move,
 a very loose encounter constraint is applied and acts on the distance of the
mass center of the two binding partners. To achieve local search for docking
refinement, the rigid-body space is restricted with respect to the initial
conformation by translation of 20Å and rotation of 90° in
UnbiasedRigidBodyPerturbNoCenterMover. 

For productive simulation, 2,000,000 Monte-Carlo steps need to be run, and 
snapshots are stored every 1,000 steps. At the end of the simulation, only
decoys generated with the reference setting ( hardrep, temperature=0.15) are
collected and analyzed as final results. To save space and computer time,
the example output in this folder are from much shorter simulation. For 
productive simulations, please also change the related parameter "trial" in
dock.xml file.

Details of each protocol please also refer to the rosetta scripts file
dock.xml in each folder.

All the four protocols need to be run with MPI. 

### Example Outputs

Take the example outputs in wte_remc_docking for example:

1. decoys_P_0001_rt.out - silent trajectory file with rotation matrices and 
   translation vectors started with “RT”, and scores
2. decoys_P_0001_traj.out - silent trajectory file with decoys and scores
3. scores.fsc - silent score file of the trajectory
4. decoys.out - final decoy of a trajectory; this file is a relict of using the 
   JD2-framework and can be generally ignored. 
5. trial.stats - acceptance rate for each mover in each replica
6. tempering.stats - exchange rate between replicas
7. we_bias.grid - well-tempered ensemble bias information at each replica, 
   including grid size, grid range, bias energy in each bin, number of 
   conformations dropped into each bin. Only exist when BiasEnergy is applied, 
   for example in wte_remc_docking and wte_h_remc_docking.
8. decoys_P_0001_m_n.out - checkpoint silent decoy files with m indicating 
   replica number and n indicating the checkpoint number, used for restarting 
   the simulation. When BiasEnergy is applied, in the checkpoint silent file, 
   WTE bias energy information is stored as REMARK started with `REMARK 
   BIASENERGY`

### Analysis

    $ ref=5 # the number of reference replica
    $ scripts/silent_data.py decoys_P_0001_rt.out Lrmsd Irms score I_sc temp_level Fnat_n bias | awk -v ref_rep=“$ref” ’$5==ref_rep’ > collected_data
    $ scripts/collect_tempering_stats.py tempering.stats # collect the average exchange rate over the whole simulation
    $ scripts/collect_trial_stats.py trial.stats # collect average acceptance for each mover at each replica over the whole simulation

MC Docking
----------

This protocol using standard monte-carlo sampling protein-protein docking
refinement. FixedTemperatureController with temperature 0.15 in Rosetta is
used with the MetropolisHastings framework. The magnitude of the step size and
sampling weight of the rigid-body and sidechain movers are fixed along the
entire simulation. Rigid-body mover has a much lower sampling weight than the
sidechain movers to serve the purpose of refinement. 

### Executable/Script

Rosetta/main/source/bin/rosetta_scripts.mpi.linuxgccrelease

### Running the Application

Run the command.sh script provided in folder mc_docking: 

    $ ./command.sh -n $N_PROC   # N_PROC should be (2 + nstruct)

REMC Docking
------------

This protocol using parallel tempering sampling protein-protein docking
refinement. Temperatures are drawn from geometric progression withe the lowest
temperature same as used in mc-docking. HamiltonianExchange is used to control
the temperatures of each replica with the Metropolis-Hastings framework.
Exchange is attempted between neighbor temperatures every 1,000 steps. The
magnitude of the step size and the sampling weight all the movers are
modulated according to the temperature in the initialization such that in the
lower levels more frequent sidechain moves and few small rigid-body moves are
applied and in higher levels less frequent sidechain moves and more bigger
rigid-body moves are applied. 

### Executable/Script

Rosetta/main/source/bin/rosetta_scripts.mpi.linuxgccrelease

### Running the Application

Run the command.sh script provided in folder remc_docking: 

    $ ./command.sh -n $N_PROC   # N_PROC should be (2 + nstruct * n_replica)

WTE REMC Docking
----------------

This protocol applied well-tempered ensemble technique with parallel tempering
to sampling protein-protein docking refinement. By increasing the tunable
factor gamma of BiasedEnergy, we can reduce the number of replicas, but
maintain an approximately same exchange rate between replicas. 

### Executable/Script

    Rosetta/main/source/bin/rosetta_scripts.mpi.linuxgccrelease

### Running the Application

Run the command.sh script provided in folder wte_remc_docking: 

    $ ./command.sh -n $N_PROC   # N_PROC should be (2 + nstruct * n_replica)

WTE H REMC Docking
------------------

In this protocol, two dimensional replica exchange, with variable of the first
dimension as temperature, and of the second dimension as the scaling of
softness of repulsive Lennard-Jones potential. In this second scaling
dimension, we have tested with five levels: standard (hard rep), soft50%,
soft55%, soft60% and soft65%. In the dimension with variable of temperature,
we have five temperatures with lowest equal 0.15. In total 25 levels are run
in parallel and exchange between neighbor levels is attempted every 1,000
steps periodically along the two dimensions.

### Executable/Script

    Rosetta/main/source/bin/rosetta_scripts.mpi.linuxgccrelease

## Running the Application

Run the command.sh script provided in folder wte_remc_docking: 

    $ ./command.sh -n $N_PROC   # N_PROC should be (2 + nstruct * n_replica)

