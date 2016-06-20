# Tips for using Rosetta
## Molecule size
This size of a macromolecule constitutes a major restriction on what you can and can not do with Rosetta (or any other simulation software). There are two aspects where the number of residues in a structure matter:  

1. The **time** for calculating energies increases with size, and make a protocol take (much) longer.  
2. The **degrees of freedom** (DOF) increase (often exponentially) with size. Exponential scaling is typical for combinatorial problems like folding simulations.  

The exact size limit depends on several factors:  

* **Resources:** Sometimes, the size limit depends on your computational resources; more DOFs require more sampling, which is fine if you have a lot of CPUs available.  
* **Type of simulation:**  Optimizing side chain conformations is much faster than predicting a strcture *de novo*. Hence, the size limits depends on your exact task and should be mentioned in the detailed documentation for a particular application.

###### Examples: 
* *Abinitio*: max 150 amino acids are cosidered possible
* *Loop modeling*: 12 amino acids is already considered quite long for giving trustworthy results
* *Docking*: As long as the backbone conformations are unchanged, there is *almost* no limit to the size of the individual subunits. 

##Sampling
Sampling means changing DOFs. The more DOFs need to be sampled, the more runs you will have to perform to get a reliable answer. There are two basic ways to increase sampling:

1. **-nstruct \< n \>** This option lets you to increase how many runs are performed by an individual CPU.
2. **Run on multiple CPUs** By starting the same simulation ( a.k.a. job ) on multiple CPUs, you can get more models in shorter time. Often, you will want to have a large nstruct *and* use multiple CPUs.  
**Note:** *You need to provide -out:prefix or -out:suffix options. Otherwise Rosetta would try to overwrite files produced by a parallel job -- in which case it will not do anything!*  


### Parallel runs using MPI 
To avoid the problem of file naming and for not having to start multiple jobs with different output file names, you can also compile Rosetta for usage with MPI (see "**Compiling Rosetta**"). **MPI** (short for Message Passing Interface) executables can launch many parallel jobs and handle file naming for you.  
However, the MPI version is intended for use on supercomputer clusters. There you can specify the number of processors for a job. Rosetta will then divide your *-nstruct* value by the number of CPUs.

Example: A snipped from a cluster submit file:

    #SBATCH --ntasks=4
    srun relax.mpi.linuxgccrelease -s input.pdb -nstruct 100
    
This will perform 25 relax runs on each of 4 CPUs. This is 4 x faster then *-nstruct 100* on a single CPU.


## Combine output models in single "silent files"
Instead of producing a large number of pdb files, Rosetta can provide the strctural and energetic information for any number of structures into a so-called silent file. Additionally, these silent files can be in binary format, which reduces the required disc space. Rosetta provides tools for working with these files. Look at the "Analysis" section for information of how to work with silent files. A couple of applications, e.g. scoring, are able to handle silent files directly as input. 


