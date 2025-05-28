#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=32:00:00
#SBATCH --no-requeue
#SBATCH --job-name="drug_design"
#SBATCH --array=1-100
#SBATCH --output=slurm/%A_%a_slurm.out
#SBATCH --error=slurm/%A_%a_slurm_error.err
 
ID=$SLURM_ARRAY_TASK_ID

protein=<protein>
LIG=<ligand>

mkdir -p production/$LIG/
path/to/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease @drug_design_options.txt \
	-parser:protocol drug_design.xml \
	-s inputs/${protein}_${LIG}.pdb \
	-extra_res_mol inputs/${LIG}.sdf \
	-out:path:all production/$LIG/ \
	-out:prefix ${ID}_ \
	-nstruct 1
