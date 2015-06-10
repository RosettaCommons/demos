#!/bin/bash

# Variable flags
################
# Please change the executable paths according to your operating system.
# For arg1 and arg2, give input 9 and 3 residue fragment library files.
# For arg3, give the target protein starting model in pdb format.
# For arg4, give input broker file.
# For arg5, give number of structures to generate, typically 1,000 structures are needed. Atleast generate 100.
# For arg6, give the ouput silent file name.
# For arg7, give output score file name.
# Protein system dengue virus protease (dvp)                                                                                                       

protein="dvp"
Executable="/Users/julialeman/Documents/julia/git_final/Rosetta/main/source/bin/minirosetta.macosclangrelease"
Database="-database /Users/julialeman/Documents/julia/git_final/Rosetta/main/database"
arg1="-frag9 ../setup/frag9_$protein.tab"
arg2="-frag3 ../setup/frag3_$protein.tab"
arg3="-native ../setup/ideal_model.pdb"
arg4="-broker:setup ../setup/rigid_plus_pcs_broker.txt -run:protocol broker -overwrite"
arg5="-nstruct 5"
arg6="-out::file::silent $protein.silent"
arg7=" -out:file:scorefile ${protein}_allatom.fsc"

# PCS weights patch
###################

patch1="-abinitio::stage1_patch ../1_weights/pcs.wts"
patch2="-abinitio::stage2_patch ../1_weights/pcs.wts"
patch3="-abinitio::stage3a_patch ../1_weights/pcs.wts"
patch4="-abinitio::stage3b_patch ../1_weights/pcs.wts"
patch5="-abinitio::stage4_patch ../1_weights/pcs.wts"                                                                                                                                                      
patch_flag="$patch1 $patch2 $patch3 $patch4 $patch5"                                                                                                

# Default flags
###############
arg8="-abinitio::increase_cycles 1.0"
arg9=" -skip_stages 1 2 -seq_sep_stages 1 1 1 -short_frag_cycles 1  -ramp_chainbreaks -sep_switch_accelerate 0.8 -skip_convergence_check  -fail_on_bad_hbond false"

#Perform relax after abinitio
#############################
arg10="-abinitio:relax -relax:fast"

# optional flags
################
#replace the executable path with mpi compiled one to run in parallel
# mpirun="mpirun --hostfile ./mpi_hostfile -np 96"

run="$Executable $Database $arg1 $arg2 $arg3 $arg4 $arg5 $arg6 $arg7 $arg8 $arg9 $patch_flag $arg10"
echo $run
$run
