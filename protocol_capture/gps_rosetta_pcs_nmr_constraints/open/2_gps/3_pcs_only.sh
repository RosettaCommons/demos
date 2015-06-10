#!/bin/sh                                                    

# flags                                                                                                                                             
################                                                                                                                                    
# Please change the executable paths according to your operating system                                                                             
# For arg1, give the silent file generated using 1_Rosetta_wts.sh
# For arg2, give the target protein starting model in pdb format                                                                                    
# For arg3, give input broker file                                                                                                                  
# For arg4, give the ouput silent file name                                                                                                         
# For arg5, PCS energy with unity weights
# For arg6, exclude unstructured flexible region from RMSD calculation
# Protein system dengue virus protease (dvp)                                                                                                       

protein="dvp"                                                                                                                                       
Executable="/Users/julialeman/Documents/julia/git_final/Rosetta/main/source/bin/score.macosclangrelease"                                      
Database="-database /Users/julialeman/Documents/julia/git_final/Rosetta/main/database"                                                              
arg1="-in:file:fullatom -in:file:silent ./$protein.silent"
arg2="-in:file:native ../setup/ideal_model.pdb"                                                                                                     
arg3="-broker:setup ../setup/rigid_plus_pcs_broker.txt -run:protocol broker -overwrite"                                                             
arg4="-out:file:scorefile ${protein}_pcs_only.fsc"
arg5="-score:weights ../1_weights/pcsweight.patch"
arg6="-native_exclude_res 1 2 3 4 5 6 7 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247"
                                                                                                                                                    
# optional flags                                                                                                                                    
################                                                                                                                                    
#replace mpi compiled executable path to run in parallel
# mpirun="mpirun -np 2"                                                                                                                             
                                                                                                                                                   
run="$Executable $Database $arg1 $arg2 $arg3 $arg4 $arg5"                                                                                 
echo $run                                                                                                                                           
$run                                                                                                                                                
       

