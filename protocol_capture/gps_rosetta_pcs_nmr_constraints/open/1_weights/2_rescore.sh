#!/bin/sh                                                    

# flags                                                                                                                                             
################                                                                                                                                    
# Please change the executable paths according to your operating system                                                                             
# For arg1, give the silent file generated using 1_Rosetta_gps.sh
# For arg2, give the target protein starting model in pdb format                                                                                    
# For arg3, give input broker file                                                                                                                  
# For arg4, give the ouput silent file name                                                                                                         
# For arg5, score all the decoys using the unity PCS weights patched via 'pcsweight.patch' file
# Modify the pcsweight.patch file with the number of tag data that you have
# Protein system dengue virus protease (dvp)                                                                                                       

protein="dvp"                                                                                                                                       
Executable="/Users/julialeman/Documents/julia/git_final/Rosetta/main/source/bin/score.macosclangrelease"                                      
Database="-database /Users/julialeman/Documents/julia/git_final/Rosetta/main/database"                                                              
arg1="-in:file:silent ./$protein.silent"                                                                                                            
arg2="-in:file:native ../setup/ideal_model.pdb"                                                                                                     
arg3="-broker:setup ../setup/rigid_plus_pcs_broker.txt -run:protocol broker -overwrite"                                                             
arg4="-out:file:scorefile ${protein}_pcs.csc"                                                                                                       
arg5="-score:weights pcsweight.patch"                                                                                                             
                                                                                                                                                    
# optional flags                                                                                                                                    
################                                                                                                                                    
#replace mpi compiled executable path to run in parallel
# mpirun="mpirun -np 2"                                                                                                                             
                                                                                                                                                   
run="$Executable $Database $arg1 $arg2 $arg3 $arg4 $arg5"                                                                                 
echo $run                                                                                                                                           
$run                                                                                                                                                
       

