/path/to/Rosetta/main/source/bin/mpdocking.macosclangrelease \
-database /path/to/Rosetta/main/database \			#path to Rosetta database
-in:file:s input/1AFO_AB_noMEM.pdb \					#input PDB
-in:file:native input/native.pdb \					#native for comparison for RMSD calculation
-mp:setup:spanfiles input/1AFO_AB.span \				#input spanfile
-out:file:scorefile score_mpdock_1AFO.sc \			#name of the output scorefile
-nstruct 1 \											#number of models to build, at least 1000 for production run
-docking:partners A_B \								#chains of the partners to dock
-dock_pert 3 8 \										#docking perturbation: 3A translation, 8degrees rotatin 
-show_simulation_in_pymol 0 \						#output every frame to pymol, take out for production run!!!
>& output/log_mpdock_1AFO.log &						#run in background and log the output
