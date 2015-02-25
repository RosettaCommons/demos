/path/to/Rosetta/main/source/bin/docking_prepack_protocol.macosclangrelease \
-database /path/to/Rosetta/main/database \			#path to Rosetta database
-in:file:s input/1AFO_AB.pdb \						#input PDB
-mp:setup:spanfiles input/1AFO_AB.span \				#input spanfile
-out:file:scorefile score_ppk_1AFO.sc \				#name of the output scorefile
-nstruct 1 \											#number of models to build
-score:weights mpframework_smooth_fa_2012.wts \		#weights file for scoring
-packing:pack_missing_sidechains 0 \					#don't pack sidechains until membrane is set up
>& output/log_ppk_1AFO.log &							#run in background and log the output
