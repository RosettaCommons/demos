#Protein-Protein Docking 

KEYWORDS: DOCKING INTERFACES  

[[_TOC_]] 

**Bold text means that these files and/or this information is provided.**

*Italicized text means that this material will NOT be conducted during the workshop*

    fixed width text means you should type the command into your terminal

If you want to try making files that already exist (e.g., input files), write them to a different directory! (mkdir my_dir)

#Tutorial

This tutorial presents a cross-docking benchmark experiment. Antibody CR6261 binds to multiple sub- types of influenza antigen hemagglutinin (HA). It has been crystallized with H1 and H5 HA sub-types. Antibody from one crystal structure will be docked to the antigen from the other crystal structure. This type of experiment is useful for protocol optimization and development.

This tutorial was updated on 29 May 2017 by Vikram K. Mulligan (vmullig@uw.edu) to make it compatible with Rosetta's new ref2015 energy function.            

1. Prepare the input template for docking
    1. Download the PDB files. **The 3GBN.pdb and 3GBM.pdb files are provided in the input_files directory.**
        1. Download 3GBN from the Protein Data Bank.
            1. Go to rcsb.org and type '3gbn' in the search bar.
            1. Click on 'Download Files' on the right side of the page, then 'PDB File (Text)'.
            1. Save the PDB file in the my_files directory as `3GBN.pdb`.
        1. Repeat for '3GBM'

    1. Clean the PDBs. 

        1. We want the hemagglutinin (chains A and B) from 3GBM

                $> $ROSETTA_TOOLS/protein_tools/scripts/clean_pdb.py input_files/3GBM AB

                $> $ROSETTA_TOOLS/protein_tools/scripts/pdb_renumber.py --norestart input_files/3GBM_AB.pdb input_files/3gbm_HA.pdb

        1. We want the antibody (chains H and L) from 3GBN. (Note that chains A and B are not the antibody "Ab"). We only need the variable domain which is actually involved with binding HA. The crystal structure also contains a partially resolved portion of the constant domain. You should manually edit the PDB file with a text editor to remove the unnecessary portions. You should be able to see them in a structure viewer. It should be residues 121-160 of the heavy chain (chain H) and residues 268-311 of the light chain (chain L) in the cleaned structure.

                $> $ROSETTA_TOOLS/protein_tools/scripts/clean_pdb.py 3GBN HL

                cp 3GBN_HL.pdb 3GBN_trim.pdb
                pymol 3GBN_trim.pdb
                gedit 3GBN_trim.pdb

                $> $ROSETTA_TOOLS/protein_tools/scripts/pdb_renumber.py --norestart input_files/3GBN_trim.pdb input_files/3gbn_Ab.pdb

   -> Note: you may not have gedit in your computer. Try using other text editors.
    
    1. Close the chain break between Ser-127 and Val-128 in chain L of the antibody. 
        
        Chain breaks can cause unexpected behavior during docking. Because this chain break is small and away from the expected interface, it can be quickly and easily fixed. The goal is simply to close the chain break within the secondary structure element and not to rigorously build this loop. Build ten models (a minimal computational effort) and pick one with a good score and a good representative structure.
        
        1. Open 3gbn_Ab.pdb with PyMol to identify the chain break between residues 127 and 128 of chain L.
            
                pymol 3gbn_Ab.pdb

            type '`as cartoon`' or '`as ribbon`'
            
            type '`show lines, resi 125-130`'

        1. Prepare a loops file for closing the chain break. Use a text editor such as gedit. Several amino acids must be mobile for the loop to close successfully. Select several residues on each side of the chain break. (Note that you may not have gedit on your computer in which case you need to use a different text editor).

                    gedit chainbreak_fix.loops
            
            1. Type '`LOOP 125 130 0 0 1`', save the file, and close gedit.
            1. Copy the prepared options file from the input_files directory 
            
                    $> cp input_files/chainbreak_fix.options .
            1. If you don't have the prepared files, you can copy them into the directory
            ```
                    $> cp input_files/3gbn_Ab.pdb input_files/chainbreak_fix.loops .
            ```
        1. Run the Rosetta loopmodel application to close the loop.

                $> $ROSETTA3/bin/loopmodel.default.linuxgccrelease @chainbreak_fix.options 

                pymol 3gbn_Ab*pdb &
                > sort -nk 2 chainbreak_fix.fasc
        1. Copy the best scoring model (3gbn_Ab_0010.pdb in the example output_files directory) to 3gbn_Ab_fixed.pdb. Alternatively you can copy the provided structure into your directory:

                > cp 3gbn_Ab_0010.pdb 3gbn_Ab_fixed.pdb
                    or
                $> cp input_files/3gbn_Ab_fixed.pdb .

    1. Repack or relax the template structures.

        Repacking is often necessary to remove small clashes identified by the score function as present in the crystal structure. Certain amino acids within the HA interface are strictly conserved and their conformation has been shown to be critical for success in docking. RosettaScripts allows for fine control of these details using TaskOperations.
        
        1. Copy the XML [scripts]() and options file for repacking from the input_files directory.
        
                $> cp input_files/repack.xml .
                $> cp input_files/repack_HA.xml . 
                $> cp input_files/repack.options .
            1. make sure you have all the necessary files. If not, copy them from the input directory:
            ```
            $> cp input_files/3gbm_HA.pdb .
            ```
        1. Familiarize yourself with repack.xml and repack_HA.xml. Notice that repack_HA.xml is a modified version of repack.xml, representative of the versatility of RosettaScripts.
        1. Run the XML script with the rosetta_scripts application.
        
                $> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @repack.options -s 3gbm_HA.pdb -parser:protocol repack_HA.xml
                
                $> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @repack.options -s 3gbn_Ab_fixed.pdb -parser:protocol repack.xml 

        1. When the repacking runs are done (in about 3-4 minutes), copy the best scoring HA model to 3gbm_HA_repack.pdb and the best scoring antibody model to 3gbn_Ab_repack.pdb. (For brevity, we only generated a single structure. For actual production runs, we recommend generating a number of output structures, by adding something like "-nstruct 25" to the commandline. -- For the example outputs in the output_files/ directory, the lowest energy structures are 3gbm_HA_0011.pdb and 3gbn_Ab_fixed_0005.pdb.)
        
                 $> cp input_files/3gbm_HA_repacked.pdb .
                 $> cp input_files/3gbn_Ab_repacked.pdb .
    
        It can also be useful to pre-generate backbone conformational diversity prior to docking particularly when the partners are crystallized separately. The Rosetta FastRelax algorithm can be accessed through RosettaScripts. XML scripts, input files, options files and a command are available in the input_files directory. Backbone conformational diversity will not be explored in this tutorial due to time constraints.
        
    1. Orient the antibody in a proper starting conformation.
    
        Use available information on the participating interface residues to decrease the global conformational search space. This improves the efficiency of the docking process and the quality of the final model. In this benchmark case we will use the ideal starting conformation.
        
        1. Align the structures with pymol. (Note that you need to have your environment parameters set correctly in order to access pymol from terminal. You may need to open pymol manuall)
        
                > cp ../input_files/3gbm_native.pdb .
                pymol 3gbm_native.pdb 3gbm_HA_repacked.pdb 3gbn_Ab_repacked.pdb

            1. Type 'align 3gbn_Ab_repacked, 3gbm_native'
            1.  Type 'save 3gbm_HA_3gbn_Ab.pdb, 3gbm_HA_repacked + 3gbn_Ab_repacked'
        1. Renumber the pdb from 1 to the end without restarting. You can also copy a renumbered_version that is pre-generated for you:
        
                > $ROSETTA_TOOLS/protein_tools/scripts/pdb_renumber.py --norestart 3gbm_HA_3gbn_Ab.pdb 3gbm_HA_3gbn_Ab.pdb
                or
                $> cp input_files/3gbm_HA_3gbn_Ab.pdb .

1. Perform docking utilizing the RosettaScripts application.
    1. Prepare a RosettaScripts XML file for docking. This file outlines a protocol that performs docking. It then further minimizes the interface. Familiarize yourself with the contents of the script.
        1. Copy docking_full.xml from the input_files directory. Go through the xml script to understand the protocol.

                $> cp input_files/docking_full.xml .
                cat docking_full.xml

        1. Prepare an options file for docking.
            1. Copy the options file (docking.options) from the input_files directory. Familiarize yourself with the options in the file.

                    $> cp input_files/docking_full.options . 
                    $> cp input_files/docking_minimize.options .
                    cat docking.options

        1. Generate fifty models using the full docking algorithm. (The command below only generates one structure, as fifty will likely take a while - adjust the command and move on to the next steps while this is running.)

                $> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @docking_full.options -nstruct 1

        1. To have a good comparison for the native structure, we want to minimize the experimental structure into the Rosetta energy function. To do this, we run a similar protocol, but skipping the coarse and fine resolution search stages, keeping only the minimization stage. This provides us with a like-to-like comparison of the native structure.

            1. Copy docking_minimize.xml from the input_files directory.

                The docking_minimize.xml file differs from docking_full.xml only in the PROTOCOL section. The movers dock_low, srsc, and dock_high have been turned off by deleting the angle bracket at the beginning of these lines.

                    > cp input_files/docking_minimize.xml . 
                    > cat docking_minimize.xml

            1. Generate ten models using only the minimization refinement stage of docking.

                    > $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @input_files/docking_minimize.options -nstruct 10

1. Characterize the models and analyze the data for docking funnels.

    There are many movers and filters available in RosettaScripts for characterization of models. The InterfaceAnalyzerMover combines many of these movers and filters into a single mover. The RMSD filter is useful for benchmarking studies.

    The native structure used in this step (**3gbm_native.pdb**) has been cleaned as above. A complete structure is necessary for comparison to models. Missing density has been repaired through loop modeling or grafting of segments from 3gbn.pdb.
    
    **If your docking run is not finished yet, you can try out this step with pre-generated results. Make a new directory and copy the file docking.silent from the  output_files/ directory into this new directory.** You can also simply copy a smaller set of extracted files using:
    
            > cp output_files/*full*pdb output_files/*minimize*pdb .

    1.  Characterize your models using the InterfaceAnalyzer mover in RosettaScripts and calculate the RMSD to the native crystal structure with the RMSD filter.

            > cp input_files/docking_analysis.xml .
            > cp input_files/docking_analysis.options . 
            > cp input_files/3gbm_native.pdb .
            > $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @docking_analysis.options -in:file:s *full*pdb *minimize*pdb 

            sort -nk 7 docking_analysis.csv
            pymol 3gbm_native.pdb 3gbm_HA_3gbn_Ab_full_0001.pdb

    1. Plot various scores against rmsd for total_score, dG_separated, etc., to identify a binding funnel.
        1. Open docking_analysis.csv as a spreadsheet and create a scatter plot. (Check the boxes for space separation and merging delimiters.)

                ooffice docking_analysis.csv &

        1. Or use the provided R script to make score vs rmsd plots

                cp input_files/sc_vs_rmsd.R .
                Rscript ./sc_vs_rmsd.R docking_analysis.csv total_score 
                Rscript ./sc_vs_rmsd.R docking_analysis.csv dG_separated 
                Rscript ./sc_vs_rmsd.R docking_analysis.csv dG_separated.dSASAx100 
                gthumb *png &
 
