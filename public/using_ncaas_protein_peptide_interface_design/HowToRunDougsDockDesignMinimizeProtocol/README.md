## Part IV    INSTRUCTIONS FOR RUNNING THE DOUGSDOCKDESIGNMINIMIZE PROTOCOL


KEYWORDS: DOCKING INTERFACES STRUCTURE_PREDICTION

The DougsDockDesignMinimize (DDDM) protocol was used in the accompanying manuscript to redesign the protein/peptide interface of Calpain and a fragment of its inhibitory peptide calpastatin. The protocol was written for this specific protein/peptide interaction and modifications to the code will be necessary to run the protocol on a different system. A modified form was used as the example protocol in the advanced section of the Rosetta 3.0 release manual.

Minor modifications to the protocol have been made from the version used to produce the designs in the accompanying publication. The interfaces of the Rosetta libraries have changed since the initial implementation of the protocol and these modifications were necessary to allow the protocol to work with the current release. 

The protocol has two loops referred to as the inner and outer loops. The outer loop controls the number of structures generated, outputting structures only if certain filters are met. The inner loop iterates between two phases: a perturbation phase and a design phase. During the perturbation phase three types of perturbations are used. Two types of perturbations are applied to the peptide: rotational and translational rigid body permutation of the peptide in the binding pocket of the protein, and small and shear perturbations of the backbone dihedral angles in all residues of the peptide. The third type of perturbation is preformed on residues 1-45 of the N-terminus which comprise the peptide binding site and surrounding residues of the protein are perturbed using small and shear perturbations of the backbone dihedral angles. One of the three perturbations is randomly chosen to be used each time the function in called and is followed by a coarse side chain optimization (RotamerTrials). During the design phase of the protocol the side chains positions are optimized using a more intensive side chain optimization routine (PackRotamers) followed by minimization of the backbone and side chain dihedrals angles as well as the jump between the peptide and protein. The number of perturbations per design as well as the magnitude of perturbations are controlled using command line options described bellow.

There are three basic steps to running the protocol as it was used in the accompanying manuscript.
 - generating the input resfiles and folders
 - running the DougsDocDesignMinimize protocol for each mutation at each position 
 - running the analysis scripts

Example output is provided for a small version of the full run.

------------------------------------------
   GENERATING INPUT RESFILES AND FOLDERS
------------------------------------------
The Protein Databank code for the Calpain/Calpastatin structure used in the designs is 1NX1. The files contain two copies of calpain (chains A and B) and two copies of the calpastatin inhibitory peptide (chains C and D). To reduce the size of the system, designs were done on chain A and chain C only because chains A and C have lower B-factors than chains B and D. The calpain/calpastatin interface is distal to the calpain/calpain interface and the resides that make up the calpain/calpain interface were held fix during the protocol. Additionally the calcium atoms and water molecules were removed. Rosetta doesn't work with with water molecules and the calcium ions are not near the protein peptide interface. This input pdb was repacked and both the side chain and backbone dihedrals minimized using a modified version of the fixbb app. The input PDB file is called 1NX1_clean_repack_min_all.pdb

Each of the NCAAs added in the accompanying publication was tried at each position in the peptide. To do this and to keep all of the output pdb organized a script is provided that creates folders and generates resfiles based on templates for each sequence. To generate the resfiles and folders run the following commands...

```
$> ln -s HowToRunDougsDockDesignMinimizeProtocol/scripts .
$> ln -s HowToRunDougsDockDesignMinimizeProtocol/run_dir/* .
$> scripts/make_folders_resfiles.tcsh
```

NOTE: To save space, the script has been modified to only produce resfiles and folders for position 610 in the peptide and residue type MPA (4-methyl-phenylalanine). To modify the script to produce folders and resfiles for each position simply uncomment the lines that read "#foreach i ( 601 602 603 604 605 606 607 608 609 610 611 )" and "#foreach j ( ABA APA HLU ..... C92 C93 C94 )" comment out the line that reads "foreach i ( 610 )" and "foreach j ( MPA )".

Using the templates in the run_dir the make_folders_resfiles.tcsh script makes folders and resfiles to preform all of the DDDMI runs.



###   RUNNING DOUGSDOCKDESIGNMINIMIZEINTERFACE PROTOCOL


To run the protocol modifications need to made to files in the database. In the file rosetta_database/chemical/residue_type_sets/fa_standard/residue_types.txt all of the paths to the residue type parameter files under the L-NCAA heading need to be uncommented by removing the "#" from the front of the line. Additionally the rotamer libraries for the NCAA are not provided in the default Rosetta database because they are more than 400MB. The rotamer libraries for the NCAAs added in the accompanying publication are provided as supplemental information. 

NOTE: Turning on all of the additional residue types dramatically increases the number of residue types and the memory footprint of Rosetta. The memory foot print can be reduced by commenting out unnecessary patches in the rosetta_database/chemical/residue_type_sets/fa_standard/patches.txt file. For the DougsDockDesignMinimizeProtocol all but the NtermProteinFull.txt and CtermProteinFull.txt can be safely commented out by placing "#" symbols at the beginning of each line of the patches.txt except for the lines that say "NtermProteinFull.txt" and "CtermProteinFull.txt".

To run the protocol as in the accompanying publication preform the following commands starting at the HowToRunDougsDockDesignMinimizeProtocol directory.

```
$> ln -s pos_610_MPA/* .
$> $ROSETTA/bin/doug_dock_design_min_mod2_cal_cal.linuxgccrelease -database /PATH/TO/rosetta_database -s inputs/1NX1_clean_repack_min_all.pdb -resfile resfile_pos_603_MPA -nstruct 1 -inner_num 45 -pert_num 25 -ia_ener 100 -use_input_sc -pdb_gz
```

>**change nstruct to 255!**

The above command generates 255 structures and will take approximately 5 minutes per structure depending on your hardware.

NOTE: The extension of your executable maybe different than the above. Also in the publication the nstruct command line option was 255. However to save space for the protocol capture an nstruct of 10 was used.

A script is provided that will preform the above command for each folder created by the make_folders_resfiles.tcsh script.

```
$> for i in pos_*; do cd $i $ROSETTA3/doug_dock_design_min_mod2_cal_cal.linuxgccrelease ../inputs/ \  
1NX1_clean_repack_min_all.pdb -resfile ../resfile_$i \  
-nstruct 1 -inner_num 45 -pert_num 25 -ia_ener 100 \  
-use_input_sc -pdb_gz >& $i.log; \ 
 echo "Finished $i"; cd ../; done
          
```

NOTE: You will need to set the path to your database and executable in the run_script.bash.

There are a number command line options to control the running of the application.

pert_mc_temp: The MC temperature to use for the perturbation phase of the DDDM protocol, defaults to 0.8 kT.
pert_dock_rot_mag: The rotation magnitude for the ridged body perturbation in the perturbation phase of the DDDM protocol, defaults to 0.5. 
pert_dock_trans_mag: The translation magnitude for the ridged body perturbation in the perturbation phase of the DDDM protocol, defaults to 0.25. 
pert_pep_small_temp: The temperature for the internal MC object in the small mover of the peptide perturbations, defaults to 0.8 kT. 
pert_pep_shear_temp: The temperature for the internal MC object in the shear mover of the peptide perturbations, defaults to 0.8 kT. 
pert_ter_small_temp: The temperature for the internal MC object in the small mover of the termini perturbations, defaults to 0.8 kT. 
pert_ter_shear_temp: The temperature for the internal MC object in the shear mover of the termini perturbations, defaults to 0.8 kT. 
pert_pep_small_H: The maximum angle of perturbation for helical secondary structure for the peptide small mover, defaults to 1 degree.
pert_pep_small_L: The maximum angle of perturbation for loop secondary structure for the peptide small mover, defaults to 1 degree.
pert_pep_small_E: The maximum angle of perturbation for strand secondary structure for the peptide small mover, defaults to 1 degree.
pert_pep_shear_H: The maximum angle of perturbation for helical secondary structure for the peptide shear mover, defaults to 1 degree
pert_pep_shear_L: The maximum angle of perturbation for loop secondary structure for the peptide shear mover, defaults to 1 degree.  
pert_pep_shear_E: The maximum angle of perturbation for strand secondary structure for the peptide shear mover, defaults to 1 degree.
pert_ter_small_H: The maximum angle of perturbation for helical secondary structure for the termini small mover, defaults to 1 degree
pert_ter_small_L: The maximum angle of perturbation for loop secondary structure for the termini small mover, defaults to 1 degree.  
pert_ter_small_E: The maximum angle of perturbation for strand secondary structure for the termini small mover, defaults to 1 degree.
pert_ter_shear_H: The maximum angle of perturbation for helical secondary structure for the termini shear mover, defaults to 1 degree
pert_ter_shear_L: The maximum angle of perturbation for loop secondary structure for the termini shear mover, defaults to 1 degree.  
pert_ter_shear_E: The maximum angle of perturbation for strand secondary structure for the termini shear mover, defaults to 1 degree.
pert_pep_num_rep: Number of small and shear iterations for the peptide, defaults to 100. 
pert_ter_num_rep: Number of small and shear iterations for the terminus, defaults to 100. 
pert_num: Number of iterations of perturbation loop per design, defaults to 100.
inner_num: Number of iterations of the inner loop, defaults to 100. 
ia_ener: Upper energy limit for final design/interface analysis checkpoint, defaults to 0.0. 
desn_mc_temp: The temperature to use for the design/minimization phase of the DDDM protocol, defaults to 0.8 kT.

For the most part the defaults should suffice and are what was used in the paper. The "ia_ener" command line option is dependent on the complex that the protocol is run on. Setting it unreasonably low will cause the protocol to run forever and never output any structures. Setting it to 100 above allows all but the most aggressions structures to be output. The nstruct command line option controls the outer loop. 

-----------------------------
   ANALYZING THE RESULTS
-----------------------------

Each of the output pdb files contains information about the protein peptide complex that can used to evaluate the designs.  For example at the end of the 1NX1_clean_repack_min_all_0004.pdb is shown bellow. For each filter the value is calculated for the protein and peptide together (COMPLEX) and separated by 1000 angstroms (SEPARATE) and the difference between the two (DIFF). ENERGY is the Rosetta energy. The ENERGY_COMPLEX is the primary determinant to how good a design is and the ENERGY_DIFF can give an estimate for the binding energy. SASA is the solvent accessible surface area. SASA_DIFF is indicative of sequences that make a more protein-peptide contacts and can be used for screening designs for example placing a very large side chain at a constrained interface position can cause the peptide to be pushed out of the binding pocket which would be reflected in a smaller magnitude SASA_DIFF. HB_ENER is the hydrogen bonding component of the Rosetta energy. Larger HB_ENER_DIFF values indicate that the design is making more or better hydrogen bonds across the protein peptide interface. PACK is the RosettaHoles score and is a measurement of how well the protein is packed. A PACK_COMPLEX that is larger than the PACK_SEPARATE is favorable and suggests that the complex is better packed than the protein alone. Additionally the RosettaHoles score penalizes holes that cannot be occupied by solvent so larger PACK_DIFF score indicate that the designed peptide is capable of filling cavities in the protein that are inaccessible to solvent.

```
ENERGY_COMPLEX:	   -63.3501
ENERGY_SEPERATE:   -48.599
ENERGY_DIFF:	   -14.7512
SASA_COMPLEX:	   11774.1
SASA_SEPERATE:	   13092.3
SASA_DIFF:	   -1318.15
HB_ENER_COMPLEX:   -146.728
HB_ENER_SEPERATE:  -145.037
HB_ENER_DIFF:	   -1.69085
PACK_COMPLEX:	   0.515277
PACK_SEPERATE:	   0.466995
PACK_DIFF:	   0.0482824
```

A script is provided that pulls the information out of a set of pdb files and sorts it based on the ENERGY_COMPLEX metric. 

```
$ cd pos_610_MPA
$ ../scripts/get_interface_data.tcsh
```

