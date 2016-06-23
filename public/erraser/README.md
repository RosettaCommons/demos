ERRASER: Enumerative Real-Space Refinement ASsitted by Electron density under Rosetta
=====================================================================================

KEYWORDS: NUCLEIC_ACIDS EXPERIMENTAL_DATA RNA

This demo illustrates the ERRASER (Enumerative Real-Space Refinement ASsitted 
by Electron density under Rosetta) protocol, which improves an RNA 
crystallographic model using Rosetta under the constraint of experimental 
electron density map. It was written in Mar. 2012, by Fang-Chieh Chou (fcchou 
at stanford dot edu) and based on a paper to be published:

* Chou, F.C., Sripakdeevong, P., Dibrov, S.M., Hermann, T., and Das, R. Correcting pervasive errors in RNA crystallography with Rosetta, arXiv:1110.0276. 

A preprint is available at:

* http://arxiv.org/abs/1110.0276

Setting up the demo
-------------------

The example input files are in rosetta_input; you may wish to copy them locally 
with the command:

    cp rosetta_inputs/* ./

Python codes needed to run the job are located at 
rosetta/rosetta_tools/ERRASER/

The following setup steps are required prior to running this demo.

1. Download and install PHENIX from http://www.phenix-online.org/. PHENIX is 
   free for academic users.

2. Ensure you have correctly setup PHENIX. As a check, run the following 
   command:
   
       phenix.rna_validate 

3. Check if you have the latest python (v2.7) installed. If not, go to the 
   rosetta/rosetta_tools/ERRASER/ folder and run 

        ./convert_to_phenix.python

    This will change the default python used by the code to phenix-built-in 
    python, instead of using system python.

4. Set up the environmental variable "$ROSETTA", point it to the Rosetta 
   folder. If you use bash, append the following lines to ~/.bashrc:

        ROSETTA=YOUR_ROSETTA_PATH; export ROSETTA" # Change YOUR_ROSETTA_PATH to the path in your machine!

    Also add the ERRASER script folder to $PATH. Here is a bash example:

        PATH=$PATH:YOUR_ROSETTA_PATH/rosetta_tools/ERRASER/" # Change YOUR_ROSETTA_PATH to the path in your machine!

Now you are ready to go!

Running the demo
----------------

As a fast quick run, run the following command:

    erraser.py -pdb 1U8D_cut.pdb -map 1U8D_cell.ccp4 -map_reso 1.95 -fixed_res A33-37 A61 A65

In this example, we specify the input pdb file, ccp4 map file and the map 
resolution. The -pdb and -map are required option, and -map_reso is optional but 
recommended (default is 2.0 if no input is given). In this example, the 
"1U8D_cut.pdb" file is a segment cutting from a deposited PDB file. The input 
pdb file should follow the standard PDB format, and no pre-processing is 
needed. The input map must be a CCP4 2mFo-DFc map. To avoid overfitting, Rfree 
reflection should be removed during the creation of the map file.

Note that we also manually fixed the position of residue 33-37, 61 and 65 in 
chain A, therefore we will only optimize residue 62-64. The -fixed_res argument 
is optional.

After the job finished successfully, you should see the output pdb file 
"1U8D_cut_erraser.pdb" in the current folder. A sample output is in the 
example_output folder for comparsion.

Note that the output file is in the standard PDB format and inherits all the 
ligands, metals and waters from the input pdb file (these atoms are not 
optimized in ERRASER). The user can then refine the output model using PHENIX 
or other refinement packages without any post-processing.

By inspecting the structure, you should be able to see a backbone conformation 
change at residue 63-64 after ERRASER.

Arguments for erraser.py (for your reference)
---------------------------------------------

#### Required:

* -pdb  
  Format: -pdb \<input pdb>  

  The starting structure in standard pdb format

* -map  
  Format: -map <map file>  

  2mFo-DFc map file in CCP4 format. Rfree should be excluded.

#### Commonly used:

* -map_reso  
  Format: -map_reso <float>  
  Default: 2.0  

  The resolution of the input density map. It is highly recommended to input 
  the map resolution whenever possible for better result.

* -out_pdb  
  Format: -out_pdb <string>  
  Default: \<input pdb name>\_erraser.pdb.  

  The user can output to other name using this option.

* -n_iterate  
  Format: -n_iterate <int>  
  Default: 1  

  The number of rebuild-minimization iteration in ERRASER. The user can 
  increase the number to achieve best performance. Usually 2-3 rounds will be 
  enough. Alternatively, the user can also take a ERRASER-refined model as the 
  input for a next ERRASER run to achieve mannual iteration.

* -fixed_res  
  Format: -fixed_res <list>  
  Default: <empty>  
  Example: A1 A14-19 B9 B10-13  (chain ID followed by residue numbers)  

  This allows users ton fix selected RNA residues during ERRASER. For example, 
  because protein and ligands are not modeled in ERRASER, we recommand to fix 
  RNA residues that interacts strongly with these unmodeled atoms. ERRASER will 
  automatically detect residues covalently bonded to removed atoms and hold 
  them fixed during the rebuild, but users need to specify residues having 
  non-covalent interaction with removed atoms mannually.

* -kept_temp_folder  
  Format: -kept_temp_folder <True/False>  
  Default: False  

  Enable this option allows user to examine intermediate output files storing 
  in the temp folder. The default is to remove the temp folder after job 
  completion.

#### Other:

* -rebuild_extra_res  
  Format/Default: Same as -fixed_res  

  This allows users to specify extra residues and force ERRASER to rebuild 
  them. ERRASER will automatically pick out incorrect residues, but the user 
  may be able to find some particular residues that was not fixed after one 
  ERRASER run. The user can then re-run ERRASER with -rebuild_extra_res 
  argument, and force ERRASER to remodel these residues.

* -cutpoint_open  
  Format/Default: Same as -fixed_res  

  This allows users to specify cutpoints (where the nucleotide next to it is 
  not connected to itself) in the starting model. Since ERRASER will detect 
  cutpoints in the model automatically, the users usually do not need to 
  specify this option.

* -use_existing_temp_folder  
  Format: -use_existing_temp_folder <True/False>  
  Default: True  

  When is True, ERRASER will use any previous data stored in the existing temp 
  folder and skip steps that has been done. Useful when the job stopped 
  abnormally and the user try to re-run the same job. Disable it for a fresh 
  run without using previously computed data.

* -rebuild_all  
  Format: -rebuild_all <True/False>  
  Default: False  

  When is True, ERRASER will rebuild all the residues instead of just 
  rebuilding errorenous ones. Residues in "-fixed_res" (see below) are still 
  kept fixed during rebuilding. It is more time consuming but not necessary 
  leads to better result. Standard rebuilding with more iteration cycles is 
  usually prefered.

* -native_screen_RMSD  
  Format: -native_screen_RMSD <float>  
  Default: 2.0  

  In ERRASER default rebuilding, we only samples conformations that are within 
  2.0 A to the starting model (which is the "native" here). The user can modify 
  the RMSD cutoff. If the value of native_screen_RMSD is larger than 10.0, the 
  RMSD screening will be turned off.
