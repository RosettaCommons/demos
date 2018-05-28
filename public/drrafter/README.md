DRRAFTER: De novo RNP modeling in Real-space through Assembly of Fragments Together with Electron density in Rosetta 
=====================================================================================

KEYWORDS: NUCLEIC_ACIDS EXPERIMENTAL_DATA RNA DENOVO STRUCTURE_PREDICTION  

Written in March 2018 by Kalli Kappel (kappel at stanford dot edu). Updated May 2018.    

This demo shows how to use DRRAFTER to build a structure of an RNA-protein complex into a cryoEM density map and how to estimate the error in the resulting models.  

## Installing DRRAFTER:  
1. Download Rosetta [here](https://www.rosettacommons.org/software/license-and-download). You will need to get a license before downloading Rosetta (free for academic users). DRRAFTER is available in the Rosetta weekly releases starting with 2018.12. **DRRAFTER is NOT available in Rosetta 3.9.** 
2. Install Rosetta following the instructions available [here](https://www.rosettacommons.org/docs/latest/build_documentation/Build-Documentation).  
3. Make sure that you have python (v2.7) installed.
4. Install Rosetta RNA tools. See instructions and documentation [here](https://www.rosettacommons.org/docs/latest/application_documentation/rna/RNA-tools).
5. Check that the ROSETTA environmental variable is set (you should have set this up during RNA tools installation). Type `echo $ROSETTA`. This should return the path to your Rosetta directory. If it does not return anything, go back to step 4 and make sure that you follow the steps for RNA tools setup.  
6. Set up the executables for DRRAFTER. Type:
```
ln -s $(ls $ROSETTA/main/source/bin/rna_denovo* | head -1 ) $ROSETTA/main/source/bin/rna_denovo
```
Then type: 
```
ln -s $(ls $ROSETTA/main/source/bin/drrafter_error_estimation* | head -1 ) $ROSETTA/main/source/bin/drrafter_error_estimation
``` 
7. Add the path to the DRRAFTER script to your $PATH (alternatively, you can type the full path to the DRRAFTER.py script each time that you use it). It is found in `main/source/src/apps/public/DRRAFTER/` in your Rosetta directory. An example for bash:
```
export PATH=$PATH:$ROSETTA/main/source/src/apps/public/DRRAFTER/
```

## Brief explanation of input files:  

All of the necessary files for this demo are available in `$ROSETTA/demos/public/drrafter/`, where `$ROSETTA` is the path to your Rosetta installation.   

`fasta.txt`: The FASTA file listing the full sequence of the complex being modeled. It should contain at least one line that starts with '>' and lists chains and residue numbers for the sequence, e.g. A:136-258 E:1-23. Here we are modeling chain A residues 136-258 and chain E residues 1-23. The subsequent lines should list the full sequence of the complex. Protein residues are specified by uppercase one-letter codes. RNA residues are specified with lowercase one-letter codes ('a', 'u', 'g', and 'c').  

`secstruct.txt`: A file containing the secondary structure of the complex in dot-bracket notation. Secondary structure for the protein should be specified by dots. The secondary structure should be the same length as the sequence found in the fasta file. For RNA residues, this secondary structure will be enforced during the DRRAFTER run.  

`1wsu_simulated_7A.mrc`: The density map file in mrc format (ccp4 is also acceptable). For this demo, this map has been simulated from PDB ID 1WSU at 7Å.   

`RNA_helix.pdb`: This is an ideal RNA helix corresponding to residues E:1-4 and E:20-23. Ideal RNA helices can be generated using `rna_helix.py` (documentation on the [RNA tools page](https://www.rosettacommons.org/docs/latest/application_documentation/rna/RNA-tools)).  

`protein_fit_into_density.pdb`: The protein structure that has been fit into the density map. For this demo, this is a crystal structure of the unbound protein (PDB ID 1LVA). The sequence of this protein structure should exactly match the sequence in the fasta file.  

`protein_and_RNA_helix_fit_into_density.pdb`: Coordinates of both the protein structure and the RNA helix (from `RNA_helix.pdb`) fit into the density map. This will be the starting conformation for the DRRAFTER run. (The relative positions of the protein and RNA helix will be allowed to change though.)  


## Running DRRAFTER:

1. Use `DRRAFTER.py` to set up the run. Type:  
```
DRRAFTER.py -fasta fasta.txt -secstruct secstruct.txt -start_struct protein_and_RNA_helix_fit_into_density.pdb -map_file 1wsu_simulated_7A.mrc -map_reso 7.0 -residues_to_model E:1-23 -include_as_rigid_body_structures protein_fit_into_density.pdb RNA_helix.pdb -absolute_coordinates_rigid_body_structure protein_fit_into_density.pdb -job_name demo_run -dock_into_density -demo_settings -rosetta_directory /your/path/to/rosetta/executables
```

**Note** that `/your/path/to/rosetta/executables` needs to be replaced with the actual path to your Rosetta executables. Or if the path to the Rosetta executables is already in your system PATH, then you can omit the `-rosetta_directory` flag.   
**Note** also the `-demo_settings` flag: this flag is designed to make the DRRAFTER run finish quickly, and should not be used for normal runs.  

This will create the following files:  

`fasta_demo_run.txt`: The FASTA file for the region that will be included in the DRRAFTER run.   

`secstruct_demo_run.txt`: The file specifying the secondary structure for the region that will be included in the DRRAFTER run.  

`coord_csts_demo_run.txt`: A file describing coordinate restraints that will be applied during the DRRAFTER run. By default, residues in the starting structure will be restrained to be within 10 Å of their initial coordinates. This can be turned off with -no_csts and the distance at which the restraints will be activated (default 10 Å) can be controlled with -cst_dist.   

`flags_demo_run`: A file listing all of the Rosetta options for the run.  

`init_struct_demo_run.pdb`: The starting structure for the DRRAFTER run. In this case, this structure is identical to `protein_and_RNA_helix_fit_into_density.pdb`.  

`DRRAFTER_command`: This file contains the command to run the DRRAFTER job.  

2. Run the `DRRAFTER_command`. This can be done either by typing:  

```
source ./DRRAFTER_command
```

OR by copying the line in the `DRRAFTER_command` file to the command line:

```
/your/path/to/rosetta/executables/rna_denovo @flags_demo_run
```

This will take several minutes to run and it should create a file called `demo_run.out`, which contains all of the structures from the run.

3. Extract PDB files from the compressed output file created in the previous step. Type:  

```
extract_lowscore_decoys.py demo_run.out 10
```

This extracts the 10 best scoring structures from the run and will create 10 PDB files named `demo_run.out.1.pdb`, `demo_run.out.2.pdb`, etc. Note that in this case we only generated 10 structures total, but for a real run it is recommended that you build at least 2000-3000 structures.  

4. Look at the structures! Open them in pymol or Chimera. You will see that all of the missing RNA residues have been built. The protein and RNA helix have also moved slightly from their initial positions.  

5. Estimate the error in the DRRAFTER models. Again, use the `DRRAFTER.py` script to do this. Type:  

```
DRRAFTER.py -final_structures demo_run.out.1.pdb demo_run.out.2.pdb demo_run.out.3.pdb demo_run.out.4.pdb demo_run.out.5.pdb demo_run.out.6.pdb demo_run.out.7.pdb demo_run.out.8.pdb demo_run.out.9.pdb demo_run.out.10.pdb -estimate_error -rosetta_directory /your/path/to/rosetta/executables
```

This will print out information about the error estimation to your screen. For example:  

```
apps.public.DRRAFTER.drrafter_error_estimation: #############################################
apps.public.DRRAFTER.drrafter_error_estimation: Mean pairwise RMSD (convergence): 5.41042
apps.public.DRRAFTER.drrafter_error_estimation: Estimated minimum RMSD: 3.65777
apps.public.DRRAFTER.drrafter_error_estimation: Estimated mean RMSD: 5.39297
apps.public.DRRAFTER.drrafter_error_estimation: Estimated RMSD of median structure: 4.73551
apps.public.DRRAFTER.drrafter_error_estimation: Median structure: demo_run.out.1.pdb
apps.public.DRRAFTER.drrafter_error_estimation: #############################################
```

All numbers have units of Å. The mean pairwise RMSD describes the "convergence" of the run, i.e. how similar the final structures are to each other. The estimated RMSD (root mean square deviation) values to the "true" coordinates are based on this convergence value. The estimated minimum RMSD predicts the best accuracy of the final structures. The estimated mean RMSD predicts the average RMSD accuracy of the final structures. The median structure is determined to be the final structure with the lowest average pairwise RMSD to the other final structures. The accuracy estimate of this model is also printed to the screen.  

For reference, example output is provided in the `example_output/` directory.   


## Additional information

See the DRRAFTER documentation [here](https://www.rosettacommons.org/docs/latest/application_documentation/rna/drrafter).
