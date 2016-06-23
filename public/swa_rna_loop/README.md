# Application of the Stepwise Assembly method to RNA loop modeling (SWA_RNA_LOOP)

KEYWORDS: NUCLEIC_ACIDS LOOPS RNA

## Authors
- Parin Sripakdeevong (sripakpa [at] stanford.edu)
- Rhiju Das (rhiju [at] stanford.edu)

Last updated on March 18, 2012.
 
# Purpose and Algorithm

This demo illustrates a protocol to built single-stranded RNA loops using a deterministic, enumerative sampling method called Stepwise Assembly. The modeling situation considered here is the lock-and-key problem. Given a template PDB that contains nucleotides surrounding a missing RNA loop, the Stepwise Assembly method finds the loop conformation (the key) that best fits the surrounding structure (the lock). Details of this method is described in "An enumerative stepwise ansatz enables atomic-accuracy RNA loop modeling" by P. Sripakdeevong, W. Kladwang, and R. Das (2012), Proc Natl Acad Sci USA.
	As detailed in the paper, the Stepwise Assembly method constructs full-length RNA loops through the recursive building of each individual RNA nucleotides over multiple steps. The enumerative nature of the method makes the full-calculation quite computationally expensive, requiring for example 15,000 CPU hours to build a single 6-nucleotides RNA loop. While this full-calculation is now feasible on a high-performance computer clusters, performing the full-calculation in this demo would be too excessive. 
	Instead we will illustrate in this demo, the Stepwise Assembly protocol to build the first (5' most) nucleotide of a 6-nucleotides RNA loop. Performing the individual building step takes roughly 15 minutes on an Intel Core i7 2.66 GHz processor. The same Stepwise Assembly protocol can then be recursively applied to build the remaining 5 nucleotides of the loop (see SI of the referenced paper for details).

# Required Tools and Input files

There are two required files: 
The template_PDB file (template.pdb): A PDB file containing the coordinates of surrounding nucleotides in the vicinity of the missing RNA loop to be build. We recommend including all surrounding nucleotides within a 10-Angstrom vicinity of the missing RNA loop. Supplied PDB file must be in the Rosetta RNA PDB format.

The fasta file (fasta): this is the sequence file of the full-length RNA. The fasta file has the RNA name on the first line (after >), and the sequence on the second line. Valid letters are a, c, g and u. 

# Optional Tools and Input Files

## Optional additional files:
The native_PDB file (native.pdb): A PDB file containing the 'native' crystallographic or NMR structure. The PDB file should contain the coordinates of the native loop nucleotides plus the coordinates of the surrounding nucleotides inherited from template_PDB. The supplied native_PDB file is not used to guide the modeling process and only used for reporting the RMSD of the generated rosetta models to the native loop. Supplied PDB file must be in the Rosetta RNA PDB format.

# How to run the job

The SWA_RNA_python package located at `rosetta_tools/SWA_RNA_python/` contains the scripts necessary to setup and run the Stepwise Assembly protocol. Instructions are provided in steps 1)-4) below: 

1. Specify the location of the rosetta bin folder and rosetta database folder by editing the file `rosetta_tools/SWA_RNA_python/SWA_dagman_python/utility/USER_PATHS.py`

    For example, if the main rosetta folder is located at `~/rosetta/`, then the file should look as follow:
    ```
	    #!/usr/bin/python

	    USER_ROSETTA_BIN_FOLDER="~/rosetta/rosetta_source/bin/"
	    USER_ROSETTA_DATABASE_FOLDER="~/rosetta/rosetta_database/"
    ```

2. Add the SWA_RNA_python package location to the PYTHON path. For bash shell users, the location can be directly added to the `~/.bashrc` file:

    ```
	export PYTHONPATH=$PYTHONPATH:~/rosetta/rosetta_tools/SWA_RNA_python/
    ```

3. After the paths are correctly specified, the following command is used to setup everything needed run the Stepwise Assembly job:

    ```
	rosetta_tools/SWA_RNA_python/SWA_dagman_python/SWA_DAG/setup_SWA_RNA_dag_job_files.py -s template.pdb -fasta fasta -sample_res 3-8 -single_stranded_loop_mode True -local_demo True -native_pdb native.pdb
    ```

    The "-s" flag specifies the template_PDB file

    The "-fasta" flag specifies the fasta file

    The "-sample_res" flag specifies the sequence number of nucleotides in the missing loop. In this demo case, this correspond to nucleotides at sequence number 3 4 5 6 7 and 8.

    The "-single_stranded_loop_mode" flag specifies that the job involve modeling a single-stranded loop (i.e. the lock-and-key problem).

    The "-local_demo" flag indicate that this is demo to be run on a local laptop or desktop. The calculation perform here is to only build the first (5' most) nucleotide of the 6-nucleotides RNA loop.

    The "-native_pdb" flag specifies the native_PDB file and is optional. 

4. Type `source LOCAL_DEMO` to execute the Rosetta protocol.

The provided instruction will allow the user to build the first (5' most) nucleotide of a N-nucleotide loop. As previously stated, the full-calculation to build the full-length RNA loops is quite computationally expensive and is beyond the scope of this demo. The SWA_RNA_python package is, however, equipped to run this recursive full-calculation on a high-performance computer clusters. The package utilize concept familiar from the Map/Reduce Direct Acyclic Graph framework to order the calculation steps and allocate resources to recursive build the full-length RNA loop over multiple steps, one individual RNA nucleotide at a time. If any user is interested, please contact Parin Sripakdeevong (sripakpa [at] stanford.edu) and we will be happy to provide additional instructions.

# Expected Outputs

The expected outputs are two silent_files:
A) region_0_1_sample.out: This silent_file contain 108 structures, corresponding to the 108 lowest energy conformations.
B) region_0_1_sample.cluster.out: Same as A) but after clustering of the models to remove redundant conformations.

In both silent_files, the total energy score is found under the 'score' column. If the "native_pdb" flag was included, then the RMSD (in angstrom units) between the native_pdb and each Rosetta model is found under the 'NAT_rmsd' column.

Finally, use the following command to extract the top 5 energy cluster centers:

```
	rosetta_tools/SWA_RNA_python/SWA_dagman_python/misc/SWA_extract_pdb.py -tag S_0 S_1 S_2 S_3 S_4  -silent_file region_0_1_sample.cluster.out
```

After running the command, the extracted PDB files should appear in the pose_region_0_1_sample.cluster.out/ subfolder.

