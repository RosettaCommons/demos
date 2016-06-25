# Commonly used options (generic options)

KEYWORDS: CORE_CONCEPTS GENERAL

[[_TOC_]]

There are a hand full of options, that are useful most of the time you run Rosetta. Here, we will briefly go through those.

#### 1. Controlling input
###### Database (option normally not required) 
  
    -database /.../Rosetta/main/database

During a Rosetta run, a lot of information of pulled from the Rosetta database. It contains e.g. score files (weights for individual energy terms), bond lengths and atom connectivity for amino acids. **By default, Rosetta finds it's database automatically, and this option is NOT required**. If you have a non-standard intallation, specify the path to your local Rosetta database  Following the *-database* flag.
   
###### Provide input structures
Here is how you pass structure files to be used as input. This can be a structure where you want to re-design the sequence, re-design a loop, or two subunits that you want to dock together.   
   
    -in:file:s <pdb file>  
    -in:file:l <list of multiple pdbfiles>
    -in:file:native <pdb file>
    
 *-in:file:s* is the way to provide a single input structure. A *list* is a text file with pdb files (and their path, if located in a different directory). Rosetta will go through the files consecutively. The *native* flag gives you the opportunity to provide a reference structure for comparison. E.g if you have a guess how two proteins interact, you can compare to this structure or when re-designing the backbone of an existing protein.
    
     -ignore_unrecognized_res
     
This option is very important if Rosetta fails to read in your pdb file because of non-standard chemical entities in the file.

###### Provide input sequence
For *ab initio* structure predictions, you will have to provide a sequence in fasta format. Here is how you do that: 
 
    -in:file:fasta < e.g. my_sequence.fasta >
    
#### 2. Controlling output
###### Change default file names
Rosetta lets you to modify the names and format of output files as well as the amount of information given to you.  
Often you will see a *structure file* and a *score file* output. Structure files are named after the input file and automatically get numbers appended (e.g. MyFile_0001.pdb). You can change score file names and add a prefix or suffix to all output file names. All of these are optional.
 
    -out:prefix < e.g. trial2_ >
    -out:suffix < e.g. _trial2 >
    -out:file:scorefile < e.g. abinitio_scores.dat >
    -overwrite
    
The *-overwrite* option is needed if you already has a structure file (like MyFile_0001.pdb), and you re-run Rosetta with the same input file. In that case, Rosetta will fail if this option is missing.

###### Silent files - combine many output pdb files in a single compressed file
When running Rosetta, you will often get hundereds or thousands of individual pdb files. It can be useful to instead have a single file, where all structural information is combined in a binary file (mostly unreadable for humans). They take up less space, and let you extract only a small number of structures (e.g. best energy structures) after a run is finished. They get the *.out* file extension

    -out:file:silent_struct_type binary
    -out:file:silentfile < e.g. design_1.out >

Look at the **Analysis** chapter to find out how to extract pdb files from silent files.

###### Control the number of runs
    -nstruct < e.g. 1000 >

Typically, you run a simulation in Rosetta multiple times. This is important, because the sampling in an individual run is normally much too limited. *-nstruct 1000* will run 1000 iterations of a simulation and output 1000 files:
   
     MyFile_0001.pdb
     MyFile_0002.pdb  
     ...   
     MyFile_1000.pdb.  

    
    
