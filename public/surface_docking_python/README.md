# PyRosetta RosettaSurface

KEYWORDS: CORE_CONCEPTS GENERAL

## Authors
Written by Emily Koo and Michael Pacella (Graylab), mpacella88@gmail.com

## General Description
This demo will describe how to run the PyRosetta version of the RosettaSurface 
algorithm.  The test case shown here will analyze the adsorption of the 
biomineralization protein osteocalcin to the 100 surface of the mineral
hydroxyapatite

## Algorithm
Simultaneous optimization of protein rigid-body orientation, backbone and 
side chain conformations on a solid surface.  

## Commands
1. `cd rosetta_inputs/Osteocalcin demo`
2. `surface_docking.py@flags`

## Input Files
- 1Q8H100.pdb = input pdb file, extended osteocalcin positioned above HAp 100
- 1Q8H_disulf = disulfide bond specification file
- 1Q8H_native.pdb = native crystal structure of osteocalcin
- aa1Q8H_03_05.200_v1_3 = 3mer fragments for osteocalcin
- aa1Q8H_09_05.200_v1_3 = 9mer fragments for osteocalcin
- flags = arguments for rosetta surface


## PyRosetta RosettaSurface Setup and Pre-processing

Before running the scripts for the first time, make sure that all the following details are addressed.
 
1. When PyRosetta is installed, the install directory may be different from user to user, so adding 
the scripts directories to the PATH environment variable in .bashrc file will allow the scripts to be 
run from any directory, independent from the actual location of the scripts. The following 
instructions show how the directories are added:

    1. Create/open .bashrc in home directory
        ```
        vi ~/.bashrc
        ```
        
    2. If file exists, skip to 3. Else, add the following lines to the file:
        ```
        # .bashrc

        # Source global definitions
        if [ -f /etc/bashrc ]; then
            . /etc/bashrc
        fi

        export PATH=$PATH
        ```
        
    3. Add all directories containing scripts to end of PATH statement, delimited by a colon
        ``` 
        export PATH=$PATH:/path/to/scripts:/path/to/scripts/pdb_objects/
        ``` 
        where /path/to/scripts and /path/to/more/scripts should be modified to the correct directories.

    4. Make sure that the SetPyRosettaEnvironment.sh is sourced every session. Add 
        ```
        source /path/to/PyRosetta/SetPyRosettaEnvironment.sh
        ```
                
    5. Save and close file, then source it.
        ```
        :wq
        source ~/.bashrc
        ```

2. Make sure the python and bash scripts are given the permission to be executed. Run the following 
command in the directory containing the scripts
    ```
        chmod +x *.py *.sh ./pdb_objects/*
    ```
        
3. Gnuplot is the plotting program used to generate the plots, so the program path has to be modified 
in each Plot\* file for the plots to be generated automatically. The version used is 4.2.

4. To generate plots, simply run:
    ```
        PostProcessRS.sh
    ``` 
    Make sure that the folder contains either only adsorbed state PDBs (and native.pdb) or solution state PDBs. 
    
## Post-processing
(in the directory with output decoys)

1. GetTop.sh Ads 4
2. cd TOP4.Ads
3. PostProcessRS.sh

These commands will extract the top 4 adsorbed-state decoys for analysis and generate
secondary structure, protein-protein contact maps, and protein-surface contact maps

## Example output:
- Ads\_SecStruct.png = adsorbed state secondary structure
- Sol\_SecStruct.png = solution state secondary struccture
- Ads\_ContactMap = adsorbed state contact map
- Sol\_ContactMap = solution state contact map
- Surface\_ContactMap.png = protein-surface contact map
- Diff\_ContactMap.png = difference map between contacts in the ads/sol state

## Limitations
This app requires a single protein positioned above a solid surface whose parameters are 
present in Rosetta.  The protein needs to appear before the surface in the pdb file
and the two need to be separate chains

