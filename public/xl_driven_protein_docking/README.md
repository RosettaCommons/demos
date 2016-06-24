#Crosslink Guided Protein-Protein Docking

KEYWORDS: DOCKING EXPERIMENTAL_DATA STRUCTURE_PREDICTION

Original author Abdullah Kahraman (abdullah.kahraman@gmail.com).  This file was updated on 22 June 2016 by Vikram K. Mulligan, Ph.D. (vmullig@uw.edu), as part of the 2016 Documentation eXtreme Rosetta Workshop (XRW).

## Contents
-----------
- Introduction
- Demo Run
- System Requirements
- Required Software
- Input Files
- Commandline options for the excutable
    - docking_protocol
    - Xwalk
    - NACCESS
    - Superimpose
    - AvgInterface

## Introduction
---------------
Distance restraints from chemical cross-linking experiments combined with mass-spectrometry can
guide protein-protein docking calculations and significantly improve the accuracy of the
simulations. This protocol uses cross-link information for the prediction of protein complexes using
the standard ROSETTA docking protocol. Central to the protocol is the application of Xwalk, a
software able to simulate cross-links on protein surfaces, NACCESS for calculating the size of
predicted binding interfaces and a quality threshold clustering approach. The protocol requires
the structural coordinates of two proteins and distance restraints from chemical cross-link
experiments and predicts the structure of binary protein complexes and their binding interface.

The cross-link guided protein docking workflow consists of following steps:

1. Relax the receptor and ligand PDB structures (using ROSETTA's relax protocol) or prepack their side chains (using ROSETTA's docking_prepack_protocol).

2. Perform global docking in low resolution centroid mode producing ~100,000 models (using ROSETTA's docking_protocol), while applying the cross-link data as distance restraints in the scoring function.

3. Post-dock filter all predicted models (using Xwalk) in order to determine which models satisfy most cross-links.

4. Select the 500 models with the lowest ROSETTA energy scores and compute binding interface sizes (using the NACCESS software) to select only for those models that have a sufficiently large binding interface.

5. Perform a Quality-Threshold (QT) clustering (using the Superimpose application) on all models having passed the interface size filter and choose the lowest scoring model from the largest three clusters.

6. Perform a local refinement docking calculation on each of the three cluster representative ~3 x 5000 models (using ROSETTA's docking_protocol), while applying the cross-link data as distance constraints in the scoring function.

7. Repeat step 3.

8. Repeat step 4.

9. Perform an all vs all RMSD calculation combined with hierarchical clustering to select the lowest scoring model from the three largest clusters as the best prediction from the entire docking run.

10. Calculate the contact frequency for each amino acid with all local refined models that satisfy to most cross-links and have a sufficiently large interface size in order to predict the interface between the receptor and ligand structure.

## Demo Run
--------
This directory has following file structure:
-   example_output -> A sample of a finished cross-link guided docking workflow run.
-   inputs         -> Input files required for performing the cross-link guided docking workflow.
-   output         -> The directory in which a new cross-link guided docking workflow will be
                      executed.
-   README.txt     -> This read me file.
-   run_demo.sh    -> A master-SHELL-script executing automatically each of the scripts in the
                      scripts directory. Please type ./run_demo.sh for more help.
-   scripts        -> A set of SHELL and PERL scripts performing each of the steps above.

The set of scripts performing each of the cross-link guided docking steps above are provided in the
scripts directory and can be sequentially executed with the ./run_demo.sh master script. To run the
master script, you will need to provide the various paths different software tool. The list of the
required software tools together with their download link can be found in the table below. Please
type ./run_demo.sh for a list of the application and their order at the commandline.

As demo scripts are intended to only provide a quick example of the individual steps in the
cross-link guided docking workflow, by default
   - only prepacking of the relax and ligand structures are performed at step 1.
   - only up to 5 models in the global docking run are produced at step 2.
   - only up to 5 models per cluster representative are produced in the local docking refinement at
     step 6.
   - models having a too small interface are not removed at step 4 and 8.

These parameters can be changed within the run_demo.sh script.
The total run time of the cross-link guided docking demo is around 10 minutes on a modern computer.

The result will be a bestModel subdirectory within the output directory. The subdirectory should
hold up to three best predictions from the cross-link guided docking workflow and a PDB file with
the predicted interface for that particular protein complex. The latter PDB hold the absolute and
relative contact frequencies in the occupancy and temperature factor column, respectively. It can
be loaded into PyMOL and colored e.g. with following commands:
   - absolute contact frequency: cmd.spectrum("q","green_white_red","predicted_interface",-1,1,0,0)
   - relative contact frequency: cmd.spectrum("b","green_white_red","predicted_interface",-1,1,0,0)

## System Requirements
----------------------
- Unix like operating system, recommended and tested are Mac OS X 10.8.2 and Linux CentOS 6.
- bash or csh shell
- JAVA version 1.6
- PERL version 5

## Required Software
--------------------
Name           |   Description                                                          |   Version   |   URL
------------------------------------------------------------------------------------------------------------------------------------------
ROSETTA        |   Molecular modeling suite. Protocols for structure relax, side-chain  |   3.4       |   http://www.rosettacommons.org
               |   prepacking and protein-docking are called relax,                     |             |
               |   docking_prepack_protocol and docking_protocol, respectively.         |             |
               |                                                                        |             |
Xwalk          |   Predicts, validates and visualizes chemical cross-links.             |   0.5       |   http://www.xwalk.org
               |   Applications for distance calculation is names Xwalk. The archive    |             |
               |   contains also a JAVA class called AvgInterface, with which contact   |             |
               |   frequencies can be calculated for predicting binding interfaces.     |             |
               |                                                                        |             |
NACCESS        |   Calculates solvent accessible surface areas. Requires a license,     |   2.1.1     |   http://www.bioinf.manchester.ac.uk/naccess
               |   which is free for academic users.                                    |             |
               |                                                                        |             |
CleftXplorer   |   Primarily written for the geometrical and physicochemical analysis   |   0.1       |   http://code.google.com/p/cleftxplorer
               |   of protein-small molecule binding sites. Holds different type of     |             |
               |   JAVA classes, e.g. an executable for quality threshold alignments    |             |
               |   and RMSD calculations called Superimpose. The archive contains also  |             |
               |   a simple script for RMSD calculation called rmsd.pl.                 |             |
               |                                                                        |             |
colt.jar       |   JAVA library, which holds data structures and algorithms for         |   1.2.0     |   http://acs.lbl.gov/software/colt/colt-download/releases/colt-1.2.0.zip
               |   scientific computing. Is required by CleftXplorer. It is located in  |             |
               |   the colt/lib directory of the downloaded archive file.               |             |
               |                                                                        |             |
cdk-1.0.2.jar  |   Chemistry Development Kit, which holds libraries for structural      |   1.0.2     |   http://sourceforge.net/projects/cdk/files/cdk/1.0.2/cdk-1.0.2.jar
               |   chemoinformatics. Is required by CleftXplorer.                       |             |
               |                                                                        |             |
R              |   Software environment for statistical computing and graphics.         |   2.15      |   http://www.r-project.org
               |                                                                        |             |
PyMOL          |   Molecular visualization software.                                    |   1.5       |   http://www.pymol.org


## Input Files
-------------
The cross-link guided docking protocol requires only four different text files with a specific file
format as input files:
1. receptor.pdb: 3D coordinates of the first protein in Protein Data Bank (PDB) format. The
        PDB file format and in particular the ATOM record describing atomic coordinates of
        proteins can be found at http://www.wwpdb.org/documentation/format33/sect9.html#ATOM. The
        first protein is usually the larger protein and referred to as receptor.
2. ligand.pdb: 3D coordinates of the second protein in PDB format. The second protein is
        usually the smaller protein and referred to as ligand.
3. inter-protein-cross-link.cst: Experimental inter-protein cross-link information in ROSETTA
        constraint file format (See also http://www.rosettacommons.org/manuals/archive/rosetta3.4_user_guide/de/d50/constraint_file.html).
        The constraint file is a space delimited text file and should have following structure for
        this protocol:

        AtomPair {atom name1} {residue number1, chain ID1} {atom name2} {residue number2, chain ID2} FLAT_HARMONIC {x0} {standard deviation} {tolerance}.

    The flat harmonic function guarantees that models are penalized only if the Euclidean distance
    between two cross-linked atoms exceeds 30.0 Å. The function takes a distance measurement dist
    and three parameters that are set for the DSS cross-linking guided docking protocol to x0 = 15,
    tolerance = 15 and standard deviation = 1. The flat harmonic function returns 0,
    if dist - x0 ≤ tolerance and otherwise the square of
    ((dist - x0 - tolerance)/ standard deviation).

4. cross-link.dist: Experimental intra-protein and inter-protein cross-link and mono-link
        information in Xwalk distance file format
        (see also http://www.xwalk.org/cgi-bin/help.cgi#vXLtable), which is a tab delimited text
        file with following information. 
        I.   Incremental cross-link index.
        II.  Arbitrary protein name.
        III. Dash delimited PDB information about first cross-linked atom in the following format:
             residue name – residue number – chain ID – atom name.
        IV.  Dash delimited PDB information about second cross-linked atom in the following format:
             residue name – residue number – chain ID – atom name.
    Mono-links are described with the first three columns only.

## Commandline options for the excutable
----------------------------------------

###docking_protocol

-database {rosetta-database}
    - Path to the ROSETTA database.
-in:file:s relaxed-receptor-ligand.pdb
    - Input PDB structure of the receptor ligand complex in arbitrary conformation.
-constraints:cst_file inter-protein-cross-link.cst
    - ROSETTA constraint file listing all chemical cross-link data.
-docking:low_res_protocol_only
    - Run the docking protocol in low resolution with centroid representation of amino acid side
      chains.
-docking:spin
    - Spin the second docking partner around axes that is span between the centre of mass of the
      receptor and ligand structure.
-docking:randomize1
    - Randomize the orientation of the receptor structure.
-docking:randomize2
    - Randomize the orientation of the ligand structure.
-out:nstruct {N}
    - Number of PDB output files to be generated.

###Xwalk

-infile global-docking_models.pdb
    - Input file corresponding to a ROSETTA docking_protocol output file. 
-dist cross-link.dist
    - List of chemical cross-links in Xwalk file format (see Input file section). 
-max {max-distance}
    - Maximum distance threshold which should correspond to the maximum distance a cross-linker
      can span between two amino acids.
-xSC
    - Removes only side chain atoms of cross-linked amino acids except for CB atoms for SAS
      distance calculations.
-mono
    - Assess the solvent accessibility of mono-linked amino acids.
-out global-docking_models.dist
    - List of chemical cross-links in Xwalk file format (see Input file section). The 5th, 6th and
      7th column hold the sequence distance, Euclidean distance and SAS distance between the
      cross-linked amino acids.

###NACCESS

- 1st parameter of the application is the input file, which corresponds to a ROSETTA
  docking_protocol output file.

###Superimpose

- colt.jar: JAVA library, which holds data structures and algorithms for scientific computing
  (e.g. linear algebra and multi-dimensional arrays). It can be downloaded from
   http://acs.lbl.gov/software/colt/colt-download/releases/colt-1.2.0.zip, and is located in the
   colt/lib directory of the downloaded archive file.
- cdk-1.0.2.jar: JAVA library, which holds libraries from the Chemistry Development Kit for
  structural chemoinformatics. It must be in version 1.0.2 and can be downloaded from
  http://sourceforge.net/projects/cdk/files/cdk/1.0.2/cdk-1.0.2.jar
-dir <dir>
    - Directory with all selected PDB files that satisfy most cross-links and have sufficiently
      large binding interfaces
-chain {chain ID}:{chain ID}
    - Two times chain ID of the receptor
-dock
    - Calculates the RMSD value on the ligand structure excluding coordinates from the receptor
      structure
–ta {translational threshold}:{rotational threshold} 
    - Maximum threshold for the translational and rotational difference between the ligand
      structures of two PDB structures.
-xseq
    - Does not run a sequence alignment to determine the atom mapping between two PDB structures.
-ca
    - Calculates the RMSD value only on C-alpha backbone carbon atoms
–r
    - Outputs only the RMSD value.

AvgInterface

- 1st parameter is a TAR-GNU-ZIP archive file holding the best models from the ROSETTA
  docking_protocol.
- 2nd parameter is a selected PDB file. Ideally THE best model from the ROSETTA docking_protocol.
- 3rd parameter is the path to the NACCESS application
- The output is a PDB file with occupancy values and temperature factors corresponding to the
  absolute and relative frequency, respectively, of an atom to occur at the interface.

