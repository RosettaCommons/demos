Antibody Docking
================

KEYWORDS: ANTIBODIES DOCKING

The entire workflow for this demo should be described in this file.
It should describe an entire work flow, with command lines, tested if possible.

Authors:
* Jianqing Xu (xubest at gmail dot com)
* Christine Tinberg
* Jeff Gray
* Angela Loihl

Demo files
----------

`starting_files/`
* directory in which the raw input files are given - these are provided for you and serve as the base for your tutorial

`rosetta_inputs/`
* empty in this demo
* directory in which the modified starting files should be placed which will be used as inputs to Rosetta.  You may need to make modifications like stripping extra chains from the input PDB; store the modified PDB here and leave the unaltered one in starting_files 

`scripts/`
* empty in this demo
* python scripts, shell scripts, used in the workflow
* awk, grep, sed, etc. command lines

`README.md`
* A prose or list description of how to perform the protocol

`FOR_AUTHORS.txt`
* A description for the demo creators of what their demo should achieve.
* Most of what starts in this file should end up in the README file as well.

Running the demo
----------------

1.  The first thing a user should notice is that there’s an antibody 
    modeler in Rosetta 3, but still under development.  As of August 
    2011, Rosetta3 should be used for camelid antibody modeling, but 
    Rosetta2 should be used for other antibody modeling and for antibody 
    docking via SnugDock.  The current stable version of Rosetta2 
    (Rosetta++) is the released Rosetta-2.3.1

2.  Obviously, you need structures of both antibody and antigen in order 
    to do antibody-antigen docking. If you don't have antibody structures, 
    but have antibody sequences, you can use Gray lab antibody homology 
    modeling server (http://antibody.graylab.jhu.edu/) and input the sequence 
    of the light and heavy chain. You will get best 10 structures.

    If you want to manually run the scripts yourself, you can download 
    the scripts source code, example, and instructions from the link below: 
    https://svn.rosettacommons.org/source/trunk/antibody/.
    Again, please realize that the H3 loop modeling is still from Rosetta++.

3.  Download the released version of Rosetta++: 
    https://svn.rosettacommons.org/source/branches/releases/rosetta-2.3.1/ 

    For people outside of the rosetta community: 
    https://www.rosettacommons.org/software/academic/2.3.1/RosettaSnugDock-2.3.1.tgz

4.  Compile rosetta:
    ```
    tar –zxvf RosettaSnugDock-2.3.1.tgz (if you download the second link)
    cd rosetta++
    scons mode=release –j12    (assuming you can use 12 CPUs)
    ```

5.  The SnugDock example in Rosetta++ can be found at:
    https://svn.rosettacommons.org/source/branches/releases/rosetta-2.3.1/example/
    * Besides the original example shown above, we made a new example in the current directory.
      Please be careful with different flags used in the command line.
      The documentations of SnugDock options can be found at:
      http://www.rosettacommons.org/guide/SnugDock.
      Please also be careful with each "paths.txt" file.
    * Do Ensemble Prepack:
      ```
      cd ./PrePack_input # you need AB_model*.pdb, ABRM.fab, ABRM.pdb, ABRM.unbound.pdb, Antigen.pdb, pdblist1, pdblist2
      cd ../Prepack # you need `paths.txt` and `prepack.bash` file available
      ./prepack.bash
      ```
      After prepack, you will see `*.ppk` files in the PrePack_input directory.
      Your original pdblist1 and pdblist2 files will be modified as well, please see the README file inside that folder.
    * Do SnugDock+Ensemble:
      ```
      cd ../SnugDock # you need EnsembleDock_plus_SnugDock.bash and paths.txt file
      ./EnsembleDock_plus_SnugDock.bash # please realize that the example we used here is a little slow, due to the protein size
      ```

6.  Some extra information you may need, besides the example tutorial linked above:

    * Make fab file:
      The CDR loops of antibody should point to the antigen.
      By specifying the antibody loops in the fab file, one can reduce the computational cost for global docking.
      The scripts to make fab file can be found at:  
      https://svn.rosettacommons.org/source/branches/releases/rosetta-2.3.0/rosetta_scripts/docking/  
      Run makefab.pl on your pdb of choice.
      ```
      ./makefab.pl `input pdb` <heavy and/or light chain i.e. HL>
      ./makefab.pl AB_model1.pdb HL
      ```

    * Make `FR02.pdb` complex file:
      Use pymol to open both the antibody and antigen in one session and save both into one pdb file.
      In the example: `ABRM.pdb`.
      It's better to point the antibody CDRs to the antigen, and keep them at a certain distance.

