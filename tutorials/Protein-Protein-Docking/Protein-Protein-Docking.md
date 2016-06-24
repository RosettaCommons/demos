Protein-Protein Docking
=======================
KEYWORDS: DOCKING GENERAL STRUCTURE_PREDICTION    
Written by by Sebastian Rämisch (raemisch@scripps.edu). Edited by Shourya S. Roy Burman (ssrb@jhu.edu)   
Edited Jun 24 2016   

[[_TOC_]]

Summary
-------
Rosetta can be used to predict the bound structure of two proteins starting from unbound structures. By the end of this tutorial, you should be able to understand:

* How to prepare structures for docking
* How to locally dock two proteins
* Hoe to refine an already docked structure
* How to dock two proteins whose interface region is unknown
* How to dock flexible proteins
* How to dock a flexible peptide to a protein
* How to dock symmetric proteins
* How to analyse the best docked model

Navigating to the Demos
-----------------------
The demos are available at `$ROSETTA3/demos/tutorials/Protein-Protein-Docking`. All demo commands listed in this tutorial should be executed when in this directory. All the demos here use the `linuxgccrelease` binary. You may be required to change it to whatever is appropriate given your operating system and compiler.

Compare your output files to the ones present in `output_files/expected_output`.

Preparing Structures for Docking
--------------------------------
This tutorial will introduce you the main steps required for predicting the bound structure of two interacting proteins starting from the unbound structures. For this example, we will dock Colicin-D with its inhibitor, IMM. You are provided with the two refined input files `COL_D.pdb` and `IMM_D.pdb`, and a native file `1v74.pdb` in the folder `input_files`.

To prepare structures for docking, be sure to refine them as described in [[preparing inputs tutorial|input_and_output#controlling-input_preparing-a-structure-by-refinement]].

Local Docking
-------------
Rosetta is most accurate when docking locally. **In local docking, we assume that we have some information about the binding pockets of the two proteins.** First, we must manually place the two proteins (within ~10 Å) with the binding pockets roughly facing each other as shown in this figure:
![unbound](images/COL_IMM_unbound.png)

We will pass the following options to indicate that i) chain B is being docked to chain A, ii) we want to randomly perturb the ligand of the input strucure (chain B) by 3 Å translation and 8° rotation before the start of every individual simulation, and iii) we want to spin the ligand around the receptor (chain A).

```
-partners A_B
-dock_pert 3 8
-spin
```

We will also compare the input with the bound structure 1v74.pdb by passing it as native. Now to start docking, run:

    $>$ROSETTA3/main/source/bin/docking_protocol.linuxgccrelease @flag_local_docking
    
This should take ~30 seconds to run and produce a structure file and a score file in `output_files`. The structure file might not dock well with just one attempt (as the `-dock_pert` flag will move the ligand away more often than towards). Make sure you use `-nstruct 500` or more in production runs.

Local Refinement of Docked Structures
-------------------------------------

3. During a docking simulation, side chains conformations will be optimized. This optimization can improve the energy, even if the two proteins hardly interact and they make it harder to indentify the best-docked model by comparing energies.

 To minimize the noise from docking-independent changes in energy, you should allways perform a side-chain optimization step (**packing**).

        $> ../../../main/source/bin/docking_prepack_protocol.default.linuxgccrelease -s combined.pdb
        
  You will see several output files. The prepacked combined pdb is combined_0001.pdb. Change that name:
  
        $> mv combined_0001.pdb combined_ppk.pdb
        
   This is your input pdb file for the docking run.
   
##2. Run a docking simulation

### Global docking
Here is a minimal commandline to start a docking simulation (here we include also the known complex structure, 1V74, as a reference):

     $> ../../../main/source/bin/docking_protocol.linuxgccrelease -s combined_ppk.pdb -ex1 -ex2 -nstruct 2
     
1. Run the docking protocol  
     
 Rosetta will internally connect the centers of the two chains with a so-called jump (see [[Fold Tree|fold_tree]]). Along this jump the chains are being pulled together (slide into contact). The [[Monte Carlo|Core_Concepts#monte-carlo-sampling]] moves, which are selected randomly, are:   
 
 * Translations (in x,y or z direction)
 * Rotations (around x,y, z axis)  
 
 > You can test the docking protocol by running with:  
 > -nstruct 20000  
 > -out:file:silent\_struct_type binary  
 > -out:file:silent dock.out  
 > -native 1V74.pdb                 
 
 1V74.pdb is the native complex.

2. Open the ouput files combined\_ppk\_0001.pdb and combined\_ppk_0002.pdb in Pymol. The two chains should now be in contact and you can see the different orientations that the chains have now.

### Local refinement only:

If you already know where two proteins interact, but not exactly how, you can do a **local docking** by including:

        -docking_local_refine
        
This will move the chains only by a few degrees and Angstroms. An example case would docking antibodies to antigens, when the sequence of th epitope is known. The required steps would be:

1. Fix chain chain breaks, if present. 
2. Prepare a combined file, where the two chains are in conformation that you anticipate to be close to the native one.
3. Run a docking simulation using the above option. The required nstruct is usually smaller than for global docking, as the sampling space is much smaller. 5000 is a reasoable number to start with. 

	
### More stuff that is good to know
1. **Dock a group of chains** against a protein. There is an options that allows you tell Rosetta to leave a group of chains unchanged, relative to each other:

		$ -partners LH_A
		
 this would dock chain A onto the complex of chain L and H

2. **Symmetric docking**  
Homomeric protein complexes are often symmetric. Rosetta provides the opportunity to use symmetry information during docking (and many other applications). Find the symmetry documentation and demos for mor information.

> **NOTE: Take a closer look at the [docking documention](https://www.rosettacommons.org/docs/latest/application_documentation/docking/docking-protocol). There is a lot more information to be found, on how to manipulate the behaviour of this protocol.**

##3.Analysis
You can repeat the run with e.g. *-nstruct 20000* and outputting silent files. Then you could extract the best model and create a score vs. rmsd plot.


Look at [[Analysis|Analysis]] for informationon how to analyze your results.


