#Protein-Protein-Docking

This tutorial will introduce you the main steps required to predict how two proteins interact.

######Input files
* COL_D.pdb
* IMM_D.pdb

##1. Preparation

1. First, you have to check that there are no missing residues in any of the structures. If a protein chain is interrupted, the docking protocol might start moving around two parts of the same chain, instead of moving the two chains to be docked.

2. Then, you have to combine the protein structures into a single file. 

        $ cp COL_D.pdb combined.pdb
        $ cat IMM_D.pdb >> combined.pdb
      
 You can open the new pdb file in Pymol to see, whether both proteins are present.

3. During a docking simulation, side chains conformations will be optimized. This optimization can improve the energy, even if the two proteins hardly interact and they make it harder to indentify the best-docked model by comparing energies.

 To minimize the noise from docking-independent changes in energy, you should allways perform a side-chain optimization step (**packing**).

        $> ../../../main/source/bin/docking_prepack_protocol.linuxgccrelease \
        -s combined.pdb \
        -ex1
        -ex2
        
  You will see several output files. The prepacked combined pdb is combined_0001.pdb. Change that name:
  
        $> mv combined_0001.pdb combined_ppk.pdb
        
   This is your input pdb file for the docking run.
   
##2. Run a docking simulation

### Global docking
Here is a minimal commandline to start a docking simulation (here we include also the known complex structure, 1V74, as a reference):

     $> ../../../main/source/bin/docking_protocol.linuxgccrelease \
     -s combined_ppk.pdb
     -ex1
     -ex2
     -nstruct 2
1. Run the docking protocol  
     
 Rosetta will internally connect the centers of the two chains with a so-called jump (see [Fold Tree]()). Along this jump the chains are being pulled together (slide into contact). The Monte Carlo moves, which are selected randomly, are:   
 
 * Translations (in x,y or z direction)
 * Rotations (around x,y, z axis)

2. Open the ouput files combined\_ppk\_0001.pdb and combined\_ppk_0002.pdb in Pymol. The two chains should now be in contact and you can see the different orientations that the chains have now.

### Local refinement only:

* If you already know where two proteins interact, but not exactly how, you can do a **local docking** by including:

        -docking_local_refine
        
### More stuff that is good to know
1. **Dock a group of chains** against a protein. There is an options that allows you tell Rosetta to leave a group of chains unchanged, relative to each other:

		$ -partners LH_A
		
 this would dock chain A onto the complex of chain L and H

2. **Symmetric docking**  
Homomeric protein complexes are often symmetric. Rosetta provides the opportunity to use symmetry information during docking (and many other applications). Find the symmetry documentation and demos for mor information.

> **NOTE:** Take a closer look at the [docking documention](https://www.rosettacommons.org/docs/wiki/application_documentation/docking/docking-protocol). There is a lot more information to be found, on how to manipulate the behaviour of this protocol.

##3.Analysis
You can repeat the run with e.g. *-nstruct 20000* and outputting silent files. Then you could extract the best model and create a score vs. rmsd plot.


Look at [Analysis](Analysis) for informationon how to analyze your results.


