Model a Missing Loop
====================

KEYWORDS: STRUCTURE_PREDICTION LOOPS

Authors: Roland Pache, Michal Sperber, Steven Combs, George Rosenberger  
Last updated: August 2011 (RosettaCon9)

---

This demo shows how missing electron densities of several consecutive residues 
can be modeled using the loop modeling application (loopmodel) and the 
KInematic Closure algorithm (KIC).

The starting structure (1tr2_missing_density.pdb) is based on vinculin (1TR2). 
5 loop residues have been removed (32-36) and should be replaced by the 
sequence VDGKA for loop modeling (simulating missing electron density). 
Afterwards, this PDB structure can be used to model the loop. For this demo, 
the water molecules (HOH) have been removed and the structure was truncated to 
the first 132 residues.

Running the demo
----------------

1.  Insert the new residues into the structure file

    Open the file 1TR2_missing_density.pdb in the text editor of your choice. 
    Search the first gap line (residue 32). Search for the first residue Valine 
    in the file and copy all atoms to the new line. Repeat this step for all 
    other residues (DGKA) and insert the coordinates below Valine. The file 
    should then look like 1TR2_manually_added_dummy_residues.pdb.   Renumber 
    the residues you copied from another place to 32-36 and remove all 
    eventually inserted new lines. Save this file as 
    1TR2_manually_added_dummy_residues_renumbered.pdb.

2.  Create the loop file.

    Create a new file, called 1TR2.loop, and open it in your text editor. 
    Insert the following line:

        LOOP 31 37 37 0 1

    This excerpt from the loopmodel documentation describes the meaning of the 
    6 columns in that line:

        column1  "LOOP":     Literally the string LOOP, identifying this line as a loop
                             In the future loop specification files may take other data.
        column2  "integer":  Loop start residue number
        column3  "integer":  Loop end residue number
        column4  "integer":  Cut point residue number, >=startRes, <=endRes.
        column5  "float":    Skip rate. default - never skip (0)
        column6  "boolean":  Extend loop. Set to 1

    For this example, we select the one residue before and after the loop to 
    have real coordinates that can be used as anchor points by the KIC loop 
    modeling algorithm. The cut point residue number is set to the last loop 
    residue, since it must be inside the loop. The skip rate is set to 0 for 
    this short example (since we want to model this loop) and the extend loop 
    setting is set to true to idealize all bond lengths, bond angles and 
    torsion angles of the loop residues before modeling.

3.  Execution of the algorithm and definition of the flags

    Assuming that Rosetta 3.3 is installed and all paths are set correctly, 
    open your shell and change the directory to the one where the demo files 
    are stored.

		$> cp rosetta_inputs/* .  
    	$> $ROSETTA3/bin/loopmodel.linuxgccrelease -s 1TR2_manually_added_dummy_residues_renumbered.pdb -loops:loop_file 1TR2.loop -loops:remodel perturb_kic -loops:refine refine_kic -ex1 -ex2 -nstruct 1 -loops:max_kic_build_attempts 100 -in:file:fullatom

    Brief descriptions for all the components of this command-line:

        loopmodel.linuxgccrelease : loopmodel application (linuxgccrealease or macosgccrelease)
        -database : path to your Rosetta 3.3 DB
        -loops:input_pdb 1TR2_manually_added_dummy_residues_renumbered.pdb : name of your edited pdb
        -loops:loop_file 1TR2.loop : name of your loops file
        -loops:remodel perturb_kic : kinematic closure based loop modeling low resultion stage (side chains: centroids)
        -loops:refine refine_kic : kinematic closure based loop modeling high resultion stage (side chains: fullatom)
        -ex1 : extra chi rotamers for chi-1 angle (+/- 1 stddev from the optimal rotamer for better loop reconstruction)
        -ex2 : extra chi rotamers for chi-2 angle (+/- 1 stddev from the optimal rotamer for better loop reconstruction)
        -nstruct 1 : number of structures to generate (set to at least 1000 for real application; computationally expensive)
        -loops:max_kic_build_attempts 100 : the maximal number of trials the algorithm should do to find a closed confirmation for the loop (default: 100); can be increased for difficult problems.
        -in:file:fullatom : keep native amino acid side chain confirmations of the non-loop residues (residues within 10 Angstrom of the loop will be remodeled by default)

4.  Analysis of the results

    If you set -nstruct > 1, look at the Rosetta energy score in the standard 
    output and identify the model with the lowest energy. Compare them visually 
    using your favorite molecular visualization application (pymol, etc).
    Sample output files (incl. visualization using pymol) can be found in the 
    output directory.
    
