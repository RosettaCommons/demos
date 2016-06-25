Analyzing Interface Quality
===========================

KEYWORDS: ANALYSIS INTERFACES STRUCTURE_PREDICTION

Outline
-------
You have been provided with two PDBs to use for demonstration: 1U6E, a simple homodimer, and 3R2X, Sarel's hemaglutinin/designed protein structure.  In the latter case, we will calculate the interface between the designed protein and the HA - the interface between HA chains is unimportant.

Preparing the inputs
--------------------
We usually score the input pdbs to make sure they are able to be read by Rosetta and to replace any missing sidechain atoms.
From the main demo directory run the following command to score 

    $> cd rosetta_inputs
    $> $ROSETTA3/bin/score_jd2.default.macosgccrelease -s ../starting_files/*.pdb.gz -no_optH false -ignore_unrecognized_res -out:pdb

There are a few options here that need to be described:
* `-no_optH false`: This will make rosetta consider Q,N,H ring flips (this helps remove some buried unsatisfied polar atoms).
* `-ignore_unrecognized_res`: This will drop any residues and waters from the PDB that rosetta does not recognize.
* `-out:pdb`: This forces the score application to output a pdb to use later

This will also output file called score.sc which is the scored input structures.
Now we'll rename these PDBs to make this a bit easier to keep track of:

    $>  mv 1u6e_0001.pdb 1u6e_scored.pdb
    $>  mv 3r2x_0001.pdb 3r2x_scored.pdb

Running InterfaceAnalyzer
-------------------------
Now we have all of our inputs set to go to run InterfaceAnalyzer.
If you are interested you should read the full documentation for this application located in `/path/to/rosetta/doc/apps/public/analysis/interface_analyzer.dox`:

    $> cd ../example_output

First we are going to look at the 1U6E homodimer interface.  It only has two chains so the default way of defining an interface will work. 
You will need to modify the options files to contain the path to your rosetta database.

To analyze this interface run the following commands:

This one will use only the input sidechains from the input structure and not repack them.  This probably only a good idea if you have already designed/minimized this structure, but is provided as an example here anyway

	$> $ROSETTA3/bin/InterfaceAnalyzer.default.macosgccrelease -s ../rosetta_inputs/1u6e_scored.pdb @../rosetta_inputs/no_pack_input_options.txt

This one will repack the input sidechains from the input structure before analyzing anything about the interface. This is a better idea in this case because we are using a raw pdb as input.

    	$> $ROSETTA3/bin/InterfaceAnalyzer.default.macosgccrelease -s ../rosetta_inputs/1u6e_scored.pdb @../rosetta_inputs/pack_input_options.txt

You can look in the options files for a full description of each of the command line flags.

Now to analyze the 3R2X interface between the designed protein (Chain C) and Hemaglutinin (Chains A & B)  we need to consider chains A and B as one monomer and chain C as the binding partner. To do this we have the option available to keep any given chains together in this case the option needed is 

    -fixedchains A B

We will add this to the command lines given above; so now run these:

    $> $ROSETTA3/bin/InterfaceAnalyzer.default.macosgccrelease -s ../rosetta_inputs/3r2x_scored.pdb -fixedchains A B @../rosetta_inputs/no_pack_input_options.txt
    $> $ROSETTA3/bin/InterfaceAnalyzer.default.macosgccrelease -s ../rosetta_inputs/3r2x_scored.pdb -fixedchains A B @../rosetta_inputs/pack_input_options.txt

Now in the current directory (`example_output`) we have two different score files, one for when we repacked the interface (`pack_input_score.sc`), and one for when we kept the input rotamers fixed (`no_pack_input_score.sc`)

From here we can get some important info.  We will only concentrate on a few here. For a full description of what everything means see the documentiation page for this application: `/path/to/rosetta/doc/apps/public/analysis/interface_analyzer.dox`

Looking at the results
----------------------
Open the output files in some form of text editor or spreadsheet application.  Below we describe how to find a few important output values.

Binding energy (dG_separated): this is computed deltaG of binding. Notice that in no_pack_input the value for 3R2X is unrealistic.  This is probably due to clashes in the input structure, which is why it is a good idea to relax your structure somehow before calculating these values.

Number of buried unsatisfied polar atoms (delta_unsatHbonds): note that the number is higher for 1U6E when you pack the input, this sometimes happens when rosetta tries to relieve clashes and thus leaves some polars without a hydrogen bond partner

Packing score (packstat): This is a measure of how well packed the interface is with 0.0 being as poor as possible and 1.0 being perfect shape complementarity (usually values above 0.65 are good)

Buried Surface Area (dSASA_int): change in exposed surface area upon formation of an interface.

Binding energy per unit area (dG_separated/dSASAx100): This is the dG_separated binding energy divided by the total interface surface area (dSASA_int). We multiply by 100 to scale it up the value so it fits better in the score file.  We like using this to make sure that rosetta is making high quality contacts instead of making a lot of low quality contacts across the interface.  Usually values below -1.5 are pretty good. 

