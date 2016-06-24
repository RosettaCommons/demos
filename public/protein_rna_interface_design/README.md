#Protein-RNA Interface Design

KEYWORDS: NUCLEIC_ACIDS INTERFACES RNA

This tutorial explains how to use rosetta scripts to redesign amino acids in a protein-RNA interface.

This demo was last modified by Parisa Hosseinzadeh in 2016.

This tutorial provides for you a simple script (protein_rna_interface_design.xml) to design protein residues around certain RNA residues. Luckily, Rosetta can now handle RNA residues as well as DNA residues but in the script folder, and in the Appendix, you can find older scripts and requirements for working with RNA. These files remained for historic reasons.

Step1: Prepare PDB
In this tutorial we are using the 1a9n.pdb file. It is in the starting_files directpry or you can directly download it from pdb. Now, run this command to clean the PDB and save one of the dimers interacting with RNA:

```
> <path-to-Rosetta>/Rosetta/tools/protein_tools/scripts/clean_pdb.py starting_files/1a9n.pdb AB
```
The file is also provided for you in the starting_file directory.

Now, clean the RNA chain by running:
```
> scripts/clean_pdb_rna_prot.py starting_files/1a9n.pdb 
```
This will only outputs the RNA part of the sructure in a cleaned manner. open the file with an editor and choose all of chain Q. 
```
> cp 1a9n_AB.pdb 1an9_RNA.pdb
```
Add the selected chain to the end of 1a9n_RNA.pdb. Open the file in pymol and look at it. You should see one RNA chain and two protein chains that interact with the RNA chain as a dimer. There are other ways to prepare this file, for example using pymol that you may choose. 

Step 2: Running the script

Take a look at the protein_rna_interface_design.xml. It uses a collection of residue selectors and task operations to tell Rosetta to only design protein around 10A of certain residues in RNA. For more information on the rosetta scripts, check the documentation and tutorials. You need to change these numbers for your purpose.

If you don't have the prepared pdb, copy it to your directory:

```
$> cp starting_files/1a9n_RNA.pdb .
```

Now run the script: (where `$ROSETTA3`=path-to-Rosetta/main/source)
```
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease -parser:protocol scripts/protein_rna_interface_design.xml -in:file:s 1a9n_RNA.pdb -ex1 -ex2
```
You should see a file named 1a9n_RNA_0001.pdb and score.sc. Take a look. You can compare them to the ones pre-generated in the output directory.

-------------------------------------------------------------------
APPENDIX

Step 1:
	For complicated reasons, rosetta can only handle DNA and proteins in the same pdb.  Since we're dealing with RNA and protein, this first script makes changes in the rosetta database to allow use to work with RNA and protein.  Note that once you run this script you'll have to revert to trunk rosetta database in order to deal with DNA and protein.
	The script is located in the scripts folder and is named fix_rna_db.sh.  When you run this command, include the location of your rosetta database WITHOUT a trailing slash.  Proper command line is as follows :
		location_of_script/fix_rna_db.sh location_of_database (Remember, no trailing slash)

Step 2:
	The pdb file won't work with rosetta as is, so you need to make a couple of modifications.  The first modification is completed by the script named clean_pdb_rna_prot.py.  This script places all amino acids in chain A and all nucleic acids in chain B.  The script is run by the command python (pathname to scripts)/clean_pdb_rna_prot.py
	The second modification to the pdb is to rename the nucleic acid residues to the names that rosetta uses.  This scripts is called with the command pathname_to_script/make_pdb_rosetta_compatible.py

Step 3:  Congratulations!  You're ready to run rosetta!  We use an xml scripting language called rosetta_scripts to run specific protocols.  See the command script and the protein_rna_interface_design.xml file to see exactly how the protocol is run.  There's a command script as well called command.sh.  Use this script to run rosetta_scripts and the .xml file.  Note that to change protein-RNA design to your favorite pdb, you need to change the xml (see command and xml for more details).

Step 4: ???

Step 5: Dinner!
