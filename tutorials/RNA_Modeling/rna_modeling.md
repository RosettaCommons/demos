# RNA 3D Modeling Protocol
Tutorial by Fang-Chieh Chou, Dec. 2015.
Adapted by Kristin Blacklock, Jun. 2016.

## Setup

*Assumes you have already downloaded and built the current version of Rosetta.*

### Set up the Rosetta environment 

Set up the necessary environment variables. Edit your ~/.bash_profile and add the following lines:

```
#Rosetta Paths
-export ROSETTA=$HOME/Rosetta
export ROSETTA3_DB=$ROSETTA/main/database
-export ROSETTA_TOOLS=$ROSETTA/tools
-export PATH=$PATH:$ROSETTA/main/source/bin
export RNA_TOOLS=$ROSETTA_TOOLS/rna_tools/
export PATH=$PATH:$RNA_TOOLS/bin/
export PYTHONPATH=$RNA_TOOLS/bin/
```

Here export `ROSETTA=$HOME/Rosetta/" points to your Rosetta root path. If you installed Rosetta elsewhere (other than `~/Rosetta`) you need to modify it.

Restart your shell window. As a quick test, try running
```bash
rna_helix.py -h
```
It should print a help message for the `rna_helix.py` application.

### Create a fasta file containing the sequence

Example:
```bash
> SAM I-IV riboswitch, RNA puzzle 8
ggaucacgagggggagaccccggcaaccugggacggacacccaaggugcucacaccggagacgguggauccggcccgagagggcaacgaaguccgu
```
###### Note: Must have **lowercase** letters.
Here we name this file as "fasta". 

### Create a secondary structure file

This file contains the sequence and and the parenthesis notation for secondary structure. You can use '[]' and '{}' in addition to represent pseudoknots.

Example:
```bash
(((((....((((....))))(((..((((((.[[[[[.))).))).)))..(((((....)))))))))).((((....))))......]]]]].
GGAUCACGAGGGGGAGACCCCGGCAACCUGGGACGGACACCCAAGGUGCUCACACCGGAGACGGUGGAUCCGGCCCGAGAGGGCAACGAAGUCCGU
```
Here we name this file as "secstruct"

## Threading

1. Open the PDB `2gis.pdb` file in PyMOL. This will be our input PDB.

2. Cut out the useful region from the input PDB in PyMOL and copy that structure to the working folder (an example has already been done for you - see the `2gis_cut.pdb` file).

3. Convert it to proper Rosetta format
```bash
$> make_rna_rosetta_ready.py 2gis_cut.pdb
```

4. Align the sequence for threading in a fasta file  
Copy the fasta file of the target RNA to the working folder. Convert the uppercase letters to lowercase (should make this automatic in the future). Generate the sequence of the template by running
```bash
$> pdb2fasta.py 2gis_cut_RNA.pdb >> fasta
```

Edit the `fasta` file to adjust the alignment. Use '-' to fill missing residues.

Example:
```bash
> SAM I-IV riboswitch, RNA puzzle 8
ggaucacgagggggagaccccggcaaccugggacggacacccaaggugcucacaccggagacgguggauccggcccgagagggcaacgaaguccgu

>2gis_cut_RNA.pdb  A:1-8
----caga----------------------aaug--------------------------------------------------------------
```
This file is provided for you, called `fasta.aligned`.

5. Run the `rna_rethread` application

Command line example:
```bash
$> rna_thread.default.linuxgccrelease -fasta fasta.aligned -s 2gis_cut_RNA.pdb -o core.pdb
```
This script mutates the nucleotides in 2gis_cut_RNA.pdb to the SAM I-IV riboswitch identities that they match to in the fasta.aligned file. The output of this command should be new threaded model named `core.pdb`.

Open the output pdb file in PyMOL and make sure it looks correct.

## RNA Helix

Sometimes you may want to put in a rigid, ideal helix for the FARFAR modeling. In this case, you can use the `rna_helix.py` code to build it.

Example command:
```bash
$> rna_helix.py -seq ggac gucc -resnum 1-4 45-48
```

This command build an ideal RNA duplex named `helix.pdb` with sequence and residue number as:
```
1  G G A C 4
48 C C U G 45
```

You can now use this helix as one of your input pdb models in FARFAR.

## FARFAR

Examine the contents of the `README_SETUP` script:
```bash
rna_denovo_setup.py -fasta fasta -secstruct_file secstruct \
 -tag full1 \
 -working_res 1-96 \
 -s core.pdb P2-loop.pdb P4-loop.pdb P4p-loop.pdb pseudoknot.pdb \
 -cycles 10000 \
 -fixed_stems \
 -filter_lores_base_pairs_early -autofilter -filter_chain_closure_halfway
```

The meaning of the options are:

| Flag            | Meaning |
| :-------------- | :---------------------------------------------------------------------------- |
| -fasta          | The fasta file for target RNA.                                                |
| -secstruct_file | The secondary structure file for target RNA.                                  |
| -tag            | Name tag for the output files.                                                |
| -working_res    | Residues to be built. Using the numbering in fasta.                           |
| -s              | Input pdb files. These input residues are assumed to be fixed during FARFAR.  |
| -cycles         | Total FARFAR cycles to be executed before it ends.                            |


Run the `README_SETUP` file to generate the Rosetta command and params files:
```bash
$> source README_SETUP
```

This should generate a `README_FARFAR` script, which you can run with:
```bash
$> source README_FARFAR
```

To get the output pdbs, use the following command:
```bash
$> extract_lowscore_decoys.py full1.out 9
```
This will extract nine decoys with the lowest scores from the output pdb file. You can then view those models with PyMOL.

