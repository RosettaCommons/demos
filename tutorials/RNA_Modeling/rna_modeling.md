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
GGAUCACGAGGGGGAGACCCCGGCAACCUGGGACGGACACCCAAGGUGCUCACACCGGAGACGGUGGAUCCGGCCCGAGAGGGCAACGAAGUCCGU
```
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

1. Open the PDB `2gis.pdb` in this directory in PyMOL. This will be our input PDB.

2. Cut out the useful region from the input PDB in PyMOL and copy that structure to the working folder (this has already been done for you - see the `2gis_cut.pdb` file).

3. Convert it to proper Rosetta format
```bash
make_rna_rosetta_ready.py 2gis_cut.pdb
```

4. Align the sequence for threading in a fasta file  
Copy the fasta file of the target RNA to the working folder. Convert the uppercase letters to lowercase (should make this automatic in the future). Generate the sequence of the template by running
```bash
pdb2fasta.py 2gis_cut_RNA.pdb >> fasta
```
