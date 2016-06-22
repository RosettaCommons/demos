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


