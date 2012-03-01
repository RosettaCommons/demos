#!/bin/bash

#usage: ./command.sh <starting PDB> <path_to_mini>

mini_path="/work/shawnyu/mini"
minirosetta_path="/work/shawnyu/rosettacon_demo/minirosetta_database/"
release="static.linuxiccrelease"

$mini_path/bin/rosetta_scripts.$release -s $1 -use_input_sc -nstruct 1 -database $minirosetta_path -ex1 -ex2 -parser:protocol protein_rna_interface_design.xml 
# -s $pdb # Input PDB file
# -use_input_sc # Tells Rosetta not to force sidechain rotamers into Dunbrack rotamers (i.e. uses rotamer of input PDB instead) # -nstruct 1 # Generates one output structure; sufficient for this design without rigid body movements
# -parser:protocol protein_rna_interface_design.xml # This xml file runs the main packing and design protocol; see xml for documentation for specific design steps used.  Note that the residue numbers for the RNA side are hardcoded and need to be changed when using different input PDBs (see xml).

