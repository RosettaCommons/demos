Rosetta Membrane Framework: Visualizing the Membrane in PyMOL  
==============================================================

## Data
Author: Rebecca F. Alford (rfalford12@gmail.com)
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)
Last Updated: September 2014
Rosetta Revision #<XXX>
Citation: Modeling membrane proteins in Rosetta 3

Corresponding Documentation Link: 

## Application Description
This demo describes how to combine the membrane framework with the PyMOL viewer 
tool to visualize membrane geometry. Here we describe a specific application for
visualization. 

## Executable/Script
The pymol viewer can be used directly with the membrane framework using view_membrane_protein<.exe>
However, these set of flags can be used with any membrane framework protocol. 

## Generating Inputs
Viewing the membrane planes only requires a spanfile as input. 

  1. Generating a Spanfile
  A spanfile describing transmembrane spanning regions can be generated using the OCTOPUS server
  (http://octopus.cbr.su.se/). This file must be converted to a Rosetta spanfile format using octopus2span.pl

    cd mpframework-ddG/scripts/
    ./octopus2span.pl octopus_pred.out > spanfile.txt

## Useful Scripts
This demo contains a script directory with: 
  - octopus2span.pl: Convert OCTOPUS topology prediction to Rosetta spanfile format

## Running the Application

## Example Outputs

