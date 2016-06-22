Fragment Picking for CS-Rosetta
===============================

Presenting author: Dominik Gront (dgront at chem dot uw dot edu dot pl)  
Protocol name: fragment picker : CS-Rosetta style  
Brief description: The protocol substitutes CS-Rosetta application

Source code location
--------------------

* Fragment picker is located in Rosetta/main/source/src/core/fragment/picking
* The source lives in Rosetta/main/source/src/apps/public/fragmentpicker.cc

Running the protocol 
----------------------------

The required inout files are in `./input_files`.

1. In the cs-rosetta.flags, change the path to your vall database
2. Run the picker:

   ```
   $> <pth-to-Rosetta>/main/source/bin/fragmentpicker.default.linuxgccrelease @cs-rosetta.flags
   ```
