Fragment Picking for CS-Rosetta
===============================
KEYWORDS: STRUCTURE_PREDICTION UTILITIES

Presenting author: Dominik Gront (dgront at chem dot uw dot edu dot pl)  
Protocol name: fragment picker : CS-Rosetta style  
Brief description: The protocol substitutes CS-Rosetta application

Source code location
--------------------

* Fragment picker is located in Rosetta/main/source/src/core/fragment/picking
* The source lives in Rosetta/main/source/src/apps/public/fragmentpicker.cc

Running the protocol 
----------------------------

### 1. Standart CS Rosetta

The required inout files are in `./input_files`.

1. In the cs-rosetta.flags, change the path to your vall database
2. Run the picker:

   ```
   $> $ROSETTA3/bin/fragment_picker.default.linuxgccrelease @cs-rosetta.flags -in::file::vall $ROSETTA3_DB/sampling/small.vall.gz
```

**IMPORTANT**
The *small.vall.gz* used here for fragment picking is only used to speed up the demo. You have to change this to the vall database on your system!

### 2. CS Rosetta with TALOS L1 rama
```
   $> $ROSETTA3/bin/fragment_picker.default.linuxgccrelease @cs-rosetta_TALOS_rama.flags -in::file::vall $ROSETTA3_DB/sampling/small.vall.gz
```
**IMPORTANT**
The *small.vall.gz* used here for fragment picking is only used to speed up the demo. You have to change this to the vall database on your system!


