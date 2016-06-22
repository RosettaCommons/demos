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

### 1. Standart CS Rosetta

The required inout files are in `./input_files`.

1. In the cs-rosetta.flags, change the path to your vall database
2. Run the picker:

   ```
   $> <path-to-Rosetta>/main/source/bin/fragment_picker.default.linuxgccrelease @cs-rosetta.flags -in::file::vall 
    $ROSETTA3/main/database/sampling/small.vall.gz
```

### 2. CS Rosetta with TALOS L1 rama
```
   $> <path-to-Rosetta>/main/source/bin/fragment_picker.default.linuxgccrelease @cs-rosetta_TALOS_rama.flags -in::file::vall
    $ROSETTA3/main/database/sampling/small.vall.gz
```

