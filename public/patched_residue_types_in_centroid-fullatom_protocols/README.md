# Tutorial for using modified residue types in centroid-level applications.

KEYWORDS: CORE_CONCEPTS NONCANONICALS

In this tutorial, you will describe how to wrangle modified residue types through combined centroid/fullatom protocols in Rosetta.

##  Background

1. Many Rosetta protocols use a centroid phase followed by a fullatom phase.  However, the centroid residue set is compatible with the canonical 20 amino acids, and very little else.  For example, ligands and post-translational modifications are incompatible with centroid mode; Rosetta generally crashes when trying to switch these residues to centroid (or silently drops the post-translational modification).

2. You have been provided with a PDB file containing a two phosphoresidues as input (one p-SER and one p-TYR). This tutorial will show you how run the docking\_protocol application using this PDB model as input without Rosetta exiting with an error.

## Clean up the PDB file

1. This is an NMR model. Rosetta will read this as a giant complex of superimposed structures, which is bad. You need to find the line labeled "ENDMDL" and delete all lines below it.
    ```
    gunzip -c starting_files/2lax.pdb.gz > rosetta_inputs/2lax_edited.pdb
    <your favourite text editor> rosetta_inputs/2lax_edited.pdb
    ```

2. Rosetta wants the phosphorylated residues to be named the same as their canonical countertypes. Open the pdb file rosetta\_inputs/2lax\_edited.pdb for editing. For residue 202 and residue 206, rename the three-letter code "TPO" to "TYR" and "SEP" to "SER". Rosetta will identify these as phosphorylated based on their atom names.

    - note: If at this point you were to run any protocol that utilizes centroid mode, you would see this error when Rosetta crashed:
        ```
        can not find a residue type that matches the residue TYR_p:phosphorylatedat position 38

        ERROR: core::util::switch_to_residue_type_set fails
        ```

## Create centroid residue parameter patch files

For modified residues (e.g. phosphorylation or acetylation), Rosetta uses "patch" files to modify the pre-existing residue parameter files. Unfortunately, a centroid-level patch file for phosphorylated residues does not exist. You need to create two, one for P-TYR and one for P-SER.
1. Go to the database directory ```<my_rosetta_directory>/rosetta_database/chemical/residue_type_sets/centroid/patches/``` and copy the file tyr_sulfated.txt to tyr_phosphorylated.txt. Open that file for editing. Change these lines:
```
NAME sulfated
TYPES SULFATION
```
to this:
```
NAME phosphorylated
TYPES PHOSPHORYLATION

NOT VARIANT_TYPE PHOSPHORYLATION 
```
Now copy this new file tyr_phosphorylated.txt to ser_phosphorylated.txt. Open ser_phosphorylated.txt for editing. Change this line:
```
AA TYR
```
to this:
```
AA SER
```
Now, we have patch files for our centroid-level phosphorylated residues. All we have to do now is point Rosetta to these files. Open the file  ```<my_rosetta_directory>/rosetta_database/chemical/residue_type_sets/centroid/patches.txt```. Add these two lines to the bottom of the file:
```
patches/tyr_phosphorylated.txt
patches/ser_phosphorylated.txt
```
 Now run the docking_protocol application like this:
 (where `$ROSETTA3`=path-to-Rosetta/main/source)
```
$> $ROSETTA3/bin/docking_protocol.default.linuxgccrelease -s rosetta_inputs/2lax_edited.pdb
```

This application should now work correctly, converting the input pdb coordinates to a centroid model, performing rigid-body docking, then converting back into a full-atom model ( which should now *correctly convert the phosphorylated residue types* ) before performing a final docking rotamer packing optimization (see the docking_protocol documentation for more information). 
