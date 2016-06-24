# Prepare PDBs for Rosetta

KEYWORDS: UTILITIES GENERAL

## Introduction
For the most part Rosetta does a pretty good job at reading
structures, however with unusual chemical modifications, sometimes
Rosetta runs into problems. Here are some of the problems with pdb
files that we have encountered that currently cause problems with
Rosetta:


## 0: Check the PDB
If you run into problems with processing a structure, it would
be a good idea to go in and look the text of the pdb file. 

## 1: Necessary Unrecognized Residue
If there is an unrecognized residue,
  1. Go look into the structure and see what it is and what's going
  on. Some structures in the PDB may be C-alpha only like, 2rdo.pdb or
  3htc.pdb.

  2. If it is a HET group that you do not need, such as a non-covalent
  small molecule, like a NAD, HEME group or glycerol that is included
  with the structure, then you can include -ignore_unrecognized_res
  option to the command line.

  3. If there is a covalently linked ligand or modified amino acid or
  base, either there is a parameter file in the
  database/chemical/residue_type_sets and it is looked up
  correctly. If it is important, (for instance if it is part of the
  covalent chain like the CSP residue in 1qu9.pdb) but not in the
  database, it is possible to convert a MDL, MOLfile, or SDF file to a
  Rosetta residue topology file.
    ```
     src/python/apps/public/molfile_to_params.py --help 
    ```
  Please see documenation on preparing ligands for using
  molfile_to_params correctly.

## 2: Unnecessary Unrecognized Residue
If there are unrecognized residues but you don't actually care about them,

  1. Manually remove the HETATM atom records and checking in PyMOL
  that the structure is still reasonable (no broken chains etc.)

  2. Use a script to clean the PDB appropriately. For example
  tools/protein_util/scripts/clean_pdb.py can handle some basic
  cleanup proceedures. This script is a little hacky so use with
  caution.

  3. Run with `-ignore_unrecognized_res` or just `-ignore_water`. Simply
  ignoring unrecognized residues can coding assumptions that are not
  rigerously checked (i.e. expose bugs in Rosetta). So be prepared
  Rosetta to fail in unusual ways: 
     1. if the first residue of a chain is ignored (e.g. with 1mhm.pdb)
    then Rosetta may fail with this error:

        ```
        ERROR: ( anchor_rsd.is_polymer() && !anchor_rsd.is_upper_terminus() ) && ( new_rsd.is_polymer() && !new_rsd.is_lower_terminus() )
        ERROR:: Exit from: src/core/conformation/Conformation.cc line: 428
        ```

    2. if the residue type three letter code mis identified, such as
    SUC (sucrose in 3KB8.pdb) for the SUCK residue type, used in a
    score term to reduce buried cavities.

        ```
        core.io.pdb.file_data: [ WARNING ] discarding 23 atoms at position 207 in file /home/momeara/scr/data/pdb/kb/pdb3kb8.ent.gz. Best match rsd_type:  SUCK
        core.io.pdb.file_data: [ WARNING ] discarding 23 atoms at position 750 in file /home/momeara/scr/data/pdb/kb/pdb3kb8.ent.gz. Best match rsd_type:  SUCK
        core.io.pdb.file_data: [ WARNING ] discarding 23 atoms at position 980 in file /home/momeara/scr/data/pdb/kb/pdb3kb8.ent.gz. Best match rsd_type:  SUCK
        core.conformation.Conformation: [ WARNING ] missing heavyatom:  OXT on residue SER_p:CtermProteinFull 201
        core.conformation.Conformation: [ WARNING ] missing heavyatom: ORIG on residue SUCK 202
        ...
        ERROR: too many tries in fill_missing_atoms!
        ERROR:: Exit from: src/core/conformation/Conformation.cc line: 2590
        ```

## 3: If You Care About Hydrogens
If you care about the coordinates of hydrogens in the PDB file, you
should be careful that they are named correctly. Rosetta does not use
the most uptodate hydrogen atom naming convention. If the hydrogen
atom names are not recognized, then they will be stripped off and
replaced automatically by the -optH protocol when the pdb file is
processed.

To make sure H-atoms are named correctly you can use the
tools/convert_hatom_names.py script.

```
   python tools/convert_hatm_names.py --data_dir input_files --output_dir prepared_input_files
```

