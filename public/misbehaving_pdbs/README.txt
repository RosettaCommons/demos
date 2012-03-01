
For the most part Rosetta does a pretty good job at reading structures, however with unusual chemical modifications, sometimes Rosetta runs into problems. Here are some of the problems with pdb files that we have encountered that currently cause problems with Rosetta:

For several of the structures that gave us problems, we tried to
isolate the problems were not able to find a consistent issue.



## 0 ##
If you run into problems with processing a structure, it would
be a good idea to go in and look the text of the pdb file. 

## 1 ##
If there is an unrecognized residue,
  1) Go look into the structure and see what it is

  2) If it is a HET group that you do not need, such as a non-covalent
  small molecule, like a NAD, HEME group or glycerol that is included
  with the structure, then you can include -unrecognized_res option to
  the command line.

  3) If there is a covalently linked ligand or modified amino acid or
  base, either there is a parameter file in the
  minirosetta_database/chemical/residue_type_sets and it is looked up
  correctly. If it is important, (for instance if it is part of the
  covalent chain like the CSP residue in 1qu9.pdb) but not in the
  database, it is possible to convert a MDL, MOLfile, or SDF file to a
  Rosetta residue topology file.

     src/python/apps/public/molfile_to_params.py --help 

## 2 ##
A bug that should be fixed is that if there is an unrecognized residue
at the beginning of a chain, for example 1mhm.pdb, it leads to an
assert failure like this:

  ERROR: ( anchor_rsd.is_polymer() && !anchor_rsd.is_upper_terminus() ) && ( new_rsd.is_polymer() && !new_rsd.is_lower_terminus() )
  ERROR:: Exit from: src/core/conformation/Conformation.cc line: 428

If this is a modified amino acid or base, one way to "hack" is to change the amino acid to its common/unmodified form.


## 3 ##
Structures with C-alpha only chains such as, 2rdo.pdb or 3htc.pdb, causes an error message like this.
ERROR:: Exit from: src/protocols/jobdist/standard_mains.cc line: 602
 

## 4 ##
There is a bug in Rosetta that sometimes when it is matching residue types it gets into infinite loops.  For example in 3kb8.pdb:

  core.io.pdb.file_data: [ WARNING ] discarding 23 atoms at position 207 in file /home/momeara/scr/data/pdb/kb/pdb3kb8.ent.gz. Best match rsd_type:  SUCK
  core.io.pdb.file_data: [ WARNING ] discarding 23 atoms at position 750 in file /home/momeara/scr/data/pdb/kb/pdb3kb8.ent.gz. Best match rsd_type:  SUCK
  core.io.pdb.file_data: [ WARNING ] discarding 23 atoms at position 980 in file /home/momeara/scr/data/pdb/kb/pdb3kb8.ent.gz. Best match rsd_type:  SUCK
  core.conformation.Conformation: [ WARNING ] missing heavyatom:  OXT on residue SER_p:CtermProteinFull 201
  core.conformation.Conformation: [ WARNING ] missing heavyatom: ORIG on residue SUCK 202
  .... < many lines like this > .....
  
  ERROR: too many tries in fill_missing_atoms!
  ERROR:: Exit from: src/core/conformation/Conformation.cc line: 2590

Another example is starting_structures/1ql0AH. Note the pdb in the protein databank 1ql0.pdb is pocessed correctly.
  

## 5 ##
We recommend Rosetta community:
 
  1) adding an n-terminal ACE as to the available residue params
files.  For example see 1bbzEFH in the starting_inputs folder.

  2) add an -ignore_waters flag

  3) Finding pdb preparation scripts in the rosetta labs, collect them
  either put them in one place or incoperate the behavior into how
  Rosetta process structures. Be sure to ask the Richardsons because
  they have a pdb preparation script for MolProbity.

