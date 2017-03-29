# Relax Around A Chemically Bound Ligand

KEYWORDS: LIGANDS INTERFACES STRUCTURE_PREDICTION

This document was originally written by Andrew Leaver-Fay, and completed, expanded and modified for Doxygen by Ron Jacak. It was later modified during documentation XRW 2016 by Parisa Hosseinzadeh and Vikram K. Mulligan (vmullig@u.washington.edu) to enable automatic testing and to ensure compatibility. 

# Quick Guide
1. Create parameter files for the ligand, and modified residue (if one doesn't already exist)

2. Replace ligand HETATM lines with Rosetta ATOM lines, and change modified residue lines to proper type

3. Make constraints file

4. Run relax protocol

# Introduction

The purpose of this demo is to relax a protein with a chemically attached small molecule. Currently, one example on the LOV2 domain is provided. 

Rosetta requires a "params" file for each of the residue types that are used in a trajectory. These are already available for all of the amino acids, nucleic acids, and most of the metals.  These parameter files live in Rosetta/main/database/chemical/residue_type_sets/fa_standard/residue_types/. Before relax can be run on a liganded structure, a parameter file for that ligand needs to be created.

The "LOV2" domain binds a flavin and, when exposed to blue light, forms a thiol bond at CYS 450 to this flavin.  One question that would be nice to answer is what effect this chemical bond has on the stability of the LOV2 domain, since it is known that after this bond forms, the C-terminal "J-alpha" helix unfolds.  To address this with Rosetta, we would like to run relax on the starting structure while preserving this chemical bond.  The crystal structure 2V0W was made by crystalizing the LOV2 domain in the dark, and then exposing the crystals to light.  The thiol bond between CYS 450 and the flavin is visible in the structure. 

# Creating a Ligand Parameter File

We downloaded the .sdf file for the flavin attached to LOV2 from the PDB (http://www.pdb.org/pdb/explore/explore.do?structureId=2V0W), Ligands_noHydrogens_withMissing_1_Instances.sdf (an awful name for this file). Also, the file doesn't contain hydrogens, which we do want. So, first, we need to add hydrogens to this file.  The program REDUCE (which can be obtained from the [Richardson lab website](http://kinemage.biochem.duke.edu/software/reduce.php)) can be used to add hydrogens.  To run REDUCE, you have to put a copy (or symlink) of the reduce_het_dict.txt file in the /usr/local directory. You can also use other programs. For more information please refer to tutorials/prepare_ligand.

```bash
> cd /usr/local
> sudo ln -s /Users/andrew/reduce/reduce_het_dict.txt
> ~/reduce/reduce.3.03.061020.macosx.i386 2V0W.pdb > blah.pdb

> grep FMN blah.pdb | tail -n 9
HETATM    0 3HM8 FMN A1547      15.201   6.333   6.801  1.00 11.04           H   new
HETATM    0 3HM7 FMN A1547      13.095   4.009   5.709  1.00 11.55           H   new
HETATM    0 2HM8 FMN A1547      14.999   5.775   8.496  1.00 11.04           H   new
HETATM    0 2HM7 FMN A1547      12.825   3.407   7.380  1.00 11.55           H   new
HETATM    0 1HM8 FMN A1547      16.624   6.243   7.893  1.00 11.04           H   new
HETATM    0 1HM7 FMN A1547      13.426   5.078   7.114  1.00 11.55           H   new
HETATM    0  HN3 FMN A1547      19.774  -2.545   4.778  1.00 12.52           H   new
HETATM    0  H9  FMN A1547      18.030   4.577   7.483  1.00 12.40           H   new
HETATM    0  H6  FMN A1547      14.354   1.560   6.108  1.00 12.59           H   new
```

Now open up the PDB file of the ligand in pymol and have PyMOL write it out as a .mol file.

```bash
> grep FMN blah.pdb > FMN_w_h.pdb
> pymol FMN_w_h.pdb
```
--> save molecule

--> save as one file

--> file format: MOL

produces FMN_w_h.mol

We can convert the MOL file to a .params file with a script named molfile_to_params.py (located in <path_to_Rosetta_directory>/main/source/scripts/python/public/). This script can be run with the "-h" flag to list all the flags that are applicable.

```bash
$> cp rosetta_inputs/FMN_w_h.mol .
$> $ROSETTA3/scripts/python/public/molfile_to_params.py -n FMN -p FMN FMN_w_h.mol
```

(where `$ROSETTA3`=path-to-Rosetta/main/source)

The -n flag specifies what 3-letter name will used to identify this ligand in PDB files. The -p flags specifies a name for the .params file. We chose "FMN" which stands for flavin mononucleotide.  The python script has the possibility of mis-assigning atom types, so we opened up the FNM.params file to look and see that we got the chemistry right.

Example output from a molfile_to_params.py run:
```
Centering ligands at (  19.791,    2.871,    7.017)
Atom names contain duplications -- renaming all atoms.
WARNING:  atom  P1  has valence > 4
WARNING:  structure contains double bonds but no aromatic bonds
  Aromatic bonds must be identified explicitly --
  alternating single/double bonds (Kekule structure) won't cut it.
  This warning does not apply to you if your molecule really isn't aromatic.
Total naive charge -7.250, desired charge 0.000, offsetting all atoms by 0.234
WARNING: fragment 1 has 31 total atoms including H; protein residues have 7 - 24 (DNA: 33)
WARNING: fragment 1 has 31 non-H atoms; protein residues have 4 - 14 (DNA: 22)
WARNING: fragment 1 has 7 rotatable bonds; protein residues have 0 - 4
Average 31.0 atoms (31.0 non-H atoms) per fragment
(Proteins average 15.5 atoms (7.8 non-H atoms) per residue)
WARNING:  no root atom specified, using auto-selected NBR atom instead.
Wrote params file FNM.params
Wrote PDB file FNM_0001.pdb
```

The .params file looks good, and the PDB file FMN_0001.pdb looks good. The other thing we have to do is to tell Rosetta that there is a chemical bond between this flavin and some other residue.  A chemical connection must be placed between the flavin atom C8 and SG atom of the modified CYS
residue. To do this, add a "CONNECT C8" line to the FMN params file (anywhere after the ATOM records is fine). 

Finally, we need to add an "ICOOR" line for this connection at the bottom of the file (and to get this, we'll measure the distance to the CYS SG in the 2V0W pdb). The format for this line is:

```
ICOOR_INTERNAL CONN1 <torsion> <180 - N3:C8:SG angle> <C8:SG distance> <parent atom name> <grandparent name> <great-grandparent name>
```

The dihedral isn't particularly useful, but we do need to know the C8-SG distance, C8-SG-CB angle, and the last angle from the "parent atom" of C8 to C8 to SG. we can find the parent atom in the "ICOOR" line for C8 in the .params file.

```
ICOOR_INTERNAL C8 27.068913 64.302080 1.415450 N3 C9 C16
```

This says "N3" is the parent, C9 is the grandparent, and C16 is the great-grandparent.  So we need the N3-C8-SG angle.

```
C8-SG: 1.9 A
C8-SG-CB: 111.2 degrees
N3-C8-SG: 110.7 degrees
```

We'll need some of these parameters for the ICOOR, and some of them later for the constraints we'll use to enforce good geometry for this chemical bond (Rosetta does not do that for you!).

Thus, the line to be added should be the following:<br>

```
ICOOR_INTERNAL  CONN1  180.000000   69.300000    1.900000   C8    N3    C9
```

The name "CONN1" says that this is the location for the atom that is covalently connected to this residue. The number is taken as the order in which the connection was listed in the .params file, and in this case, there is only one CONNECT connection. This modified .params file is FMN_modded.params.

# Using the Newly Created Ligand and Modified Residue Parameter Files

The next step is to swap the newly generated FMN with the existing FMN residue: the point is to use the newly generated names for the atoms.  Atom names have to agree (perfectly!) between the input PDB file and the .params file.  Careful: the .params file has a strange format. 
The atom names given on the ATOM lines are column formatted, so don't add extra whitespace. The rest of the line (not the names) is whitespace delimited.  Paste the contents of the FMN_0001.pdb file to the bottom of the original 2V0W file, and remove the HETATM lines using the original flavin molecule atom names. To add the newly created ligand residue to the database of residue types Rosetta uses, use the flag -extra_res_fa on the command line:

```
-extra_res_fa FMN_modded.params
```

# Making a constraints file

Now we need to generate a set of constraints to make Rosetta preserve the bond geometry between the CYS and C8.  Constraint files rely on Rosetta numbering which starts at 1 and does not use the PDB numbering. The first residue in this PDB is 401, so CYS 450 will be Rosetta-residue 50. The last residue in the protein is 546, so the ligand (i.e. last) will be Rosetta-residue 147. The distance constraint is given by this line:
```
AtomPair SG 50 C8 147 HARMONIC 1.9 0.01
```

where HARMONIC is the functional form, 1.9 is the x0 (the ideal value, in Angstroms), and 0.01 is the standard deviation, i.e. the inverse of the spring constant. The two angle constraints are given by these lines:
```
Angle SG 50 C8 147 N3 147 HARMONIC 1.93207948 0.034906585
Angle CB 50 SG 50 C8 147 HARMONIC 1.94080613 0.034906585
```

where for the first constraint, 1.93207948 represents the ideal angle (x0) of 110.7 degrees, given in radians, and 0.034906585 represents a standard deviation of 2 degrees, given in radians.

```
From Google:
110.7 * pi / 180 = 1.93207948
2.0 * pi / 180 = 0.034906585
```

For more information about the constraints, you can check the [constraints tutorial](../../tutorials/Tutorial_8_Constraints/Tutorial_8_Constraints.md).

# Running the relax protocol

To run the relax protocol, we need to pass in a PDB with the correct FMN lines, the parameter files for the modified CYS and the FMN, and the constraint file.  We also need to activate constraints during scorefunction evaluation, which can be done using the score:weights flag on the command line.  
You can copy the necessary files that are pre-generated from rosetta_inputs directory:

```bash
$> cp rosetta_inputs/2V0W.pdb .
$> cp rosetta_inputs/FMN_modded.params .
$> cp rosetta_inputs/chemical_bond.cst .
```

A complete command line for the protocol would be as follows:

```bash
$> $ROSETTA3/bin/relax.default.linuxgccrelease -s 2V0W.pdb -in:file:fullatom -extra_res_fa rosetta_inputs/FMN_modded.params -overwrite -cst_fa_file chemical_bond.cst -score:weights talaris2014_cst.wts -mute basic core.init core.scoring
```

Truncated output from the command:
```
protocols.relax.FastRelax: ================== Using default script ==================
protocols.jd2.PDBJobInputter: Instantiate PDBJobInputter
protocols.jd2.PDBJobInputter: PDBJobInputter::fill_jobs
protocols.jd2.PDBJobInputter: pushing 2V0W.pdb nstruct index 1
protocols.jd2.PDBJobInputter: PDBJobInputter::pose_from_job
core.chemical.ResidueTypeSet: Finished initializing fa_standard residue type set.  Created 6480 residue types
core.conformation.Conformation: Connecting residues: 50 ( CYS ) and 147 ( FMN ) at atoms  SG  and  C8 
core.conformation.Conformation:  with mututal distances: 2.63043 and 2.09338
core.import_pose.import_pose: Can't find a chemical connection for residue 147 FMN
core.pack.task: Packer task: initialize from command line() 
protocols.jd2.PDBJobInputter: filling pose from PDB 2V0W.pdb
core.io.constraints: read constraints from chemical_bond.cst
HARMONIC 1.93208 0.0349066

HARMONIC 1.94081 0.0349066

core.pack.task: Packer task: initialize from command line() 
core.pack.dunbrack: Dunbrack library took 0.028403 seconds to load from binary
protocols.relax.FastRelax: CMD: repeat  -80.6309  0  0  0.44
core.pack.interaction_graph.interaction_graph_factory: Instantiating DensePDInteractionGraph
core.pack.pack_rotamers: built 1514 rotamers at 147 positions.
core.pack.pack_rotamers: IG: 1208136 bytes
core.pack.pack_rotamers: pack_rotamers_run(): simulated annealing took 0.080246 seconds
```

The relax protocol should output a structure named 2V0W_0001.pdb which has lower energy than the starting structure.

