# Symmetric Docking of Insulin Trimer of Dimers - Demo

KEYWORDS: DOCKING SYMMETRY

## Needed input files

Input files:

1. one pdb file
2. one symmetry definition file

(1) The pdb file represents the monomer entity; regardless if this file has multiple physical chains or its own internal symmetry this, is the monomer. However, for safe running of this script the pdbfile must be rewritten to have a single chain label.

(2) The symmetry definition file informs Rosetta about the symmetry type, number of monomers to dock, how to calculate the score and some more things.

## Example 1: Docking of an Insulin Dimer (chains A, B, C, D in 1ZEH) into a Trimer of Dimers (C3 symmetry)

### Clean Up the PDB File

IF the pdb file is multi-chain we must normalize it to a pseudo single chain with the following convenience script "addChain.pl". You can find that script in the scripts directory of this tutorial. This step is only required when the pdb is multi-chain!

```
	$> ./scripts/addChain.pl 1ZEH.pdb A | grep '^ATOM'  >  rosetta_inputs/1ZEH_monomer_c3.pdb
```

where 1ZEH.pdb is the pdb file input and A is the the dummy new chain ID.
This simple script only normalizes the ATOM records; it removes the TER chain separators and deletes unnecessary SEQRES, REMARK, and other lines.

In the example here, 1ZEH.pdb, contains 4 chains, and the 1ZEH_monomer.pdb has a single chain joining all the ATOM records of these into one chain called A. Now we have the pdb input consisting of only our chosen monomer entity.

### Generating the Symmetry File

To generate a de novo symmetry file use script: `make_symmdef_file_denovo.py` that comes with the Rosetta code.
The minimal inputs of this script are (1) the symmetry: cn (circular symmetry) or dn (dihedral symmetry) and (2) the number of "monomers" (i.e. subunits).

```
	$ python $ROSETTA3/src/apps/public/symmetry/make_symmdef_file_denovo.py -symm_type cn -nsub 3  > rosetta_inputs/1ZEH.c3.symm
```

where "cn" is cyclic symmetry and "3" is the number of "monomers" (that gives a 3 fold single axis rotational symmetry). This produces the symmetry definition file for Rosetta. The anchor point in the symmetry def is set to center of mass (COM). This can cause problems. If it SymDock fails, this can normally be replaced by a residue numer.

**Note:** Pay no attention to the fact that 1ZEH has it's own internal symmetries. It is just considered a monomer subunit.

The output file 1ZEH.c3.symm contains the following:

```
	symmetry_name c3
	subunits 3
	recenter
	number_of_interfaces  1
	E = 3*VRT0001 + 3*(VRT0001:VRT0002)
	anchor_residue COM
	virtual_transforms_start
	start -1,0,0 0,1,0 0,0,0
	rot Rz 3
	virtual_transforms_stop
	connect_virtual JUMP1 VRT0001 VRT0002
	connect_virtual JUMP2 VRT0002 VRT0003
	set_dof BASEJUMP x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)
```

### Run Rosetta

```
	$> $ROSETTA3/bin/SymDock.default.macosclangrelease @flags.c3
```

where the provided flags file (flags.c3) looks like this:

```
	-database ~/Rosetta/main/database   # change path to your directory!!!
 	-in:file:s 1ZEH_monomer.pdb
 	-symmetry:symmetry_definition 1ZEH.c3.symm
 	-packing:ex1
 	-packing:ex2aro
 	-out:nstruct 1   # typically this should be 10,000 or more.  reduced to 1 for demo.
 	-out:file:fullatom
 	-symmetry:initialize_rigid_body_dofs
 	#-symmetry:symmetric_rmsd     # this line is only used when a reference pdb at the full C3 symmetry is available for testing purposes.
 	-ignore_unrecognized_res
```

### Results

The Rosetta run generated 1ZEH_monomer_c3_0001.pdb and score.fasc. The pdb file has the trimer decoy of the monomer docked in C3 symmetry.  The 3 monomers are entered as chains A,B and C in the trimer. The score.fasc is a standard score file as one would get for any normal docking run. The score is the score of the complete trimer. (See the PLOS one paper "Modeling Symmetric Macromolecular Structures in Rosetta3" or the symmetrical docking documentation to understand how this is computed. Conceptually you may think of it as the internal score of the monomer plus the score of it's interface; all multiplied by 3.)

Now for a nominally trickier example.

## Example 2: Docking of Six Monomers (chains A, B in 1ZEH) into a Trime of Dimers (D3 symmetry)

### Clean Up the PDB File

If you had taken a peek at the file 1ZEH.pdb beforehand you might have noticed that it already had an internal dimer symmetry.  The four chains were arranged so that chain A and B were in C2 symmetry to C and D.  So the above example just computed a trimer of dimers.  The program did not know anything about the c2 symmetry because it just considered the file 1ZEH_monomer.pdb to be the subunit for the trimer.

Now we are going to run this again but first removing chains C and D from the 1ZEH.pdb file.

```
	$> perl -wane 'print if m/^ATOM/ and ($F[4] eq "A" or $F[4] eq "B")'  1ZEH.pdb  > 1ZEH_monomer_d3.pdb 
```

>1ZEH_monomer_d3.pdb can be removed after the next step!  

And re-chain it to a single pseudo chain.

```
	$> ./scripts/addChain.pl 1ZEH_monomer_d3.pdb A > rosetta_inputs/1ZEH_monomer_d3.pdb
```

### Generating the Symmetry File

Then we will use d3 symmetry for docking.  As you know, d3 symmetry is c3 symmetry with a mirror plane. In other words this has a total of 6 subunits.

```
	$> $ROSETTA3/src/apps/public/symmetry/make_symmdef_file_denovo.py -symm_type dn -nsub 6  > rosetta_inputs/1ZEH.d3.symm
```

**Note:** the number of subunits in d3 symmetry is 6.

Here is the symmetry definition file 1ZEH.d3.symm

```
symmetry_name d3
subunits 6
recenter
number_of_interfaces  4
E = 6*VRT0001 + 6*(VRT0001:VRT0002) + 3*(VRT0001:VRT0004) + 3*(VRT0001:VRT0005) + 3*(VRT0001:VRT0006)
anchor_residue COM
virtual_transforms_start
start -1,0,0 0,1,0 0,0,0
rot Rz_angle 120.0
rot Rz_angle 120.0
rot Rx_angle 180.0
rot Rz_angle 120.0
rot Rz_angle 120.0
virtual_transforms_stop
connect_virtual JUMP1 VRT0001 VRT0002
connect_virtual JUMP2 VRT0002 VRT0003
connect_virtual JUMP3 VRT0003 VRT0004
connect_virtual JUMP4 VRT0004 VRT0005
connect_virtual JUMP5 VRT0005 VRT0006
set_dof BASEJUMP x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)
set_dof JUMP3 z(50) angle_z(0:60.0)
```

### Run Rosetta

```
	$> $ROSETTA3/bin/SymDock.macosclangrelease @flags.d3
```

The result is again a pdbfile 1ZEH_haptomer_00001.pdb which has 6 monomers in d3 symmetry.   Note, again, the monomer in this second example had half as many atoms as the first example but since we used 6 of these instead of 3 the final structure has the same number of total atoms.  The score is appended into score.fasc

## Remarks

Not shown in this example is the optional but recommended practice of "pre-packing".  The docking protocol packs the rotamers of all residues in the interfaces.  However it does not pack any residues that are not in the interface.  Therefore a recommended practice is the pack the monomer ahead of time.  This is simply running the packer on the monomer file with any desired packing options.  Then follow the above protocol as written to the symmetric docking.  A possible point of confusion here is that there is a flag for docking named "prepack" however this does not pack the interior as desired here.  Instead use the Docking_prepack application (see below).


There are other flags associated with both "docking" and "symmetry" available beyond the ones show in this tutorial.  These flags are documented in their respective formal documentation section.

## Links

- https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/symmetry
- https://www.rosettacommons.org/docs/latest/application_documentation/utilities/make-symmdef-file
- https://www.rosettacommons.org/docs/latest/application_documentation/docking/sym-dock
- https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/SymDofMover

# Appendix

## Summary of the Commands

```
./scripts/addChain.pl 1ZEH.pdb A | grep '^ATOM'  >  1ZEH_monomer.pdb
src/apps/public/symmetry/make_symmdef_file_denovo.py -symm_type cn -nsub 3  > 1ZEH.c3.symm
bin/SymDock.macosgccrelease @flags.c3  >& 1ZEH.c3.log
perl -wane 'print if m/^ATOM/ and ($F[4] eq "A" or $F[4] eq "B")'  1ZEH.pdb  > /tmp/1ZEH_haptomer.pdb
./scripts/addChain.pl /tmp/1ZEH_haptomer.pdb A > 1ZEH_haptomer.pdb
src/apps/public/symmetry/make_symmdef_file_denovo.py -symm_type dn -nsub 6  > 1ZEH.d3.symm
bin/SymDock.macosgccrelease @flags.d3  >& 1ZEH.d3.log
```

## Authors
- Sebastian Rämisch
- Charlie E. M. Strauss
- Bruno Correia

**Rosetta Revision: 43611**

**Date: 08/05/2011**
