# Supercharge: Reengineer proteins for high net surface charge

KEYWORDS: DESIGN GENERAL

If you want to run supercharge now, the application is called 'supercharge' in `src/apps/public/supercharge.cc`.

Here are four examples:

(where `ROSETTA3`=path-to-Rosetta/main/source)

```
$> $ROSETTA3/bin/supercharge.default.macosgccrelease @rosetta_inputs/options1 
```
// Rosetta-mode, positive-charge, fixed surface cutoff and input ref energies
```
$> $ROSETTA3/bin/supercharge.default.macosgccrelease @rosetta_inputs/options2
```
// Rosetta-mode, negative-charge, fixed surface cutoff and target net charge
```
$> $ROSETTA3/bin/supercharge.default.macosgccrelease @rosetta_inputs/options3 
```
// AvNAPSA-mode, negative-charge, target net charge
```
$> $ROSETTA3/bin/supercharge.default.macosgccrelease @rosetta_inputs/options4 
```
// AvNAPSA-mode, positive-charge, fixed surface cutoff



Rosetta-mode and AvNAPSA-mode are explained below...

# Why supercharge protein surfaces?

Reengineering protein surfaces to have high net charge, called supercharging, can improve reversibility of unfolding by preventing aggregation of partially unfolded states. Aggregation is a common obstacle for use of proteins in biotechnology and medicine.  Additionally, highly cationic proteins and peptides are capable of nonviral cell entry, and highly anionic proteins are filtered by kidneys more slowly than neutral or cationic proteins.  

Optimal positions for incorporation of charged side chains should be determined, as numerous mutations and accumulation of like-charges can also destabilize the native state.  A previously demonstrated approach deterministically mutates flexible polar residues (amino acids DERKNQ) with the fewest average neighboring atoms per side chain atom (AvNAPSA: Lawrence MS, Phillips KJ, Liu DR, 2007, Supercharging proteins can impart unusual resilience, JACS).  Our approach uses Rosetta-based energy calculations to choose the surface mutations.  Both automated approaches for supercharging are implemented in this online server.

# Two Approaches
There are two automated approaches, **Rosetta supercharge (Rsc)** and **AvNAPSA supercharge (Asc)**

**AvNAPSA supercharge philosophy (Asc):** mutate the most exposed polar residues to minimize structural change or destabilization.  Only DE-RK-NQ residues can be mutated.

**Rosetta supercharge philosophy (Rsc):** mutate residue positions that preserve and/or add favorable surface interactions.  Hydrophobic and small polar surface residues can also be mutated.

**AvNAPSA drawbacks:** mutating surface polar residues can eliminate hydrogen bonds.  Helix capping, edge-strand interaction, and loop stabilization all result from surface hydrogen bonds.  Furthermore, this automated protocol mutates N to D and Q to E, but N and Q sometimes act simultaneously as a donor and acceptor for hydrogen bonds.

**Rosetta drawbacks:** mutating less-exposed positions can lead to better computed energies, but mistakes at these positions can be destabilizing.  AvNAPSA favors charge swaps, so Rosetta requires more mutations to accomplish the same net charge.

The AvNAPSA approach varies net charge by adjusting the surface cutoff.  The Rosetta approach varies net charge by adjusting reference energies of the positive or negatively charged residues.

The supercharge server can run in four different modes:
-AvNAPSA with a target net charge
-AvNAPSA with a surface cutoff
-Rosetta with a surface cutoff and target net charge
-Rosetta with a surface cutoff and input reference energies for charged residue types


**What does AvNAPSA stand for:** average number of neighboring atoms per sidechain atom.  This is a value that measures the extent of burial/accessibility.  It's similar to the residue neighbors by distance that Rosetta typically uses to define the surface, but it's on the atom-level rather than residue-level.  AvNAPSA-mode calculates an AvNAPSA value for every residue.  'surface_atom_cutoff' indicates the cutoff AvNAPSA value that defines surface residues.  AvNAPSA values of 50-150 are typical for surface residues.  AvNAPSA values >150 are typical for core residues.  A surface_atom_cutoff of 100 will lead to moderate supercharging.  A surface_atom_cutoff of 150 will lead to heavier supercharging.


## Workflow of Each Mode

### AvNAPSA-mode, target charge
1. Define surface.  sort NQ and RK/DE residues by AvNAPSA value (low to high)
2. Next residue in sorted list: Positive: mutate DENQ-->K, Negative: mutate RKQ-->E and N-->D
3. If net charge = target net charge, output pdb

### AvNAPSA-mode, surface cutoff
1. Define surface by AvNAPSA value (<100 default)
2. For each NQ and DE/RK residue in the surface: Positive: mutate DENQ-->K, Negative: mutate RKQ-->E and N-->D
3. Output pdb

### Rosetta-mode, surface cutoff and target charge
1. Define surface.  Neighbor by distance calculator (CB dist.), <16 neighbors default or Define surface by AvNAPSA value (<100 default)
2. Set design task
    - read user resfile, if provided
    - dont_mutate gly, pro, cys
    - dont_mutate h-bonded sidechains
    - dont_mutate correct charge residues
3. Set reference energies for RK/DE, starting at user input values
4. pack rotamers mover
5. check net charge, increment/decrement reference energies (back to step 3.)
6. Once a pack rotamers run results in the correct net charge, output pdb

### Rosetta-mode, surface cutoff and input reference energies for charged residue types
1. Define surface.  Neighbor by distance calculator (CB dist.), <16 neighbors default or Define surface by AvNAPSA value (<100 default)
2. Set design task
    - read user resfile, if provided
    - dont_mutate gly, pro, cys
    - dont_mutate h-bonded sidechains
    - dont_mutate correct charge residues
3. Set reference energies for RK/DE, using the user input values
4. pack rotamers mover
5. Output pdb


# Options

## AvNAPSA Mode

```
AvNAPSA_positive  BOOL def(false);				//run positive-charge AvNAPSA
AvNAPSA_negative  BOOL def(false); 				//run negative-charge AvNAPSA
target_net_charge  SIGNED_INT def(0);  				//residue positions will be mutated one at a time from most exposed to least exposed until target net charge is achieved
surface_atom_cutoff  UNSIGNED_INT def(100); 			// if you have no target net charge in mind, AvNAPSA will mutate all surface DE-RK-NQ residues on the surface, with this surface cutoff
```

## Rosetta Mode

```
surface_residue_cutoff  UNSIGNED_INT def(16);  //residues with <16 neighboring residues within 10 Å are considered part of the surface

include_arg  BOOL def(false);  //use arginine in Rosetta supercharge
include_lys  BOOL def(false);  //use lysine in Rosetta supercharge
include_asp  BOOL def(false);  //use aspartate in Rosetta supercharge
include_glu  BOOL def(false);  //use glutamate in Rosetta supercharge

//the reference energies of the charged residue types will govern the net charge of Rosetta designs.  Rosetta can choose between the allowed charged residue types and the native residue.  More negative reference energies will result in more charge mutations.
refweight_arg  FLOAT def(-0.98);
refweight_lys  FLOAT def(-0.65);
refweight_asp  FLOAT def(-0.67);
refweight_glu  FLOAT def(-0.81);

dont_mutate_glyprocys  BOOL def(true);		 //glycine, proline, and cysteine often serve special structural roles in proteins
dont_mutate_correct_charge  BOOL def(true);      //i.e., don’t mutate arginine to lysine
dont_mutate_hbonded_sidechains  BOOL def(true);  //don’t mutate residues with sidechains forming a hydrogen bond
pre_packminpack  BOOL def(false);                //Packrotamers is always done as the first step.  This option will go one step further and run packrotamers, sidechain+backbone minimization, packrotamers on the input structure before performing the supercharge design step.

nstruct  UNSIGNED_INT def(1);  			 //Monte Carlo sequence design of a protein surface is often convergent but it is still stochastic, multiple design runs can be performed if desired.
target_net_charge  UNSIGNED_INT def(0);  	 //a target net charge can be achieved if desired, this is done in an automated way by incrementing/decrementing charged residue reference energies until the desired net charge results from the Monte Carlo design step.
```

## AvNAPSA and Rosetta Mode

```
surface_atom_cutoff  UNSIGNED_INT def(100); // this is how AvNAPSA defines surface, can be used in either approach
compare_energies  BOOL def(false);  	 	      		//prints a full residue-by-residue energy analysis in the log file
only_compare_mutated_residues  BOOL def(false);  		//only includes mutated residues in the energy analysis
resfile  FILE;  	       	    				//this is how you can specify which residues to not mutate.  Default setting must be ALLAA, and residue-by-residue settings should be NATAA, as shown below:

ALLAA
start
  20  A  NATAA
  24  A  NATAA
  26  A  NATAA


Note: an input resfile is optional.  However, every supercharge run generates an output resfile that governs the design run.  The default of this output resfile is NATAA, which prevents core residues from mutating (see below).  The input resfile is read first, the output resfile (see below) is read second, and this is why ALLAA must be the default for the input resfile.  If the default were NATRO, for example, no design would occur!
```

# Output

As output, a log file, the residue file that governed the design run, and the output PDB are provided.  First, the log file contains the exact Rosetta command line, the residue positions identified as located on the surface, a list of charged residues in the final sequence, the net charge, a list of mutations, text for a PyMOL selection to easily view the mutations in PyMOL, and optionally, a full energetic comparison of repacked native versus supercharged structures.  Secondly, the Rosetta residue file indicates which residue positions could possibly mutate, and to what residue types.  The third output file is the atomic coordinate file of the supercharged protein, in PDB format, and the naming of the output PDB is intended to facilitate self-documentation of the inputs for a given design run.  For Rosetta designs, the name includes the final reference energies that were used and the final net charge, and for AvNAPSA designs, the name includes the net charge and the largest AvNAPSA value of the mutated residues. 


This is what an output resfile looks like for AvNAPSA-positive supercharging, which always chooses lysine:

```
NATAA
start
   6 A  PIKAA  K
   19 A  PIKAA  K
   21 A  PIKAA  K
   32 A  PIKAA  K
   34 A  PIKAA  K
   39 A  PIKAA  K
```

This is what an output resfile looks like for Rosetta positive supercharging, which allows choice between native and RK, and preserves h-bonds:

```
NATAA
start
   6 A  PIKAA ERK
   9 A  PIKAA TRK
   11 A  PIKAA VRK
   21 A  PIKAA DRK
   25 A  PIKAA HRK
   26 A  NATAA  #same charge
   30 A  NATAA  #same charge
   32 A  NATRO  #has sc hbond energy=-1.15844
   38 A  PIKAA TRK
   39 A  NATRO  #has sc hbond energy=-1.33149
   43 A  PIKAA TRK
   50 A  NATRO  #has sc hbond energy=-0.536622
   52 A  NATAA  #same charge
   76 A  PIKAA DRK
   77 A  PIKAA HRK
```

