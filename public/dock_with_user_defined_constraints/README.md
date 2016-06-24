Docking with constraints
========================

KEYWORDS: DOCKING INTERFACES

In this demo, a peptide (chain B) is docked onto a protein (chain A) with 
user-defined constraints. The purpose of this demo is to illustrate the use of 
constraint files and as many different types of constraints as possible. The 
input files are based on PDB ID (2o2m) which is a "Crystal Structure Of Human 
Gαi1 Bound To The Goloco Motif Of Rgs14". The protein is chain A. The peptide 
is chain B. 

Protocol
--------

A general overview of the protocol is available at:  
http://www.rosettacommons.org/manuals/archive/rosetta3.2.1_user_guide/app_docking.html

A write-up of the file format for Constraints is available at:  
http://www.rosettacommons.org/manuals/archive/rosetta3.2.1_user_guide/constraint_file.html

High level overview:

1. Docking Prepack
    - Useful to optimize all of the sidechains on both partners
    - http://graylab.jhu.edu/Rosetta.Developer.Documentation/all_else/de/d69/docking_prepack_protocol.html

2. Docking with constraints. In general, constraints are typically defined for CA atoms of the partner sidechains.
    - Docking will honor AtomPairContraint, AmbiguousConstraint, SiteContraint. 
      For additional information on contraints, see "Constraints" section 
      below. 

3. Analysis
    - The expected output of the docking run is (1) a decoy with both docked 
      partners and (2) a score file with score breakdown for each decoy 
    - Typical analysis would be to sort the generated score file with a command 
      line i.e. 'sort -nk 28` where the integer identifies the column number 
      for the quality to sort on.

Constraints
-----------

A good overview of constraints is available in the documentation:  
http://www.rosettacommons.org/manuals/archive/rosetta3.2.1_user_guide/constraint_file.html

The three constraints which are honored by the Docking Protocol are: 

* AtomPairConstraint
* AmbiguousConstraint
* SiteConstraint.

AtomPair and AmbiguousConstraint are described the above documentation for 
Rosetta 3.2.1.  A SiteConstraint allows you to specify that a particular 
residue should be in contact with a particular chain. An example of a 
SiteConstraint is:

    SiteConstraint CA 4A D FLAT_HARMONIC 0 1 5

This will add a FLAT_HARMONIC potential with the parameters 0 1 5 (recommended) 
around the distance between the CA of residue 4 (PDB numbering) on chain A and 
the closest CA on chain D to the ScoreFunction. 

To run these examples, use
```bash
$> <path/to/Rosetta/main/source/bin/>/docking_protocol.default.linuxgccrelease @rosetta_inputs/options
$> <path/to/Rosetta/main/source/bin/>/docking_protocol.default.linuxgccrelease @rosetta_inputs/options.constraint_ambiguous
$> <path/to/Rosetta/main/source/bin/>/docking_protocol.default.linuxgccrelease @rosetta_inputs/options_wo_constraints
```

