HBS Design Demo
===============

KEYWORDS: LIGANDS DOCKING

Written by Kevin Drew (Bonneau Lab), kdrew at nyu dot edu

This demo shows how to run the hbs_design application.  A hydrogen bond 
surrogate (HBS) is a helical mimetic scaffold used for inhibiting protein 
interactions. The demo shows the design of an hbs inhibitor for the MDM2-P53 
protein interaction.

Algorithm
---------

1. Pertubation phase: rigid body movement of hbs wrt target
2. Design phase: design user specified residues on hbs scaffold and minimize
3. Repeat 10x

Command
-------

    $> <path/to/Rosetta/>main/source/bin/hbs_design.default.linuxgccrelease @inputs/flags

Input Files
-----------

* `./input/flags`: User specified options.
* `./input/mdm2_hbs.pdb`: Input structure where target is chain 1 and hbs is 
  chain 2

Options
-------

* `-hbs_design_positions`: Residues on hbs to design (numbering is relative to 
  hbs, for example 3 is the third residue on hbs), default repacks with no 
  design.
* `-pert_num`: Number of pertubations during pertubation phase, default 10, 
  production 100.
* `-design_loop_num`: Number of pertubation + design cycles, default 10, 
  production 10.
* `-nstruct`: For production runs, use 1000.

Pre-processing
--------------

...

Post-processing
---------------

Similar to other multi chain design protocols, the ddG is computed and is a 
good indicator of a good design.  First sort by total score, take top 5 percent 
and then sort by REPACK_ENERGY_DIFF (ddG).

Limitations
-----------

This app is inflexible to adjusting Monte Carlo temperatures, score functions, 
degree of rigid body pertubations, designing noncanonical amino acids, etc. The 
app also requires the hbs is close to a plausible binding mode with respect to 
the target.

