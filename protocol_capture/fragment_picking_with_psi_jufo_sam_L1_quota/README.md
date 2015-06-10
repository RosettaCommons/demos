Fragment Picking with Psi Jufo SAM L1 Quota
===========================================

Presenting author: Dominik Gront (dgront at chem dot uw dot edu dot pl)
Protocol name: fragment picker : CS-Rosetta style
Brief description: The protocol substitutes nnmake

Source code location
--------------------

* Check out the mini SVN: https://svn.rosettacommons.org/source/trunk/mini/
* Fragment picker is located in: https://svn.rosettacommons.org/source/trunk/mini/src/core/fragment/picking
* Applications are in: https://svn.rosettacommons.org/source/trunk/mini/src/apps/pilot/dgront/fragmentpicker

Running the protocol capture
----------------------------

1. Set up the path to minirosetta database
2. Set up the path to vall database
3. Run the picker:
   ```
   picker.linuxgccrelease @2jsvX-ss_sim_l1-Q.flags
   ```
