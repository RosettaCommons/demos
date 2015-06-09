# Zinc Site Redesign

This is a demo for the mononuclear zinc metalloenzyme redesign procedure described in Khare, Kipnis, Greisen et al. Nature Chemical Biology (2012).

Authors: Sagar Khare (khares@uw.edu), Per Jr Greisen (pgreisen@gmail.com)

## Summary

Starting from a TS ensemble model and a set of zinc-containing PDBs (only one - PDBid 1A4L - is included here, list of other PDB codes used in the paper is in list\_of\_input\_pdbs file) as inputs, we generate a design model and evaluate it. 

## Input

Requires an input TS model including the metal, and a (set of) PDB file(s). The TS ensemble is generated using a starting molfile, and converted to Rosetta parameters using the script `~/rosetta/rosetta_source/src/python/apps/public/molfile_to_params.py`

## Overview of the Procedure

Each Step is included as a directory with a separate README.

1.	Analysis of ZincSite
2. 	CleanPDB after running the zinc analysis. 
3.	Align Transition State model(s) to the Zinc
4.	Minimization in polyAla pocket (optional)
5.	Matching to introduce additional catalytic residues
6.	Design to maximize TS affinity
7.	Revert to Native - Semi-automated refinement
8. 	Dock (validation)

## Details

For most Python scripts, running with a -h option will give a list of available options.

1. Given a PDB file, obtain the co-ordination sphere details of metal (how many protein-metal, and HETATM-metal interactions, in what geometry in each chain).

    `python analyze_zinc_site.py -f 1A4L.pdb`

2. Given an input PDB, "clean" it i.e. keep one chain, change MSE->MET, KCX->LYS etc. Keep HETATMS.

    `python cleanPDBfile.py -f 1A4L.pdb`

3. Align transition state model(s) onto the zinc site according to co-ordination sphere details: 

	1. Given TS model (LG_0001.pdb) and "clean" PDBfile (1A4L_clean.pdb) align them.

	    `python align.py -f 1A4L_clean_A.pdb -l LG_0001.pdb`

	2. Generate Rosetta Inputs given the aligned ligand.

	    `python generate_metal_cstfile.py -f 1A4L_clean_A.pdb -m ZN -a aligned_ligand.pdb`

4. Minimize the constraints energy for the superimposed ligand using Rosetta. Except the metal-chelating protein residues, every other ligand-proximal protein residue is temporarily converted to Ala.

    `~/rosetta/bin/enzyme_design.static.linuxiccrelease @optcst.flags -linmem_ig 10 -in:file::s rosetta_cst.pdb`

5. To introduce additional catalytic residues, run a round of RosettaMatch.

    `~/rosetta/bin/match.linuxiccrelease @general_matching.flags @scaf.flags @subs.flags  -linmem_ig 10 -in:file::s 1A4LA_clean_r.pdb`

6. Design the rest of the pocket for maximizing TS affinity.

    `~/rosetta/bin/enzyme_design.linuxiccrelease @enzdes.flags -correct -linmem_ig 10 -in:file::s <file_from_matching.pdb> > design.log&`

7. For a semi-automated refinement step, generate a so-called resfile, so that the wildtype and any other user-specified residue can be introduced at any given position. Use this resfile with similar commandline as Step 6.

8. To validate the design, perform docking of the ligand in the designed pocket to check if there are any alternative binding modes with similar energies.
