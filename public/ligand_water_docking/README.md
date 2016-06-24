Ligand-centric Water Docking
============================

KEYWORDS: LIGANDS DOCKING

Author: Gordon Lemmon  
Citation:
* Gordon Lemmon, Jens Meiler (2012). Toward ligand docking including explicit 
  interface water molecules. PLoS ONE (submitted).

---

Small molecule docking predicts the interaction of a small molecule ligand with 
a protein at atomic-detail accuracy including position and conformation the 
ligand but also conformational changes of the protein upon ligand binding. 
While successful in the majority of cases, leading docking algorithms including 
RosettaLigand fail in some cases to predict the correct protein/ligand complex 
structure. In this study we show that simultaneous docking of explicit 
interface water molecules greatly improves Rosetta’s ability to distinguish 
correct from incorrect ligand poses. This result holds true for both 
protein-centric water docking wherein waters are located relative to the 
protein binding site and ligand-centric water docking wherein waters move with 
the ligand during docking. Protein-centric docking is used to model 99 HIV-1 
protease/protease inhibitor structures. We find protease inhibitor placement 
improving at a ratio of 9:1 when one critical interface water molecule is 
included in the docking simulation. Ligand-centric docking is applied to 341 
structures from the CSAR benchmark of diverse protein/ligand complexes. Across 
this diverse dataset we see up to 56% recovery of failed docking studies, when 
waters are included in the docking simulation.

Purpose and algorithm
---------------------

The simultaneous docking of ligands and water molecules is now possible within Rosetta
This work is being submitted to the PLoS ONE Special collection from RosettaCon 2012

Waters can be docked in two ways. The first way is refered to in the paper as protein-centric 
water docking. Protein centric waters sample the protein binding pocket independent of the ligand.
Ligand-centric waters are positioned initially around the ligand. As the ligand translates and
rotates about the protein binding site, the waters move with it, maintaining their positions relative
to the protein. After low-resolution docking of the ligand. The waters undergo their own independent
but smaller movements.

Tools and Input Files
---------------------

#### Scripts:

* Use `find_waters_pymol.py` to identify waters within the protein/ligand 
  interface. The script takes 3 arguments, each specifying a PDB file 
  (protein.pdb, ligand.pdb, water.pdb). This script requires that you have 
  installed pymol as a python library

* The script `ligand_properties_from_bcl` requires that you download and 
  install BCL, available from the Meiler Lab website, www.meilerlab.org

#### Required flags:

    -in:file:s <pdb file> # a starting structure upon which docking will be performed. Should contain a protein, a ligand, and one or more waters
    -in:file:extra_res_fa # the .params file for your ligand. This is created by providing a .mol file to the script: rosetta_source/src/python/apps/public/molfile_to_params.py
    -treat_residues_with_this_chain_as_separate_chemical_entities <1 letter chain from PDB> # Useful for giving waters the same chain and having Rosetta treat them separately.

#### Optional flags:

    -ex1, -ex1aro, and -ex2 # expand the rotamer sets that are sampled during packing.
    -in:file:native # allows calculation of comparison metrics between Rosetta models and the correct pose if this is known

#### Example Rosetta Command Line (Use Rosetta3.5 or revision 48472):

The three XML files, standard.xml, protein_centric.xml and ligand_centric.xml, demonstrate how to dock waters using
RosettaLigand. Simply run the following command from the directory where this readme is found:

    <path/to/rosetta_source>/bin/rosetta_scripts.linuxgccrelease @rosetta_inputs/flags.txt -parser:protocol <one_of_the_xmls>

For example, to run all three xmls without overwriting the output from the previous xml, you can run
    $> <path/to/rosetta_source>/bin/rosetta_scripts.linuxgccrelease @rosetta_inputs/flags.txt -parser:protocol ./rosetta_inputs/standard.xml -out:suffix _standard
    $> <path/to/rosetta_source>/bin/rosetta_scripts.linuxgccrelease @rosetta_inputs/flags.txt -parser:protocol ./rosetta_inputs/protein_centric.xml -out:suffix _protein_centric
    $> <path/to/rosetta_source>/bin/rosetta_scripts.linuxgccrelease @rosetta_inputs/flags.txt -parser:protocol ./rosetta_inputs/ligand_centric.xml -out:suffix _ligand_centric
    

Expected Outputs
----------------

This demo only produces one structure. Add "-nstruct 1000" to produce the amount of sampling used in this paper.
In our paper we sort the top 1000 by total score, then the top 100 by interface score. The top model by interface
score is used as our most likely prediction of ligand pose

To sort models by total score simply:

    grep -H  total_score *.pdb | sort -nk 2 | head -n 100 | cut -d ':' -f 1 > top100.txt

To sort by interface score:

    grep -H interface_delta_X *.pdb `cat top100.txt` | sort -nk 2 | head -n 1
