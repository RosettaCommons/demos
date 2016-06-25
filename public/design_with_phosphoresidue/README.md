Design with Phosphoresidues
===========================

KEYWORDS: LIGANDS NONCANONICALS STRUCTURE_PREDICTION

Phosphoresidues are very important in biology.  Many signaling pathways depend 
on the fact that phosphorylation of residues like serine, threonine, and 
tyrosine can introduce dramatic allosteric changes.  This demo shows how to 
use phosphoresidues in Rosetta.

Make sure the phosphoresidue is read properly
---------------------------------------------

We create the parameter file for the phosphoresidue as follows:

1.  Extract a phosphotyrosine (PTR) residue from a PDB file

    We use one model of the NMR ensembles of the 2lct.pdb by deleting other 
    models and keep the first one (3lct.pdb refers to first model). From the 
    "3lct.pdb" file, extract the part containing the atoms of the 
    phosphoresidues 342 and 346 to generate "PT1.pdb" and "PT2.pdb". You can 
    use the sample command line

        grep PTR 3lct.pdb | grep 342 | grep HETATM > PT1.pdb
        grep PTR 3lct.pdb | grep 346 | grep HETATM > PT2.pdb

2.  Create a MOLFILE for PTR

    Either use avogadro(free chemical software) or go to the adress: 
    www.molecular-networks.com/online_demos/convert_demo to generate an 
    "PT1.mdl" file from the "PT1.pdb" with the original geometry of the 
    phosphoresidue.

3.  Create a Rosetta PARAMS file PTR

    Copy the "PTR.mdl" file (molfile format) into your working directory.  Then 
    run: (where `$ROSETTA3`=path-to-Rosetta/main/source)
```
        $ROSETTA3/scripts/python/public/molfile_to_params.py PT1.mdl -n PT1
        $ROSETTA3/scripts/python/public/molfile_to_params.py PT2.mdl -n PT2
```
This script creates PT1.params and PT1_0001.pdb. You'll need to delete the older PTR residues in 3lct.pdb and replace them with the co-ordinates for generated new residues PT1_0001.pdb and PT2_0001.pdb (4lct.pdb refers to 3lct.pdb with new phoshoresidues). Also, change the "HETATM" tags for phoshoresidues in the PDB file to "ATOM". We have provided example files for you in the rosetta_inputs directory.

Design around the phosphoresidue
--------------------------------

1.  Using pyMol, identify the neighbouring residues, within the desired radius 
    (5 A in our case, PRT.resfile)

2.  Write the resfile. Use the option ALLAA next to the sequence position to 
    design the selected residue to all of posibble amino acids. We prepared a sample resfile for you that includes all residues.

3.  Run the following command-line:
```
        $> $ROSETTA3/bin/fixbb.default.macosgccrelease -s rosetta_inputs/4lct.pdb  -extra_res_fa rosetta_inputs/PT1.params rosetta_inputs/PT2.params -resfile rosetta_inputs/PTR.resfile
```
This will generate a single pdb file with designed residues around the phosphoresidues. The sample outputs have been copied to /outputs. If you are using this protocol for a real design application, use the flag `-nstruct 1000` or `-nstruct 10,000` to generate enough designs to give meaningful results.



