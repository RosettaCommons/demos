# Preparing Ligands

KEYWORDS: CORE_CONCEPTS LIGANDS

Tutorial by Parisa Hosseinzadeh (parisah@uw.edu), prepared in June 2016.

[[_TOC_]]

Rosetta is residue-based. That is, Rosetta doesn't treat atoms individually, but rather as a complete chemical entity. For this to work, bond lengths, bond angles, and other properties of the residue must be known to Rosetta. This is true for both the standard amino acid residues as well as other small molecules and residues it models. Rosetta encodes the properties of a residue in **params files**. The Rosetta database contains params files for the common amino acids, as well as for a selection of other small molecules. In order to add a new ligand - one Rosetta does not recognize - you will need to generate a new params file.

#### What is a Params File

A **params** file is where the pre-computed information about the geometry and chemical features of residues and ligands are stored. In order for Rosetta to understand how to treat the residue during any run, there should be a params file associated with that residue. You can find examples of params file for many types of molecules such as 20 canonical amino acids, some non-canonical amino acids, water, metal ions, nucleic acids, etc. in Rosetta/main/database/chemical/residue_type_sets/fa_standard/residue_types. More detailed information about what each line means can be found [in the Rosetta documentation](https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/Residue-Params-file). In this tutorial, we go over how to prepare a ligand a generate a params file for it.

#### Downloading and Preparing the Ligand

Let's start from obtaining the ligand. The ligand we are working with today is going to be "serotonin", the molecule of happiness! There are many ways to download the PDB structure of a ligand. The easiest way is through the PDB site itself. Go to http://www.rcsb.org/pdb/home/home.do and choose search -> Advanced Search. Start advanced search. In choose a query type, scroll down to **chemical components** and choose **chemical name**. Type Serotonin and submit query. You can see that some options will show up. Click on "__SRO__" with formula C10 H12 N2 O. You can see that "SRO is found in 9 entries". Click on _9 entries_. It will take you to a new page, showing all the PDB entries with serotonin in them. 

Click on __Download Files__ tab on the top. Scroll down to __Download Ligands__. In *Download Options* choose : "Single SDF file". In *Instance Type* choose: "All instances from each PDB entry". Leave the rest of the filters at their default values. Make sure you only chose SRO in the ligand box. Then **Launch the Download Application**. A file named "Ligands_noHydrogens_noMissing_27_Instances.sdf" or something similar should be downloaded. Please refer to [More Advanced Preparation](#More-Advanced-Preparation) for information on how to generate conformers from a single PDB.

Navigate to the prepare_ligand directory running this command:

```
> cd $ROSETTA_DEMOS/tutorials/prepare_ligand/inputs
```

We have already provided you with the same file. You can open them and compare them, but it is renamed to serotonin.sdf. You can rename the file running this command:

```
> mv <file_name>.sdf serotonin.sdf
```

Now, at this step, we need to add hydrogens to this molecule. There are multiple ways to achieve this goal. If installed on your computer, you can use [babel](http://openbabel.org/wiki/Main_Page) to add hydrogens with the GUI interface or using command line:

```
> babel -h serotonin.sdf serotonin_withH.sdf
```

You can also use [Avogadro](http://avogadro.cc/wiki/Main_Page) and use it to add hydrogens and minimize the PDB.

The above mentioned software is free, but other software such as [mercury](https://www.ccdc.cam.ac.uk/solutions/csd-system/components/mercury/) or [omega](http://www.eyesopen.com/omega) can also be used, if you have them available.

We have provided a serotonin_withH.sdf file for you in the inputs directory if you cannot get it. 

#### Generating the Params File

Now we are at the stage of generating the params file. Rosetta provides a python script that will take an sdf file and generate a params file from it.

Use `cd ../` to get out of the inputs directory and run this command:

```
$> $ROSETTA_MAIN/source/scripts/python/public/molfile_to_params.py -n SRO -p SRO --conformers-in-one-file inputs/serotonin_withH.sdf
```

*-n* is a option that says what is the three letter name you want to give to the ligand and *-p* defines the name given to the params file. After you are done with this command you should have a file named SRO.params like the one provided in the outputs directory. Also, 27 PDB files (or any number you saw during download) are generated based on each conformer found from PDB.

>**Note1.** To learn better what the molfile_to_params.py works and what are the other options, you can run:

```
> $ROSETTA_MAIN/source/scripts/python/public/molfile_to_params.py --help
```

>**Note2**. If you have multiple files that you want to generate the params file for, you can run a batch molfile_to_params command as follows:

```
> $ROSETTA_MAIN/source/scripts/python/public/batch_molfile_to_params.py --script_path=<path-to-molfile_to_params.py> list_of_molfiles.txt
```

Now let's take a look at the generated params file:

```
NAME SRO
IO_STRING SRO Z
TYPE LIGAND
AA UNK
```

These lines are basic information about the ligand. We are telling Rosetta that this residue is a `LIGAND`, its name is `SRO` and it is shown by character `Z` in one letter code. This residue was `UNK`own to Rosetta prior to this.

```
ATOM  C6  aroC  X   -0.10
ATOM  C5  aroC  X   -0.10
ATOM  N1  Ntrp  X   -0.59
ATOM  C4  aroC  X   -0.10
ATOM  C3  aroC  X   -0.10
ATOM  C2  aroC  X   -0.10
ATOM  C1  aroC  X   -0.10
ATOM  O1  OH    X   -0.64
ATOM  H1  Hpol  X   0.45
ATOM  C8  aroC  X   -0.10
ATOM  C7  aroC  X   -0.10
[...]
```

These lines show the atoms in the residue. The first column are the atom naming that are going to be used in the PDB. The second column shows the Rosetta atom types for each atom. For example, `Hpol` means a polar hydrogen that can be involved in hydrogen bonding interactions whereas `Haro` is a hydrogen that is attached to an aromatic ring. The last column shows partial charges. Go through these columns carefully and check whether the definitions are correct. If you have calculated partial charges from some molecular dynamics sources, change these with those. Make sure atom types are correctly defined.

Next lines start with `BOND_TYPE` and define what atoms are connected to each other and whether this is a single bond or a double bond. `CHI` column defines a rotatable chi angle. A chi angle number is specified, followed by the names of the four atoms defining the angle. `PROTON_CHI 1` shows how much sampling the proton in that chi is going to go through.

Next two lines are `NBR_ATOM` and `NBR_RADIUS`

```
NBR_ATOM  C6
NBR_RADIUS 5.878178
```

`NBR_ATOM` is "neighbor atom" in Rosetta. For a ligand, it is the closest atom to geometric center of mass by default. `NBR_RADIUS` is a measure of the size of the ligand. It's an estimate of the furthest possible distance from the neighbor atom to any other heavy atom in the residue.

The final lines are `ICOOR_INTERNAL`, and define the internal coordinates of the atoms (see [here](https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/Residue-Params-file) for more information.

At the end of your file, add this line:

```
PDB_ROTAMERS SRO_conformers.pdb
```

This line tells Rosetta where the conformers are stored.

#### Using Params File during the Run

Now that you generated the params file, you should inform Rosetta where it is. This can be done by adding this option to the your command line:

```
-extra_res_fa <path-to-params-file>
```

#### More Advanced Preparation

Often times, you need to generate conformers for a ligand (i.e. different conformations a ligand can take that are energetically favorable). Rosetta cannot do this, but there are different software packages that can perform this function. After the conformers are generated, the rest of the process is the same and you can continue from [Generating the Params File](#Generating-the-Params-File).

-   [babel](http://open-babel.readthedocs.io/en/latest/3DStructureGen/multipleconformers.html)
-   [Avogadro](http://manual.avogadro.cc/content/7-optimizing-geometry/2-conformers.html)
-   [CSD Mercury](https://www.ccdc.cam.ac.uk/solutions/csd-system/components/mercury/): CSD_Discovery -> Conformer Generation
-   [Omega OpenEye](https://docs.eyesopen.com/omega/usage.html)

#### Other Residue Types
Rosetta can also be used with D-amino acids, non-canonical amino acids, non-canonical backbones, RNA and metal ions. [Here](https://www.rosettacommons.org/docs/latest/rosetta_basics/non_protein_residues/non-protein-residues) is documentation describing how to prepare them.
