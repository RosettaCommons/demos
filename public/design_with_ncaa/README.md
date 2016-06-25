Design with non-canonical amino acids (NCAA)
============================================

KEYWORDS: DESIGN NONCANONICALS 

In the first part of this tutorial we show how to make a point mutation in a 
PDB file (1l2y) including two non-canonical amino acids (NVL and HLU). 
In the second part of this tutorial we show how to incorporate novel 
non-canonical amino acid side-chains into rosetta for design.

Updated in June 2016 during documentation XRW to enable automated testing (parisah@uw.edu).

Designing with NCAA
===================

Below are two command-lines for making point-mutations including the possibility 
of incorporating non-canonical amino acids. In this procedure the backbone is 
kept fixed and the Trp6 residue is allowed to be altered. *Note: this example
makes use of PyRosetta scripts, and requires that PyRosetta is installed.*

* Trp6 → Canonical amino acids plus NVL and HLU

```bash
$> $ROSETTA3/bin/fixbb.default.linuxgccrelease -s starting_files/1l2y.pdb -use_input_sc -nstruct 1 -ex1 -ex2 -overwrite -minimize_sidechains -resfile rosetta_inputs/ncaa_resfile_pluscanon
```

where `$ROSETTA3`=path-to-Rosetta/main/source

* Trp6 → NVL or HLU

```bash``
$> $ROSETTA3/bin/fixbb.default.linuxgccrelease -s starting_files/1l2y.pdb -use_input_sc -nstruct 1 -ex1 -ex2 -overwrite -minimize_sidechains -resfile rosetta_inputs/ncaa_resfile
```

Adding new sidechains to Rosetta
================================

Generating an initial sidechain model
-------------------------------------
I usually do this in PyMOL but any molecular editor will work. PyMOL has 
rudimentary but useful editing functionality. The goal is to get the 
structure close to the final structure. Don't stress too much because the QM 
minimization will clean up most things. Make sure you double check chirality. 
To make sure the geometry for the atoms near the ends of the molecule is 
correct you will need to put capping groups on the molecule. I would recommend 
an acetyl (ACE) on the nitro-terminus of the betapeptide(BP) and a n-methyl 
(NME) on the carboxy terminus. This structure (ACE-X-NME) is usally called a 
dipeptide despite the fact that it only has one complete residue.

Pro-tips:
- Within PyMOL go to the "Mouse" drop down and select "3 button editing".
- Within PyMOL type "help edit_keys" in the command box will bring up the instructions for the editor. 
- Within PyMOL type "set valence, 1" in the command box to show single, double, triple bonds.
- Be consistent with atom order. It make future steps easier. I like to put all the atoms for the capping groups first, followed by backbone heavy atoms, side chain heavy atoms, and then then hydrogens.

See the example pdb files in the folder `stage_01_initial_structures` and shown 
below:

    <<<<< ornithine.pdb start >>>>>
    ATOM      1  C   ACE     1       0.984   0.045  -0.578  1.00  0.00           C
    ATOM      2  O   ACE     1       1.815  -0.721  -1.083  1.00  0.00           O
    ATOM      3  CH3 ACE     1      -0.445  -0.396  -0.349  1.00  0.00           C
    ATOM      4 1HH3 ACE     1      -0.468  -1.251   0.313  1.00  0.00           H
    ATOM      5 2HH3 ACE     1      -1.013   0.407   0.098  1.00  0.00           H
    ATOM      6 3HH3 ACE     1      -0.904  -0.668  -1.288  1.00  0.00           H
    ATOM     26  N   NME     3       4.937   2.266   0.684  1.00  0.00           N
    ATOM     27  CH3 NME     3       5.341   3.324   1.592  1.00  0.00           C
    ATOM     28  H   NME     3       5.534   1.511   0.316  1.00  0.00           H
    ATOM     29 1HH3 NME     3       4.689   4.177   1.478  1.00  0.00           H
    ATOM     30 2HH3 NME     3       5.286   2.977   2.614  1.00  0.00           H
    ATOM     31 3HH3 NME     3       6.355   3.627   1.379  1.00  0.00           H
    ATOM      7  N   ALA     2       1.392   1.405  -0.195  1.00  0.00           N
    ATOM      8  CA  ALA     2       2.806   1.518  -0.539  1.00  0.00           C
    ATOM      9  C   ALA     2       3.512   2.478   0.390  1.00  0.00           C
    ATOM     10  O   ALA     2       2.929   3.435   0.912  1.00  0.00           O
    ATOM     11  CB  ALA     2       2.895   1.945  -2.014  1.00  0.00           C
    ATOM     12  C01 ALA     2       2.214   0.926  -2.946  1.00  0.00           C
    ATOM     13  N01 ALA     2       1.680   0.434  -5.295  1.00  0.00           N
    ATOM     14  C02 ALA     2       2.334   1.404  -4.405  1.00  0.00           C
    ATOM     15  H   ALA     2       0.770   2.167   0.251  1.00  0.00           H
    ATOM     16  HA  ALA     2       3.283   0.529  -0.413  1.00  0.00           H
    ATOM     17 2HB  ALA     2       2.409   2.923  -2.195  1.00  0.00           H
    ATOM     18 3HB  ALA     2       3.943   2.035  -2.355  1.00  0.00           H
    ATOM     19  H01 ALA     2       2.698  -0.045  -2.840  1.00  0.00           H
    ATOM     20  H02 ALA     2       1.461  -0.405  -4.777  1.00  0.00           H
    ATOM     21  H03 ALA     2       1.161   0.837  -2.679  1.00  0.00           H
    ATOM     22  H04 ALA     2       2.300   0.206  -6.059  1.00  0.00           H
    ATOM     23  H05 ALA     2       3.387   1.492  -4.674  1.00  0.00           H
    ATOM     24  H06 ALA     2       1.850   2.375  -4.509  1.00  0.00           H
    ATOM     25  H07 ALA     2       0.828   0.835  -5.661  1.00  0.00           H
    END
    <<<<< ornithine.pdb end >>>>>

Making GAUSSIAN input
---------------------
We need to minimize the initial structure we made in PyMOL to get a good set of 
ideal bond lengths and angles. To do this we will take the coordinates from the 
pdb file and put them into a gaussian input file. Like the one shown below and 
in stage_02_gaussian_input. A discussion of the complete gaussian input 
structure is beyond the scope of this document. Documentation for gaussian can 
be found at http://www.gaussian.com/.

In short the input is as follows:

- line 1 sets the path for the checkpoint file
- line 2 describes the level of theory, options, and convergence criteria. You may need to change the basis set to something smaller if you are using stuff bellow the 4th line of the periodic table. 
- lines 3-5 are comments
- line 6 is the charge and multiplicity (this is usually 0 1 but ornithine is charged)
- line 7-37 are the elemental type and xyz coordinates of the atoms from the pdb file in the same order as the pdb file
- line 38 is blank
- line 39-40 is the modredundant input that says that we want to keep the torsion formed by atoms 1 13 14 and 15 fixed at 150.00 degrees
- line 41 is blank.

Running Gaussian is simple but the minimizations take a long time (a few hours 
per structure). The command bellow will run gaussian on the input file in the 
`stage_02` folder and put the output in the `stage_03` folder.

    $ g03 stage_02_gaussian_input/ornithine.com stage_03_gaussian_output/ornithine.log

There is a new version of Gaussian called Gaussian09. If you are using 
this version, the commands and options might be slightly different. The input 
below will probably work with Gaussian09 but was written for and tested with 
Gaussian03.

The Gaussian output is very verbose. It's not included in this document but it is 
in the demo folder.

    <<<<< ornithine.com start >>>>>
    %Chk=stage02_ace_nme_res_ordered_pdbs/ornithine.chk
    # HF/6-31G(d) Opt=ModRedundant SCF=Tight Test 

    scan rotamers

    1  1
    C        0.984   0.045  -0.578
    O        1.815  -0.721  -1.083
    C       -0.445  -0.396  -0.349
    H       -0.468  -1.251   0.313
    H       -1.013   0.407   0.098
    H       -0.904  -0.668  -1.288
    N        4.937   2.266   0.684
    C        5.341   3.324   1.592
    H        5.534   1.511   0.316
    H        4.689   4.177   1.478
    H        5.286   2.977   2.614
    H        6.355   3.627   1.379
    N        1.392   1.405  -0.195
    C        2.806   1.518  -0.539
    C        3.512   2.478   0.390
    O        2.929   3.435   0.912
    C        2.895   1.945  -2.014
    C        2.214   0.926  -2.946
    N        1.680   0.434  -5.295
    C        2.334   1.404  -4.405
    H        0.770   2.167   0.251
    H        3.283   0.529  -0.413
    H        2.409   2.923  -2.195
    H        3.943   2.035  -2.355
    H        2.698  -0.045  -2.840
    H        1.461  -0.405  -4.777
    H        1.161   0.837  -2.679
    H        2.300   0.206  -6.059
    H        3.387   1.492  -4.674
    H        1.850   2.375  -4.509
    H        0.828   0.835  -5.661

    1 13 14 15  -150.00 F
    13 14 15 7  150.00 F
    <<<<< ornithine.com end >>>>>

Converting GAUSSIAN output to molfile
-------------------------------------
The program bable we built in the first step can convert the gaussian output to 
a molfile.

Bable says it can handle output form Gaussian09. I am not sure when it was 
added, the current version is 2.3.0 and the modified version I have included is 
based on version 2.2.0. The changes to 2.2.0 were not extensive and could 
probably be ported to 2.3.0.

The command bellow will convert the gaussian output to molfile format for the 
peptide and peptoid examples:

    ```
    $ openbabel/install/bin/babel -i g03 stage_03_gaussian_output/ornithine.log -o mol stage_04_molfile/ornithine.mol
    $ openbabel/install/bin/babel -i g03 stage_03_gaussian_output/amino2.log -o mol stage_04_molfile/amino2.mol
    ```

    <<<<< ornithine.mol start >>>>>
     OpenBabel01101117083D

     31 30  0  0  0  0  0  0  0  0999 V2000
        0.4866    2.2296   -0.2354 C   0  0  0  0  0
        1.2081    1.9315   -1.1549 O   0  0  0  0  0
        0.4369    3.6268    0.3352 C   0  0  0  0  0
       -0.3187    4.1933   -0.1996 H   0  0  0  0  0
        0.1852    3.6362    1.3888 H   0  0  0  0  0
        1.3920    4.1075    0.1799 H   0  0  0  0  0
       -2.6920   -1.1572   -0.6650 N   0  0  0  0  0
       -4.0550   -1.5847   -0.3922 C   0  0  0  0  0
       -2.3500   -1.2554   -1.5934 H   0  0  0  0  0
       -4.1238   -1.9440    0.6231 H   0  0  0  0  0
       -4.7601   -0.7727   -0.5220 H   0  0  0  0  0
       -4.3066   -2.3874   -1.0709 H   0  0  0  0  0
       -0.3415    1.3293    0.3431 N   0  0  0  0  0
       -0.6023    0.0286   -0.2270 C   0  0  0  0  0
       -2.0235   -0.3660    0.1877 C   0  0  0  0  0
       -2.4547   -0.0133    1.2526 O   0  0  0  0  0
        0.3558   -1.0672    0.2829 C   0  0  0  0  0
        1.8184   -0.8036   -0.0848 C   0  0  0  0  0
        4.1559   -1.6146    0.0004 N   0  0  0  0  0
        2.7190   -1.9434    0.3678 C   0  0  0  0  0
       -0.9580    1.6065    1.0765 H   0  0  0  0  0
       -0.5197    0.1005   -1.3044 H   0  0  0  0  0
        0.2486   -1.1369    1.3608 H   0  0  0  0  0
        0.0400   -2.0197   -0.1352 H   0  0  0  0  0
        1.9103   -0.6656   -1.1566 H   0  0  0  0  0
        4.2542   -1.4797   -0.9970 H   0  0  0  0  0
        2.1487    0.1188    0.3755 H   0  0  0  0  0
        4.7999   -2.3436    0.2758 H   0  0  0  0  0
        2.4931   -2.8811   -0.1191 H   0  0  0  0  0
        2.7087   -2.0883    1.4385 H   0  0  0  0  0
        4.4568   -0.7571    0.4437 H   0  0  0  0  0
      2  1  2  0  0  0
      1  3  1  0  0  0
      1 13  1  0  0  0
      4  3  1  0  0  0
      6  3  1  0  0  0
      3  5  1  0  0  0
      9  7  1  0  0  0
      7  8  1  0  0  0
      7 15  1  0  0  0
     12  8  1  0  0  0
     11  8  1  0  0  0
      8 10  1  0  0  0
     14 13  1  0  0  0
     13 21  1  0  0  0
     22 14  1  0  0  0
     14 15  1  0  0  0
     14 17  1  0  0  0
     15 16  2  0  0  0
     24 17  1  0  0  0
     18 17  1  0  0  0
     17 23  1  0  0  0
     25 18  1  0  0  0
     18 20  1  0  0  0
     18 27  1  0  0  0
     26 19  1  0  0  0
     19 28  1  0  0  0
     19 20  1  0  0  0
     19 31  1  0  0  0
     29 20  1  0  0  0
     20 30  1  0  0  0
    M  END
    <<<<< ornithine.mol end >>>>>

Modifying the molfiles
----------------------
The `molfile2params_polymer.py` script requires some additional data to be added 
to the end of the molfile. This data is specified at the end of the file after 
the bond information. It is a list of variable names and then a list of values. 
The variable are described bellow. 

* ROOT: Single numerical value. Atom number (acording to the order above). 
  Where the atom tree is rooted for this residue type. Should probably be the 
  nitrogen of the central residue.

* POLY_N_BB, POLY_CA_BB, POLY_C_BB, POLY_CO_BB: Single numerical value. The 
  backbone nitrogen, alpha-carbon, carbonyl-carbon and carbonyl-oxygen. These 
  get special rosetta atom types and so are listed here special. Note: You will 
  need to add some additional code to work with your additional backbone atom.

* POLY_IGNORE: List of numerical values.  Atom number (acording to the order in 
  the molfile). These are the atoms for the capping groups with the exception of 
  the upper and lower connect atoms. They will not be listed in the atoms and 
  bonds in the params file but are used in determining the atom types.

* POLY_UPPER, POLY_LOWER: Single numerical value. Atom number (according to the 
  order in the molfile). These are the atoms in the capping groups that connect 
  to the residue. They will not be listed in the atoms and bonds in the params 
  file but are used in determining the atom types and they are listed in the 
  internal coordinate section.

* POLY_CHG: Single numerical value. Overall charge on the residue. 

* POLY_PROPERTIES: List of alpha-numerical values. These get used by rosetta at 
  various places in the program. You can say something like:
  ```
  if (pose.residue(10).type().is_protein() ) { /* do someting */ }
  ```
  You will want to create a new property called "BETAPEPTIDE" and corresponding 
  functions to the residue type class. 

* END: The end of the file.

Note that there are 2 spaces between the "M" and the variable name. If you have to 
make multiple residue types (which it sounds like you will), keeping the atoms 
in the same order makes assigning all these numbers easier. 

Additional info for peptide ornithine:

    M  ROOT 13
    M  POLY_N_BB 13
    M  POLY_CA_BB 14
    M  POLY_C_BB 15
    M  POLY_O_BB 16
    M  POLY_IGNORE 2 3 4 5 6 8 9 10 11 12
    M  POLY_UPPER 7
    M  POLY_LOWER 1
    M  POLY_CHG 1
    M  POLY_PROPERTIES PROTEIN POLAR CHARGED
    M  END

Additional info for amino2:

    M  ROOT 1
    M  POLY_N_BB 1
    M  POLY_CA_BB 2
    M  POLY_C_BB 3
    M  POLY_O_BB 4
    M  POLY_IGNORE 9 10 22 23 24 27 28 29 30 31 32 33 34
    M  POLY_UPPER 26
    M  POLY_LOWER 7
    M  POLY_CHG 1
    M  POLY_PROPERTIES PEPTOID POLAR CHARGED
    M  END

Running molfile2params.py
-------------------------
Finally, the `molfile2params_polymer.py` script will convert the modified molfile 
to a params file. The commands bellow work for the example files:

```bash
> ./molfile_to_params_polymer.py --clobber --polymer --no-pdb --name C40 -k ornithine.kin stage_05_modified_molfile/ornithine.mol
> ./molfile_to_params_polymer.py --clobber --polymer --peptoid --no-pdb --name P01 -k amino2.kin stage_05_modified_molfile/amino2.mol
```

There is additional tweaking that needs to happen to the params files to make 
them work. Compare these to the ones in the database for reference. 

The `molfile2params.py` script can produce a kinemage file. This is super handy 
as it lets you check how Rosetta will build the atom tree, and the atom type 
assignments, and other stuff. You can open kinemage files using the KiNG 
program from the Richardson lab at Duke 
(http://kinemage.biochem.duke.edu/software/king.php).

Using your shiny new params file
================================
Rosetta params files live in the database at database/chemical/residue_type_sets/(fa_standard)/..., but you can just pass your new parameters in on command line via `-extra_res_fa` or `-extra_res` (the latter for centroid).

Rotamer library design
======================

For questions on backbone-dependent rotamer library design please contact 
either Doug Renfrew (dougrenfrew at gmail dot com) or Tim Craven (twc254 at nyu dot edu).

