Full instructions for ligand setup, using all web tools. 
This can often be done, and is usally done, using the tools provided by openeye omega, but this tutorial assumes that you have only Rosetta and access to the web.
For instructions for ligand setup using omega (e.g. not with web tools, but with local software on your machine) use the "command" file in this directory.
For more on ligand setup with omega see the manual:
https://www.rosettacommons.org/manuals/archive/rosetta3.4_user_guide/d4/d47/ligand_dock.html 

get 7cpa.pdb from the rcsb. (I like to use fetch command in pymol)
grep out the ligand and the zn:

grep FVF 7cpa.pdb | grep HETATM > 7cpa_fvf.pdb

convert ligand to sdf here:
http://www.webqc.org/molecularformatsconverter.php
(copy/paste the text into the window)
ingput:pdb
output:sdf
no check boxes

free conformer generator from sdf:
http://bioserv.rpbs.univ-paris-diderot.fr/cgi-bin/Frog
input: sdf
output: mol2
Produce: Multi

save the output in 7cpa_fvf_multiout.mol2

python ~/rosetta/rosetta_source/src/python/apps/public/molfile_to_params.py -n CP1 -p CP1 7cpa_fvf_multiout.mol2
mv CP1_* confs
cd confs
cat *.pdb > CP1_lig_confs.pdb
cp CP1_lig_confs.pdb ..
cd ..
emacs CP1.params
put this line:
PDB_ROTAMERS CP1_lig_confs.pdb

Note that here the input conformer is *not* included by default, so even _0001.pdb is not the input.

For more on docking:
Ian Wilson rosettadock paper
docking bench:
1) self-dock, recovery of position and conformer 
2) cross-dock

grep ATOM 7cpa.pdb > 7cpa_prot.pdb
python ~/rosetta/rosetta_source/src/apps/public/relax_w_allatom_cst/clean_pdb_keep_ligand.py 7cpa_prot.pdb -ignorechain
mv 7cpa_prot.pdb_00.pdb 7cpa_prot_fix.pdb

then cat on _0001.pdb to the end to get a useful input structure: 7cpa_input.pdb
cat 7cpa_prot_fix.pdb confs/CP1_0001.pdb > 7cpa_input.pdb 
(To get a native you'd need to use omega with include input)

[path]/rosetta/rosetta_source/bin/ligand_dock.macosgccrelease -database ~/rosetta/rosetta_database/ @flags2 > log.txt &

Note: to get Zn in the docking follow these directions (slightly more advanced);

Add in the Zn line from 7cpa.pdb, but first change HETATM to ATOM and line up the atom and residue names
Now it works if the atom name is left-justified and the residue name is right-justified (check the ZN.params in rosetta_database to confirm this):
ZN__
ATOM   2438 ZN    ZN B 308      44.382  19.543  36.643  1.00 10.29          ZN

Once you have finished this setup, proceed to the ligand_dock tutorial.


