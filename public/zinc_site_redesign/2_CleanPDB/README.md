Here other chains are removed  and residues like MSE are changed to
MET, remove structural metal sites, and change KCX to LYS. The default
chain that is kept is chain A. The script is executed by 
KEYWORDS: METALS UTILITIES
```bash
python $ROSETTA_TOOLS/zinc_site_redesign/cleanPDBfile.py -f 1A4L.pdb 
```
The expected output from the above is shown in expected_output_1A4L_clean_A.pdb.

To keep chain B for example:

```bash
python $ROSETTA_TOOLS/zinc_site_redesign/cleanPDBfile.py -f 1A4L.pdb -c B
```

Other options: To remove a specific metal ion from PDB

-m: name of metal to remove from pdb
-n: residue number of metal to remove from pdb

```bash
python $ROSETTA_TOOLS/zinc_site_redesign/cleanPDBfile.py -f PDBFILE -m METAL -n NR
```

