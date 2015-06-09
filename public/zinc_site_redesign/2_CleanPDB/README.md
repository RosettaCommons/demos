Here other chains are removed  and residues like MSE are changed to
MET, remove structural metal sites, and change KCX to LYS. The default
chain that is kept is chain A. The script is executed by 

python cleanPDBfile.py -f 1A4L.pdb 

To keep chain B for example:

python cleanPDBfile.py -f 1A4L.pdb -c B

Other options: To remove a specific metal ion from PDB

-m: name of metal to remove from pdb
-n: residue number of metal to remove from pdb

python cleanPDBfile.py -f PDBFILE -m METAL -n NR

