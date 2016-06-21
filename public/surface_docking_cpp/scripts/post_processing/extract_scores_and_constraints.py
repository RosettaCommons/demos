#!/usr/bin/python2.6
# Script to extract scores from score table and constraint scores appended to end of PDB
# Emily Koo

import os
from sys import argv

def parse(pdbs):
    savefile = pdbs.name + "_constraints"
    save = open(savefile, 'a')
    for pdb in pdbs:
        decoy = pdb.strip()
        file = open(decoy, 'r')
        for line in file:
            if line.startswith('Constraint Energies'):
                newline = file.next()
                # Atom pair 4
                energy = float(newline.split()[3])
                # Dihedral 3
                newline = file.next()
                energy += float(newline.split()[2])
                save.write(decoy + "\t" + str(energy) + "\n")
                break
            
    save.close()

state = str(argv[1]).lower()

sol = 'ls SolState*.pdb > sol_pdb'
ads = 'ls AdsState*.pdb > ads_pdb'
get = "grep 'Total weighted' SolState*.pdb | awk '{print $1, $5}' | sed s/://g | sort > SolState.sorted"
get2 = "grep 'Total weighted' AdsState*.pdb | awk '{print $1, $5}' | sed s/://g | sort > AdsState.sorted"
join = 'join SolState.sorted sol_pdb_constraints > table_sol_scores'
join2 = 'join AdsState.sorted ads_pdb_constraints> table_ads_scores'
rm = 'rm sol_pdb_constraints SolState.sorted sol_pdb'
rm2 = 'rm ads_pdb_constraints AdsState.sorted ads_pdb'
rename = 'cp SolState.sorted table_sol_scores'
rename2 = 'cp AdsState.sorted table_ads_scores'

if state == "sol":
    os.system(sol)
    sol_pdb = open('sol_pdb', 'r')
    parse(sol_pdb)
    sol_pdb.close()
    os.system(get)
    if os.stat(sol_pdb.name + "_constraints").st_size != 0:
        os.system(join)
    else:
        os.system(rename)
    os.system(rm)
elif state == "ads":
    os.system(ads)
    ads_pdb = open('ads_pdb', 'r')
    parse(ads_pdb)
    ads_pdb.close()
    os.system(get2)
    if os.stat(ads_pdb.name + "_constraints").st_size != 0:
        os.system(join2)
    else:
        os.system(rename2)
    os.system(rm2)
else: #both
    os.system(sol)
    os.system(ads)
    sol_pdb = open('sol_pdb', 'r')
    ads_pdb = open('ads_pdb', 'r')
    parse(sol_pdb)
    parse(ads_pdb)
    sol_pdb.close()
    ads_pdb.close()

    os.system(get)
    os.system(get2)
    if os.stat(sol_pdb.name + "_constraints").st_size != 0:
        os.system(join)
    else:
        os.system(rename)
    if os.stat(ads_pdb.name + "_constraints").st_size != 0:
        os.system(join2)
    else:
        os.system(rename2)
    os.system(rm)
    os.system(rm2)
    
print "Done"