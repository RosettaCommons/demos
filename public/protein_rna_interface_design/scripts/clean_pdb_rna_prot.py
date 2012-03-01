#!/usr/bin/python

import sys, os, string

pdb = sys.argv[1]

# Puts the correct number of spaces before ATOMNO or RESNO fields
def spaces(str, sp):
    numspace = sp - len(repr(str))
    return ' ' * numspace + repr(str)

# Lists containing lines from original PDB
a = [] # Protein
b = [] # Nucleic acid

f = open(pdb, 'rt')
for l in f:
    l = l.rstrip('\n')
    if len(l) < 26:
	    continue
    elif l[0:4] != 'ATOM':
	    continue
	
    resn = l[17:20].strip()
    if len(resn) == 3: # If RES name has three characters, it's a protein
        a.append(l)
    else:              # Otherwise, it's assumed to be a nucleic acid.
        b.append(l)
f.close()

atomno = 1             # Running count of atom number
a_res = 0              # Running count of residue number in new chain A
cur_chres = ''         # Stores last chain ID/residue index combination in old chain

for l in a:
    if cur_chres != l[17:22]:
	    cur_chres = l[17:22]
	    a_res += 1
    print l[0:4] + spaces(atomno, 7) + l[11:21] + 'A' + spaces(a_res, 4) + l[26:]
    atomno += 1

print 'TER'

b_res = 0              # Running count of residue number in new chain B
cur_chres = ''         # Stores last chain ID/residue index combination in old chain

for l in b:
    if cur_chres != l[17:22]:
	    cur_chres = l[17:22]
	    b_res += 1
    print l[0:4] + spaces(atomno, 7) + l[11:21] + 'B' + spaces(b_res, 4) + l[26:]
    atomno += 1

print 'END'