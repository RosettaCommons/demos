#!/usr/bin/python
## make mammoth structure alignments


import string
from glob import glob
from sys import argv,stderr
from os import popen,system
from os.path import exists
from operator import add

#############################

if len(argv) <=2:
    print '\n'
    print '-'*75
    print 'USAGE: %s <pdb1> <pdb2> {... <pdbN>} > <superposition-pdb>'%argv[0]
    print '\n will superimpose pdbs 2-N onto pdb1 using maxsub, so seqs should agree'
    print '-'*75
    print '\n\n'
    assert 0==1


RENUMBER_ATOMS = 0
SHOW_MODEL_0 = 1

model_count = 0
atom_count = 0

args = argv[1:]
if args.count('-R'):
    pos = args.index('-R')
    rmsd_threshold = float(args[pos+1])
    del args[pos]
    del args[pos]
else:
    rmsd_threshold = 0.0

if args.count('-1'):
    del args[args.index('-1')]
    SHOW_MODEL_0 = 0


pdb_list = args
pdb1 = pdb_list[0]

for pdb in pdb_list[1:]:
    if rmsd_threshold:
        command = '/work/pbradley/maxsub/maxsub -R %f -p %s -e %s 2> /dev/null | grep PSI.end'\
                  %(rmsd_threshold,pdb1,pdb)
    else:
        #command = '/work/pbradley/maxsub/maxsub -p %s -e %s  | grep PSI.end'\
        #          %(pdb1,pdb)
        command = '/work/pbradley/maxsub/maxsub -p %s -e %s 2> /dev/null | grep PSI.end'\
                  %(pdb1,pdb)


    lines = popen(command).readlines()

    if not lines:
        stderr.write('empty file? %s\n'%pdb)
        continue

    l = string.split(lines[0])
    stderr.write('%s -vs- %s: %s over %s residues\n'%(pdb_list[0],pdb,l[7],l[3]))



    file = 'maxsub_sup.pdb'
    matrix = map(lambda x:map(float,string.split(x)[1:]),
                 popen('grep -A3 "Transformation Matrix" %s'%file).readlines()[1:])


    P_translation = map(float,
                        string.split(popen('grep -A1 "Translation vector (Pred" %s'\
                                           %file).readlines()[1])[1:])

    E_translation = map(float,
             string.split(popen('grep -A1 "Translation vector (Exp" %s'\
                                %file).readlines()[1])[1:])


    def E_transform(v,matrix,tP,tE):
        ans = [0.0]*3
        for i in range(3):
            for j in range(3):
                ans[i] = ans[i] + matrix[i][j]*(v[j]-tE[j])
            ans[i] = ans[i] + tP[i]

        return ans


    if model_count == 0:
        if SHOW_MODEL_0:
            print 'MODEL     %4d'%model_count
            model_count = model_count+1

            data = open(pdb1,'r')
            line = data.readline()
            while line:
                if line[:4] in ['ATOM','HETA']:
                    atom_count = atom_count + 1
                    print '%s  %5d%s'%(line[:4],atom_count,line[11:-1])
                elif line[:6] == 'ENDMDL':break
                line = data.readline()
            data.close()
            print 'ENDMDL'
        else:
            model_count = model_count + 1

    print 'MODEL     %4d'%model_count
    model_count = model_count+1
    data = open(pdb,'r')
    line = data.readline()
    while line:
        if line[:4] in ['ATOM','HETA']:
                atom_count = atom_count + 1
                pos = E_transform(map(float,[line[30:38],line[38:46],line[46:54]]),
                                  matrix,
                                  P_translation,
                                  E_translation)

                if RENUMBER_ATOMS:
                    print '%s  %5d%s%8.3f%8.3f%8.3f%s'\
                          %(line[:4],atom_count,line[11:30],pos[0],pos[1],pos[2],line[54:-1])
                else:
                    print '%s  %s%8.3f%8.3f%8.3f%s'\
                          %(line[:4],line[6:30],pos[0],pos[1],pos[2],line[54:-1])


        elif line[:6] == 'ENDMDL':break
        line = data.readline()
    data.close()

    print 'ENDMDL'


