# run from pymol: pymol -qcr <mk.py> -- <pdb-file> <resi> [<resn> ...]
# generates mutatuted pdbs of <pdb-file> at residue number <resi> for given targets
# <resi> can start with 1-char abbrev of native AA for checking
# defaults to all but native residue if no target given
import sys,pymol
from pymol import cmd,stored

aal = [ "GLY","ALA","CYS","VAL","MET","LEU","ILE","PHE","TYR","TRP","PRO","HIS","SER","THR","ASN","GLN","ASP","GLU","ARG","LYS"]
aas = [ "G"  ,"A"  ,"C"  ,"V"  ,"M"  ,"L"  ,"I"  ,"F"  ,"Y"  ,"W"  ,"P"  ,"H"  ,"S"  ,"T"  ,"N"  ,"Q"  ,"D"  ,"E"  ,"R"  ,"K"  ]

fn = sys.argv[1]    # filename
knp = sys.argv[2]   # resi, perhaps preceeded by short resn
if len(sys.argv) > 3:
    n = len(sys.argv) - 3;
    ars = sys.argv[3:3+n]
    arl = []
    for i in range(n):
        arl.append(aal[aas.index(ars[i])])
else:
    ars = aas
    arl = aal

# open file
bn = fn.replace(".pdb", "")
cmd.load(fn)
# deal with possible starting short resn
if knp[0].isalpha():
    natl = aal[aas.index(knp[0])]
    kn = knp[1:]
else:
    natl = ""
    kn = knp
# find native residue
stored.nat = []
cmd.iterate('resi %s & n. ca' % kn, 'stored.nat.append((model,chain,resn))')
if len(stored.nat) == 0:
    print "no native residue %s found, exiting" % kn
    exit(1)
elif len(stored.nat) > 1:
    for i in range(1, len(stored.nat)):
        if stored.nat[0][2] != stored.nat[i][2]:
            print "residue mismatch: %s/%s/%s != %s/%s/%s" % (stored.nat[0], stored.nat[i])
            exit(1)
elif natl != "" and stored.nat[0][2] != natl:
    print "native residue mismatch at resi %s: expected %s, found %s" % (kn, natl, stored.nat[0][2])
    exit(1)
# passed checks, show status
print "altering resi %s in:" % kn
for i in stored.nat:
    print "%s/%s/%s" % i
nr = aas[aal.index(stored.nat[0][2])]
# generate mutants
cmd.wizard("mutagenesis")
for i in range(len(ars)):
    if ars[i] == nr:
        continue
    cmd.get_wizard().set_mode(arl[i])
    for r in stored.nat:
        cmd.edit("/%s//%s/%s/CA" % (bn,r[1],kn))
        cmd.get_wizard().do_pick(0)
        cmd.get_wizard().apply()
    sn = "%s-%s%s%s.pdb" % (bn,nr,kn,ars[i])
    cmd.save(sn)
    print "saved %s" % sn
cmd.set_wizard();
