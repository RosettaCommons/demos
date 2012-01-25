#!/usr/bin/env python
import sys
from optparse import OptionParser
'''
Fix_pdbfile get chain A as well as hetatms connected to chain A.
It takes all a pdb file as input and returns
the same file.

@parameter f: PDB file to be fixed
@type       : string
@parameter c: chain in pdbfile default is chain a
@type       : character
@parameter m: metal to remove from pdb
@type       : string
@parameter n: residue number of metal to remove from pdb
@type       : string

@return     : PDB file with MSE -> MET, KCX -> LYS with name PDBFILE_clean.pdb

To run script:

python cleanPDBfile.py -f PDBFILE

or to remove metal ion

python cleanPDBfile.py -f PDBFILE -m METAL -n NR


'''               

# Write to file
# Return a cleaned pdb file
def write_pdb_file(pdblist,name,chain):
    pdb_write = open(name+'_clean_'+str(chain)+'.pdb','w')
    for line in pdblist:
        pdb_write.write(line)

# Requires PDB file, chain to collect
# Returns PDBlist
def get_chain_pdbfile(fileobject,chain):
    pdblist = []
    pdbfile = open(fileobject,'r')
    for line in pdbfile:
        if line[0:4] == 'ATOM' or line[0:4] == 'HETA':
            if line[21] == str(chain) and line[17:20] != 'MSE':
                pdblist.append(line)
            if line[21] == str(chain) and line[17:20] == 'MSE':
                nw_line = 'ATOM  '+line[6:17]+'MET'+line[20:]
                pdblist.append(nw_line)
            if line[21] == str(chain) and line[17:20] == 'KCX':
                if line[13:16] not in ['CX ','OQ1','OQ1']:
                    nw_line = 'ATOM  '+line[6:17]+'LYS'+line[20:]
                    pdblist.append(nw_line)
        elif  line[0:4] != 'ATOM' or line[0:4] != 'HETA':
            pdblist.append(line)
    return pdblist

def remove_residue_from_list(pdblist,METAL,nr):
    print type(METAL),METAL,type(nr),nr
    new_list = []
    for line in pdblist:
        if  (str(line[23:26]).strip() == nr and str(line[17:20]).strip() == METAL):
            continue
        else:
            new_list.append(line)
    return new_list

def main():
    parser = OptionParser()
    parser.add_option('-f',dest='pdbfile',
                      help='Protein file converting MSE-> MET, KCX -> LYS')
    parser.add_option('-c',dest='chain',default='A',
                      help='chain to be fixed')
    parser.add_option('-m',dest='METAL',default=False,
                      help='Metal to be removed')
    parser.add_option('-n',dest='nr',default=False,
                      help='residue number of metal to be removed')

    (options,args) = parser.parse_args()
    name = str(options.pdbfile).split('.')[0]
    pdblist = get_chain_pdbfile(options.pdbfile,options.chain)
    if (options.METAL and options.nr):
        
        pdblist = remove_residue_from_list(pdblist,options.METAL,options.nr)
    write_pdb_file(pdblist,name,options.chain)

    
if __name__ == "__main__":
    main()
