#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#
import os

slim_flags = " -chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation SpecialRotamer protein_cutpoint_upper protein_cutpoint_lower VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm "

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


def touch(fname, times=None):
    fhandle = open(fname, 'a')
    try:
        os.utime(fname, times)
    finally:
        fhandle.close()

def read_config_file( config_file ):
    import ConfigParser
    assert os.path.exists( config_file ), "config file doesn't exist."
    config = ConfigParser.RawConfigParser()
    config.read( config_file )

    return config

def read_silent_header( silentfile ):
    from os.path import exists
    assert exists( silentfile )
    fields_Dict = {}

    file = open( silentfile, "r" )
    score_terms = file.readline()
    if score_terms.find('desc') < 0:
        score_terms = file.readline()

    if score_terms.find('desc') >= 0:
        #Typical out file or score file
        labels=score_terms.split()
        i=0
        while i < len(labels) :
            fields_Dict[ labels[i] ] = i
            i = i +1

    file.close()
    return fields_Dict


def create_xyzDict( pdb ):
    ''' this one is obsolete; should call create_xyzDict_bychain instead
    '''
    #from os import popen
    #coords = []
    #lines = popen('grep ^ATOM %s | grep " %s "' %( pdb, atom )).readlines()
    #lines = popen('grep ^ATOM %s ' %( pdb )).readlines()
    #for line in lines:
    #    coords.append([line])
    #items_num  = len(coords)

    #dist = 0
    xyzDict_    = {}  # { rsn: { atom:xyz }}
    chainDict_  = {}  # { rsn: chain_id }
    pdblineDict = {}  # { rsn: pdbline }
    resDict     = {}  # { rsn: residue_name }
    #for pdbline in popen('grep ^ATOM %s ' %( pdb )).readlines():
    with open( pdb, "r" ) as f:
        for pdbline in f:
            if pdbline.startswith("ATOM"):
                rsn      = int( pdbline[22:26]) # res number
                res      = pdbline[17:20].strip() # res name
                atom     = pdbline[12:16].strip() # atom name
                chain_id = pdbline[21] # chain id
                xyz      = map(float, [ pdbline[30:38], pdbline[38:46], pdbline[46:54] ])
                if rsn not in xyzDict_.keys():
                    xyzDict_[ rsn ]    = { atom:xyz }
                    chainDict_[ rsn ]  = chain_id
                    pdblineDict[ rsn ] = pdbline
                    resDict[ rsn ]     = res
                else:
                    xyzDict_[ rsn ][ atom ]  = xyz
                    pdblineDict[ rsn ] += pdbline

    return ( xyzDict_, chainDict_, pdblineDict, resDict )



def fasta_file_reader( fasta_fn ):
    from Bio import SeqIO
    fasta_seq = str( SeqIO.read( fasta_fn, "fasta" ).seq )

    return fasta_seq
