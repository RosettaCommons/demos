#!/usr/bin/python
from optparse import OptionParser
from math import sqrt
from numpy import *
import os,shutil

'''
Fix_pdbfile get chain A as well as hetatms connected to chain A.
It takes all a pdb file as input and returns
the same file.

@parameter f: PDB file with zinc site 
@type       : string
@parameter m: Metal ion - default ZN 
@type       : string
@parameter n: Name of rosetta output pdb file 
@type       : string
@parameter l: Three ligand names with Zn and its first ligand
@type       : string
@parameter a: Coordinate file of aligned ligand and is append to the Rosetta input
@type       : string
@return     : constraint file, pdb file to use as input for rosetta

'''
# Parameters for the distance geometry interactions
# Distance between metal and protein
DISTANCEMETAL = 2.9
# Distance between metal and hetero atom
DISTANCEHET = 2.7

# Requires PDB files
# Returns list with lines in file
def get_obj_pdbfile(pdbfile):
    fl = open(pdbfile,'r')
    pdb = fl.readlines()
    fl.close()
    return pdb

# Requires list
# Returns one list
def merge(self,seq):
    merged = []
    for s in seq:
        for x in s:
            merged.append(x)
    return merged

# Calculate angle between two vectors
# from three points with vector2 being center
# Return angle in degrees
def get_calc_angle(vector1,vector2,vector3):    
    x = vector1 - vector2
    y = vector3 - vector2
    dot_p = dot(x,y)
    norm_x = sqrt((x*x).sum())
    norm_y = sqrt((y*y).sum())
    # cosinus of angle between x and y
    cos_angle = dot_p / norm_x / norm_y
    angle = arccos(cos_angle)*57.3
    return angle

# Return angle between two vectors in degrees
def get_angle(n1,n2):
    dot_p = dot(n1,n2)
    norm_x = sqrt((n1*n1).sum())
    norm_y = sqrt((n2*n2).sum())
    cos_angle = dot_p / norm_x / norm_y # cosinus of angle between x and y
    angle = arccos(cos_angle)*57.3
    return angle

# Requires 4 vectors (numpy.array([]))
# Generates plane between four vectors
# Return angle between the two planes
def get_dihedral_angle(v1,v2,v3,v4):
    v12 = v2-v1
    v23 = v3-v2
    v34 = v4-v3
    # Normal vectors are generated
    n1 = cross(v12,v23)
    n2 = cross(v23,v34)

    angle = get_angle(n1,n2)
    # Getting the sign of the angle
    sign = dot(n1,v34)
    if sign < 0:
        angle = 360 - angle
    return angle

# Requires pdb file
# Returns protein atoms
def get_atoms_pdb(PDB):
    atm = []
    for line in PDB:
        if line[0:4] == 'ATOM':
            atm.append(line)
    return atm

# Requires pdb file
# Returns heteroatoms
def get_heteroatoms_pdb(PDB):
    atm = []
    for line in PDB:
        if line[0:4] == 'HETA':
            atm.append(line)
    return atm


# Requires name of metal as string
# Only for one metal site
# Returns a list with metal coordinates
def get_metal_ion(pdbfile,METAL='ZN'):
    ht_lig = ['O','S','N']
    coor_zn = []
    hetatm = {}
    for line in pdbfile:
        if line[0:6] == 'HETATM':
            atom_name = str(line[12:15]).strip()
            het_id = line[13:14]
            if atom_name == METAL:
                coor_zn.append(array([float(line[31:39]),float(line[39:47]),float(line[47:55])]))
            elif  het_id in ht_lig:
                nid = line[7:11]+'   '+line[13:16]
                hetatm[nid] = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
    return coor_zn ,hetatm
    
# Requires metal coordinates, dictionary with heteroatoms
# Returns heteroatoms as vector
def get_hetero_around_metal(metal_vec,het_dic):
    for j in het_dic:
        value =  het_dic[j]
        metaldis = linalg.norm(metal_vec - value)
        if metaldis < DISTANCEHET:
            return value

# Requires residue name of protein ligand
# Returns names of residue required to determine
# geometry (angle, torsion etc.)
def get_list_of_atoms(resn):        
    atms = {
        'HIS' : ['NE2','CE1','ND1'],
        'HIS2' : ['ND1','CE1','NE2'],
        'GLU' : ['OE1','CG','CD'],
        'GLU2' : ['OE2','CG','CD'],
        'ASP' : ['OD1', 'CG', 'CB'],
        'ASP2' : ['OD2', 'CG', 'CB'],
        'CYS' : ['SG','CB','CA'],
        'LYS' : ['NZ','CE','CD'],
        'THR' : ['OG1','CB','CA'],
        'SER' : ['OG','CB','CA'],
        'ARG' : ['NH1','CZ','NE'],
        'ARG2': ['NH2','CZ','NE'],
        'ARG3': ['NE','CZ','CD'],
        'GLN' : ['OE1','CD','CG'],
        'ASN' : ['OD1','CG','CB'],
        'TYR' : ['OH','CZ','CE2']
    }
    return atms[resn]


# Get residue from pdb
# returns dictionary with vector necessary for
# generating constraints.
def get_residue_constraint_pdb(PDB,resid, resn):
    residue = {}
    # list of atoms 
    atoms = get_list_of_atoms(resn)
    for line in PDB:
        res = line[17:20].rstrip()
        atm = line[0:4]
        if atm =='ATOM':
            resnr = line[22:26].strip()
            if res == resn:
                if resnr == resid:
                    atom = line[13:16].rstrip()
                    if atom in atoms:
                        vec = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
                        residue[atom] = vec                        
    return residue,atoms

# Requires list with pdb lines, metal coordinates, and active site coordinates
# Returns dictionary with ligands < DISTANCE from Zn,list of protein ligands
# and their distance to Zn
def get_protein_ligand_metal(PDB,metal_vec):
    TRUE = 0
    # Distance between oxygens of phosphate and aa sidechains atoms
    # Distance between ligands and metal ions
    mt_lig = []
    active_site = {}
    coor_residues = ['THR','SER','LYS','ARG','ASN','HIS','GLN','TYR','ASP','GLU']
    lig_atm = ['NE','NH2','NH1','NZ','OG','OG1','ND1','ND2', 'NE2','SG ', 'OE1', 'OE2','OD1','OD2','OH']
    for line in PDB:
        res = line[17:20].rstrip()
        atm = line[0:4]
        if res in coor_residues and atm =='ATOM':
            atom = line[13:16].rstrip()
            if atom in lig_atm:
                vec = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
                ds_lig = linalg.norm(vec-metal_vec)
                if ds_lig < DISTANCEMETAL:
                    # assign value for ASP/GLU etc
                    if atom == 'OE2' and res == 'GLU':
                        res = 'GLU2'
                    elif atom == 'OD2' and res == 'ASP':
                        res = 'ASP2'
                    elif atom == 'ND1' and res == 'HIS':
                        res = 'HIS2'
                    iid = res+'  '+line[22:26]
                    mt_lig.append(iid)
                    chain = line[21:22]
    return mt_lig,chain

# Map of amino acids necessary for
# Rosetta internal naming
# Requires residue name
# Return rosetta 1-letter code
def get_single_residue(resn):
    res_map = {
        'CYS' : 'C',
        'HIS' : 'H',
        'HIS2': 'H',
        'ASP' : 'D',
        'ASP2' : 'D',
        'GLU' : 'E',
        'GLU2' : 'E',
        'LYS' : 'K',
        'SER' : 'S',
        'ARG' : 'R',
        'ARG2' : 'R',
        'ARG3' : 'R',
        'THR' : 'T',
        'ASN' : 'N',
        'GLN' : 'Q',
        'TYR' : 'Y'
        }
    return res_map[resn]


# Map of atom names necessary for
# Rosetta internal naming
# Requires atom name
# Return rosetta atom naming
def get_single_atom(atom):
    Atm_map = {
        'SG' : 'S',
        'NE2' : 'Ntrp',
        'ND1' : 'Nhis',
        'OD1' : 'OOC',
        'OD2' : 'OOC',
        'OE1' : 'OOC',
        'OE2' : 'OOC',
        'NZ'  : 'Nlys',
        'OG1' : 'OH',
        'OG'  : 'OH',
        'NH1' : 'Narg',
        'NH2' : 'Narg',
        'NE'  : 'Narg',
        'OH'  : 'OH'
    }
    return Atm_map[atom]

# Write constraint file with torsion A set to default value
# Requires atom names, residue names, distance etc, ligand names
# Returns string with parameters in Rosetta format
def write_constraint_file(resn,resi,atom,phosphate_atom,disAB,angA,angB,torAB,torB,LIGANDNAMES,LIGANDRESNAME):
    atm_type     = get_single_atom(atom)
    residue_type = get_single_residue(resn)
    const = '''

    # ZN - '''+str(resn)+' '+str(resi)+'\n'+'''
    CST::BEGIN
    TEMPLATE::   ATOM_MAP: 1 atom_name:  '''+LIGANDNAMES+'''
    TEMPLATE::   ATOM_MAP: 1 residue3:  '''+LIGANDRESNAME+'''
    
    TEMPLATE::   ATOM_MAP: 2 atom_type: '''+str(atm_type)+''' ,
    TEMPLATE::   ATOM_MAP: 2 residue1:  '''+str(residue_type)+'''
    
    CONSTRAINT:: distanceAB:  '''+str(disAB)+'''  0.20  100.  1
    CONSTRAINT::    angle_A:  '''+str(angA)+'''   25.0  10.0  180.  
    CONSTRAINT::    angle_B:  '''+str(angB)+'''   10.0  10.0  180.  
    
    CONSTRAINT::  torsion_A:   60.0    30.0  1.0  120.  
    CONSTRAINT:: torsion_AB: '''+str(torAB)+'''   10.0  10.00  360.  
    CONSTRAINT::  torsion_B:  '''+str(torB)+'''   10.0  10.00  360.  

    CST::END'''
        
    return const
# Computes geometry
# Return gemetrical values
def get_geometry(atoms,aminos,metal_site):
    l_a = atoms[0]
    l_s =atoms[1]
    l_t = atoms[2]
    dis =  '%.2f' %(linalg.norm(metal_site[0]-aminos[l_a]))
    angA = '%.2f' %(get_calc_angle(metal_site[1],metal_site[0],aminos[l_a]))
    angB = '%.2f' %(get_calc_angle(metal_site[0],aminos[l_a],aminos[l_s]))
    torAB = '%.2f' %(get_dihedral_angle(metal_site[1],metal_site[0],aminos[l_a],aminos[l_s]))
    torB  = '%.2f' %(get_dihedral_angle(metal_site[0],aminos[l_a],aminos[l_s],aminos[l_t]))
    return dis,angA,angB,torAB,torB


# Requires PDB
# Returns chain of first pdb
def get_chain(pdb):
    for line in pdb:
        if line[0:4] == 'ATOM':
            chain = line[21:22]
            break
    return str(chain)


# Requires residue name, type, number and chain
# Returns string with constraint remarks
def set_remarks_pdb(resn,resid,number,chain):
    
    strng = 'REMARK   0 BONE TEMPLATE X LG1    0 MATCH MOTIF '+chain+' '+str(resn)+'    '+str(resid)+'   '+str(number)+'\n'
    return strng

# Get geometrical constraint, write to file with them
# Return remarks necessary for Rosetta PDB file
def get_rosetta_constraint_files(PDB,PDBNAME,METAL,LIGANDNAMES,LIGANDRESNAME):
    remark = []
    dummy = 1
    
    rosetta_cst = open('constraint.cst','w')

    metal_site,het_dic = get_metal_ion(PDB,METAL)
    
    metal_site.append(get_hetero_around_metal(metal_site[0],het_dic))


    ligand,CHAIN = get_protein_ligand_metal(PDB,metal_site[0])
    for lig in ligand:
        
        tm = lig.split()
        resn =tm[0].rstrip()
        resid =tm[1].rstrip()


        aminos,atoms = get_residue_constraint_pdb(PDB,resid, resn)
        # Define what the different values mean
        # Fix torsion
        a,b,c,d,e = get_geometry(atoms,aminos,metal_site)

        # Adding modification
        tmp_constraint = write_constraint_file(resn,resid,atoms[0],METAL,a,b,c,d,e,LIGANDNAMES,LIGANDRESNAME)
        rosetta_cst.write(tmp_constraint)
        # Is this necessary still?
        if resn == 'ARG2' or resn == 'ARG3':
            resn = 'ARG'
        elif resn == 'HIS2':
            resn = 'HIS'
        
        chain = get_chain(PDB)
        remark.append(set_remarks_pdb(resn,resid,dummy,CHAIN))
        dummy = dummy +1
    return remark
                                


def main():
    parser = OptionParser()
    parser.add_option('-f',dest='PDB',
                      help='Cleaned pdb file with metal ion present default=ZN')
    parser.add_option('-m',dest='METAL',default='ZN',
                      help='Metal ion to get constraints for')
    parser.add_option('-n',dest='PDBNAME',default='cst.pdb',
                      help='Name for output pdb file as input for Rosetta')
    parser.add_option('-l',dest='LIGANDNAMES',default='ZN1 O2 P1',
                      help='Names from ligands involved in constraints e.g. ZN1 O2 C2')
    parser.add_option('-t',dest='LIGANDRESNAMES',default='LG1',
                      help='Residue name of ligand default=LG1')
    parser.add_option('-a',dest='LIGANDCOOR',default=False,
                      help='PDB coordinates of ligand appending it to rosetta pdb input')


    (options,args) = parser.parse_args()
    
    PDB = get_obj_pdbfile(options.PDB)
    PDBNAME = options.PDBNAME
    METAL = options.METAL
    LIGANDNAMES = options.LIGANDNAMES
    LIGANDRESNAMES = options.LIGANDRESNAMES
    
    remark_pdb = open('rosetta_'+PDBNAME,'w')
    remark = get_rosetta_constraint_files(PDB,PDBNAME,METAL,LIGANDNAMES,LIGANDRESNAMES)

    tmp_pdb = get_atoms_pdb(PDB)
    pdb_file = remark + tmp_pdb
    
    if(options.LIGANDCOOR):
        lig_file = get_obj_pdbfile(options.LIGANDCOOR)
        lig_file = get_heteroatoms_pdb(lig_file)
        pdb_file = pdb_file + lig_file

       
    for line in pdb_file:
        remark_pdb.write(line)
        


if __name__ == "__main__":
    main()
