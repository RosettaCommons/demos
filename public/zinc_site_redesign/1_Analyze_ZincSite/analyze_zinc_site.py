from numpy import *
from math import sqrt
import os,shutil
from optparse import OptionParser

'''
Analyze zinc site in PDB file and return information on its coordination
(number of residues), square planar or tetrahedral

@parameter f: PDB file containing metal ion
@type       : string
@parameter c: metal to be investigated eg. \'ZN\'
@type       : string
@return     : print information on metal site

The script is executed like this
python analyze_zinc_site.py -f PDBFILE -m METAL

tested on Zinc

'''

# Distance to hetero atoms around metal ion
DISTANCEHET = 2.7
DISTANCEPRT = 2.9
# Angles values used to calculate rmsd
TETRAHEDRAL = 109.5
TETRAPLANAR = 90

# Distance used to determine if site is multi-metal
METAL_METAL_DIS = 5

# Requires PDB files
# Returns list with lines from file
def get_obj_pdbfile(pdbfile):
    fl = open(pdbfile,'r')
    pdb = fl.readlines()
    fl.close()
    return pdb

# Calculate angle between two vectors
# from three points with vector2 being center
# Return angle in degrees
def get_calc_angle(vector1,vector2,vector3):    
    x = vector1 - vector2
    y = vector3 - vector2
    dot_p = dot(x,y)
    norm_x = sqrt((x*x).sum())
    norm_y = sqrt((y*y).sum())
    cos_angle = dot_p / norm_x / norm_y # cosinus of angle between x and y
    angle = arccos(cos_angle)*57.3
    return angle


# Calculate distance between metal ions in protein
# Requires dictionary with metal ion
# Return distances
def get_metal_metal_distance(metal_dic):
    for i in metal_dic:
        for j in metal_dic:
            if i != j:
                metal_length = linalg.norm(metal_dic[i]-metal_dic[j])
                if metal_length < METAL_METAL_DIS:
                    print 'Metal-metal distance less than '+str(METAL_METAL_DIS)+'  ',metal_length,i,j

# Requires name of metal as string
# Returns a dictionary with the coordinates
def get_metal_ions(pdbfile,mt='ZN'):
    ht_lig = ['O','S','N']
    chains = []
    coor_zn = {}
    hetatm = {}
    for line in pdbfile:
        if line[0:6] == 'HETATM':
            tmp = line.split()
            het_id = line[13:14]
            if tmp[2] == mt:
                # Add chain ID
                chainID = line[21:22]
                chains.append(chainID)
                tmp_id = line[22:26]+' '+chainID 
                coor_zn[tmp_id] = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
            elif  het_id in ht_lig:
                nid = line[0:4]+'   '+line[7:11]+'   '+line[13:16]
                hetatm[nid] = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
    return coor_zn ,hetatm,chains
    

def get_hetero_around_metal(metal_vec,het_dic):
    # Dictionary with coordinates and names
    active_site = {}
    for j in het_dic:
        value =  het_dic[j]
        metaldis = linalg.norm(metal_vec - value)
        if metaldis < DISTANCEHET:
            active_site[j]=value
    return active_site

# Requires fileobject, metal coordinates, and active site coordinates
# Returns dictionary with ligands within treshold
def get_protein_ligand_metal(pdb_obj,metal_vec,active_site):
    # Distance between ligands and metal ions
    mt_lig = []
    coor_residues = ['CYS','HIS','ASP','GLU','SER','THR']
    lig_atm = ['ND1', 'NE2','SG ', 'OE1', 'OE2','OD1','OD2','OG ']
    for line in pdb_obj:
        res = line[17:20]
        atm = line[0:4]
        if res in coor_residues and atm =='ATOM':
            atom = line[13:16]
            if atom in lig_atm:
                vec = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
                ds_lig = linalg.norm(vec-metal_vec)
                if ds_lig < DISTANCEPRT:
                    iid = res+'  '+line[23:26]
                    mt_lig.append(iid+'  '+atom+'  '+str(ds_lig))
                    active_site[iid] = vec        
    return active_site,mt_lig

# Determine distribution of coordinated residues
# Dictionary with ligands
# Prints number of ligands from protein and other
def get_type_of_coordination(ligand_dictionary):
    het = 0
    prt = 0
    for tp in ligand_dictionary:
        if tp.split()[0] == 'HETA':
            het = het + 1
        else:
            prt = prt + 1
    print 'Number of interactions between metal and protein: '+str(prt)
    print 'Number of interactions between hetero atom : '+str(het)


# Requires metal site and active site
# Returns rmsd of tetrahedral geometry
def determine_coordination(metal_vec,active_site):
    coordination = 'SQUAREPLANAR'
    th = 0
    tp = 0
    dummy = 0
    a_s = active_site.values()
    las = len(a_s)
    i = 0
    while i < las:
        coor1 = a_s[i]
        j = 1+i
        while j < las:
            coor2 = a_s[j]
            ang = (get_calc_angle(coor1,metal_vec,coor2))
            th = th + (ang - TETRAHEDRAL)**2
            tp = tp + (ang - TETRAPLANAR)**2
            dummy = dummy +1
            j=j+1
        i = i+1
    if sqrt(th)/dummy < sqrt(tp)/dummy:
        coordination = 'TETRAHEDRAL'
    return coordination


def get_number_of_ligand(ligands):
    METALSITE = 0
    # less than four ligands - remove from set
    length_ligands = len(ligands)
    if length_ligands < 4:
        METALSITE = 0
    elif length_ligands == 4:
        METALSITE = 'TH'
    elif length_ligands == 5:
        METALSITE = 'TBP'
    elif length_ligands == 6:
        METALSITE = 0
    return METALSITE
                
def main():
    parser = OptionParser()
    parser.add_option('-f',dest='pdbfile',
                      help='Charactize the metal site binding site in this pdb')
    parser.add_option('-m',dest='METAL',default='ZN',
                      help='metal to be characterized')
    (options,args) = parser.parse_args()
    
    PDB = get_obj_pdbfile(options.pdbfile)

    zinc_site,pot_ligands,chains = get_metal_ions(PDB)

    # Number of unique chains
    nr_chains = set(chains)
    # Print number of chain in PDB file
    print 'Number of chains in pdb file: ', len(nr_chains)
    
    # Test for multi zinc active site
    if len(zinc_site) > 1:
        get_metal_metal_distance(zinc_site)

    for i in zinc_site:
        active_site = {}
        # Get any heteroatom around metal site
        active_site = get_hetero_around_metal(zinc_site[i],pot_ligands)
        # containing only cysteines
        active_site,tmp = get_protein_ligand_metal(PDB,zinc_site[i],active_site)

        # Get distribution of protein ligand and heteroatoms
        get_type_of_coordination(active_site)
        
        # Check if it is cystein zinc site
        CYSTEINESITE = 0
        for j in tmp:
            if j.split()[0] == 'CYS':
                CYSTEINESITE = CYSTEINESITE +1
        # Determine type of metal site
        if(CYSTEINESITE !=4):
            GEOM = get_number_of_ligand(active_site)
            if GEOM == 0:
                print 'PDB does not contain 4 or 5 coordinated metal site'
            elif GEOM == 'TH':
                crd = determine_coordination(zinc_site[i],active_site)
                print 'Coordination type of metal site: ',crd+'  '+i
            elif GEOM == 'TBP':
                print 'PDB chain '+str(i).split()[1]+' contains 5 coordinated metal site'
        else:
            print 'PDB contains cysteine site: ',i
    
if __name__ == "__main__":
    main()
