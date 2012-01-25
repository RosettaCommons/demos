from numpy import array
from math import sqrt
import os,shutil
import scipy.linalg as linalg

# Distance cut off between metal ion and hetero atom
DISTANCE = 2.9

class get_metal_site:

    # Requires PDB files
    # Returns object of file
    def get_obj_pdbfile(self,pdbfile):
        fl = open(pdbfile,'r')
        pdb = fl.readlines()
        fl.close()
        return pdb

    # Requires name of metal as string
    # Returns a dictionary with the coordinates
    def get_metal_ions(self,pdbfile,mt='ZN'):
        ht_lig = ['O','S','N']
        coor_zn = {}
        hetatm = {}
        for line in pdbfile:
            if line[0:6] == 'HETATM':
                tmp = line.split()
                het_id = line[13:14]
                if tmp[2] == mt:
                    tmp_id = line[22:26] 
                    coor_zn[tmp_id] = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
                elif  het_id in ht_lig:
                    nid = line[7:11]+'   '+line[13:16]
                    hetatm[nid] = array([float(line[31:39]),float(line[39:47]),float(line[47:55])])
        return coor_zn,hetatm
    

    # Getting the distance between zinc ion and hetero atoms
    # returning the atoms less then 2.7 A from the ZN as well
    # as a boolean value in order to search for protein ligands
            
    def get_hetero_around_metal(self,metal_dic,het_dic):
        # Dictionary with coordinates and names
        active_site = {}
        for i in metal_dic:
            metal_vec = metal_dic[i]
            for j in het_dic:
                # Get vector
                value =  het_dic[j]
                metaldis = linalg.norm(metal_vec - value)
                if metaldis < 2.7:
                    active_site[j]=value
        return active_site
    

    # Requires fileobject, metal coordinates, and active site coordinates
    # Returns dictionary with ligands < 2.9 from Zn,list of protein ligands
    # and their distance to Zn
    def get_protein_ligand_metal(self,pdb_obj,zn_coor,active_site):
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
                    for i in zn_coor:
                        ds_lig = linalg.norm(vec-zn_coor[i])
                        if linalg.norm(vec-zn_coor[i]) < DISTANCE:
                            iid = res+'  '+line[23:26]
                            mt_lig.append(iid+'  '+atom+'  '+str(ds_lig))
                            active_site[iid] = vec
        
        return active_site,mt_lig
