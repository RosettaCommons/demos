import sys,os

class get_atoms:

    # Requires file
    # returns file object
    def read_file(self,filename):
        fl = open(file,'r')
        return fl
    
    # Only heteroatoms
    # Requires filename and list of atoms
    # Returns atoms from pdb file
    def get_hetatoms(self,filename,atoms):
        fl = open(filename,'r')
        # List of lines
        atm = {}
        for line in fl:
            lgth = len(line)
            if lgth > 5 and line[0:6] == 'HETATM':
                i = line.split()[2]
                if i in atoms:
                    atm[i] = line
        return atm

    # Define method
    
    def get_hetatoms_pdb(self,filename,atoms):
        # Need to split into atomnr and name
        atmnr,atmname = atoms.split()
        atmnr = str(atmnr)
        fl = open(filename,'r')
        # List of lines
        for line in fl:
            lgth = len(line)
            if lgth > 5 and line[0:6] == 'HETATM':
                i = line[7:11] 
                j = line[13:15] 
                j = j.strip()
                i = i.strip()
                if i ==  atmnr and j == atmname[0:2]:
                    return line



    def get_metalion_pdb(self,filename,atoms):
        fl = open(filename,'r')
        # List of lines
        for line in fl:
            lgth = len(line)
            if lgth > 5 and line[0:6] == 'HETATM':
                i = line.split()[2]
                if i == atoms:
                    return line




    # Requires filename, resname, resid, atom
    # Returns line in pdb file
    def get_lig_atoms(self,filename,resn,resid,atom):
        fl = open(filename,'r')
        # List of lines
        atm = {}
        for line in fl:
            lgth = len(line)
            if lgth > 5 and line[0:4] == 'ATOM':
                lst = line.split()
                # changed 06-01-2010
                #if lst[5] == resid:
                if str(line[23:26]).strip() == resid:
                    #if lst[3] == resn:
                    if line[17:20] == resn:
                        # if lst[2] == atom:
                        if str(line[12:16]).strip() == atom:
                            return line                            

    # Requires list of atoms
    # Write cry file with atoms
    def write_atoms_file(self,atoms):
        wr = open('cry.pdb','w')
        for atm in atoms:
            wr.write(atm)


    def main():
        fl1 = sys.argv[1]
        atoms = []
        position = 0 
        for arg in sys.argv[2:]:
            atoms.append(str(arg))
        ln = get_atoms_file(fl1,atoms)
        write_atoms_file(atoms,ln)
