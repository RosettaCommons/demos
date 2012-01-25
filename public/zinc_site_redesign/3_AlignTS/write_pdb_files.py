'''
Write and manipulate pdb files


'''
class write_pdb_files:

    
    # Requires filename
    # Returns fileobject
    def open_filename(self,filename):
        fl = open(filename,'r')
        return fl


    # Requires floating point number
    # Returns floating point with correct
    # number of digits for pdb
    def set_number_digits(self,number):	
        return '%.3f' %number

    # Requires a number
    # Set the right length for the number
    def set_length_digit(self,number):
        lngth = len(number)
        if lngth == 7:
            return ' '+number
        if lngth == 6:
            return '  '+number
        if lngth == 5:
            return '   '+number
        if lngth == 4:
            return '    '+number
        else:
            return number

        
    
    # Requires list, filename
    # write coordinates to file
    def write_coordinates(self,list_of_coordinates,filename):
        list_of_coordinates = self.merge(list_of_coordinates)
        wr = open(filename,'w')
        lngth = len(list_of_coordinates)
        while lngth > 0:
            crd = list_of_coordinates.pop(0)
            lngth = lngth - 1
            for i in crd:
                wr.write(i)

    # Requires pdb file
    # Returns protein atoms
    def get_atoms_pdb(self,filename):
        atm = []
        fl = self.open_filename(filename)
        for line in fl:
            if line[0:4] == 'ATOM':
                atm.append(line)
        return atm

    def get_metal_ion(self,filename,metalname='ZN'):
        metal = []
        fl = self.open_filename(filename)
        for line in fl:
            if line[0:4] == 'HETA':
                tmp = line.split()[2]
                if tmp == 'ZN':
                    metal.append(line)
        return metal

    def get_xyz(self,pdbline):
        splt = pdbline.split()
        # Fix 28-12-2009
        x = self.set_number_digits(float(pdbline[30:38]))
        y = self.set_number_digits(float(pdbline[38:46]))
        z = self.set_number_digits(float(pdbline[46:54]))
        x =  self.set_length_digit(x)
        y =  self.set_length_digit(y)
        z =  self.set_length_digit(z)
        return x,y,z

    # Require file
    # Input crystal coordinates 
    def set_correct_line(self,filename,coordinates):
        # File container
        fc = []
        fl = self.open_filename(filename)
        x,y,z = self.get_xyz(coordinates)
        
        for line in fl:
            if line[0:4] == 'HETA' and line[12:15] == 'ZN1':
                n_line = str(line[0:30])+x+y+z+str(line[55:])
                fc.append(n_line)
            else:
                fc.append(line)
        return fc
    
    def write_aligned_coordinates(Self,list_of_coordinates,name='aligned_ligand.pdb'):
        wr = open(name,'w')
        lngth = len(list_of_coordinates)
        while lngth > 0:
            crd = list_of_coordinates.pop(0)
            lngth = lngth - 1
            for i in crd:
                wr.write(i)

    # Requires samplefile with ligand
    # Return ligand coordinates
    def get_ligand(self,filename):
        fl = self.open_filename(filename)
        lgnd = []
        for line in fl:
            if line[0:4] == 'HETA':
                lgnd.append(line)
            else:
                break
        return lgnd

    

