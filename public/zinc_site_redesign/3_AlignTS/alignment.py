from numpy import *

'''
Superposition of two matrices using the algorithm
develpoed by Kabsch and Wolfgang

'''

class align_to_substrate:

    # Transform coordinates into numpy array
    # Requires coordinates in pdb format
    # Returns numpy array with coordinates
    def get_data(self,filename):
        # Need to add try/catch
        fname = open(filename,'r')
        # Coordinates
        dt = []
        # atom names
        a_n = []
        # Looping over line in file
        for line in fname:
            sl = line.split()
            # Bug fix 28-12-2009
            x = float(line[31:38])
            y = float(line[38:46])
            z = float(line[46:54])
            tmp = x,y,z  
            dt.append(tmp)
            a_n.append(sl[2])
        fname.close()
        return array(dt),a_n

    # Writes the transformed coordinates to file
    # Requires coordinates translated and rotated
    # Returns xyz file of rotated and translated
    def write_data(self,atomnames,coordinates):
        al = len(atomnames)
        #    assert al == len(coordinates),'Lengths of lists do not match'
        wname = open('superimposed.xyz','w')
        wname.write(str(al)+'\n Generated to Enzyme Design\n')
        for i in range(al):
            wname.write(atomnames[i]+'\t'+str(coordinates[i][0])+'   '+str(coordinates[i][1])+'   '+str(coordinates[i][2])+'\n')

    def write_pdb(self,atomnames,coordinates,fname='superimposed.pdb'):
        '''Write a pdb file of the aligned coordinates '''
        al = len(atomnames)
        # Name of chain
        ch = '  SUB'
        wname = open(fname,'w')
        for i in range(al):
            wname.write('ATOM   '+str(i).rjust(3)+'   '+atomnames[i]+' '+ch+' X '+str(i).rjust(3)+'    '+str(coordinates[i][0])[0:6]+'  '+str(coordinates[i][1])[0:6]+'  '+str(coordinates[i][2])[0:6]+'\n')
            
    # Kabsch Algorithm for superimposing
    # Superimpose matrix1 on matrix2
    def get_rotate_translate(self,matrix1,matrix2):
        assert shape(matrix1) == shape(matrix2), 'Matrices not of same dimensions'
        
        # Store number of rows
        nrows = shape(matrix1)[0]
        
        # Getting centroid position for each selection
        avg_pos1 = matrix1.sum(axis=0)/nrows
        avg_pos2 = matrix2.sum(axis=0)/nrows

        # Translation of matrices
        avg_matrix1 = matrix1-avg_pos1
        avg_matrix2 = matrix2-avg_pos2

        # Covariance matrix
        covar = dot(avg_matrix1.T,avg_matrix2)
        
        # Do the SVD in order to get rotation matrix
        u,s,wt = linalg.svd(covar)
        
        # Rotation matrix
        # Transposition of u,wt
        rot_matrix = wt.T*u.T
        
        # Insure a right-handed coordinate system
        # need to fix this!!
        # if linalg.det(rot_matrix) > 0:
        wt[2] = -wt[2]
        rot_matrix = transpose(dot(transpose(wt),transpose(u)))

        trans_matrix = avg_pos2-dot(avg_pos1,rot_matrix)
        return trans_matrix, rot_matrix
        

    # Returning the superimposed coordinates
    # Requires two arrays same length
    # Returns array of transformed coordinates
    def get_transformed_coor(self,a,b):
        tr, rt = self.get_rotate_translate(a,b)
        nw_coor = dot(a,rt) + tr
        return nw_coor
