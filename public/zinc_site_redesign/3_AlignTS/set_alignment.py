import os,shutil,sys
from alignment import *

from get_metal_site import *
from get_atoms import *
from write_pdb_files import *

# Initialization of classes
no = get_metal_site()
ga = get_atoms()
wpf = write_pdb_files()

VIZ = 'END\n'

class set_alignment:


	# Requires PDB file
	# Returns tmp.pdb, superimposed.pdb, sample_substrate.pdb, cry.pdb

	def generate_alignment_all(self,METAL,PDB,LIGANDFILE,SUBSTRATEATOMS):
		HETATM = 1
		PROTEINLIG = 1
		MOREZNATOMS = 0
		NoHetAtoms = 0
		
		# Getting list of transformed coordinates - append to a list 
		atomic_coor = []
		
		# PDB file to be aligned to - here make loop over all files
		file_object = no.get_obj_pdbfile(PDB)
		
		# ZN coordinates and heteroatoms are retrived
		zn_coor,het = no.get_metal_ions(file_object)
	
		# If more then one zinc atom in structure
		if (len(zn_coor) != 1):
			print 'More zinc atoms',PDB
			MOREZNATOMS = 1
			HETATM = 0

		# Returns heteroatom around Zn atom
		active_site = no.get_hetero_around_metal(zn_coor,het)
    
		if len(active_site) == 0:
			print 'In the PDB file '+PDB+' no hetatms ligands were found'
			HETATM = 0
			NoHetAtoms = 1

		if len(active_site) != 0:
			het_atom = active_site.keys().pop()
			# Returns dictionary with active site atoms from protein and ligand plus a list with
			# resname and resid of protein ligands

		ac_site,protein_lig = no.get_protein_ligand_metal(file_object,zn_coor,active_site)

		if (HETATM == 1 and PROTEINLIG == 1):    
			# Generating list with atoms from crystal structure to align to
			# Metal set in the begining of file, het_atom is deduced from above
			lst = []
			lst.append(ga.get_metalion_pdb(PDB,METAL))
			lst.append(ga.get_hetatoms_pdb(PDB,het_atom))

			# Atoms to use for the alignment
			# The atoms are from the tetralhedral geometry
			# Write to cry.pdb file - used to solve the SVD
			ga.write_atoms_file(lst)

			# Writing the tmp.pdb file with the atoms which should be aligned
			self.get_atoms_from_ligand(SUBSTRATEATOMS,LIGANDFILE)

			# Getting atomic coordinates
			atomic_coor.append(self.get_aligned_coor(LIGANDFILE))
			
			# Write coordinates to file sample_substrate.pdb
			wpf.write_aligned_coordinates(atomic_coor,'sample_substrate.pdb')
		

	# Requires floating point number
	# Returns floating point with correct
	# number of digits for pdb
	def set_number_digits(self,number):	
		return '%.3f' %number


	# Requires list of ligand atoms  
	# Writes file 'tmp.pdb' for the alignment
	# Need to fix this 15-09-2011
	def get_atoms_from_ligand(self,atm,ligandfile):
		pdb_lines = {}
		fileobject = open(ligandfile,'r')
		# File name
		fname = 'tmp.pdb'
		tm = open(fname,'w')
		for line in fileobject:
			if len(line.split()) > 3:
				atm_name = line.split()[2]
				if atm_name in atm:
					pdb_lines[atm_name] = line
		for i in atm:
			tm.write(pdb_lines[i])
		tm.close()



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

	# Returns a list with transformed coordinates
	def get_aligned_coor(self,templateFile):
		# Initialize object
		obj = align_to_substrate()
		# Reading data from crystal structure where one wants
		# the alignment from
		# Get data crystal structure
		cry_data,atom_names = obj.get_data('cry.pdb')
		# Name of list the transformed coordinates are collected in
		outfile = []
		fname = 'tmp.pdb'
		sub_data,atom_names = obj.get_data(fname)
		# Superimpose substrate data in crystal structure
		# getting the translation and rotation matrix
		t_m, r_m = obj.get_rotate_translate(sub_data,cry_data)
		# Getting the transformed coordinates
		nw = obj.get_transformed_coor(sub_data,cry_data)
		# We transform the original data
		sub,at = obj.get_data(templateFile)
		# The transformed coordinates
		t_c = dot(sub,r_m) + t_m
		# Writing the coordinates
		# Files name of coordinates is
		# Writing to a file called superimposed.pdb
		obj.write_pdb(at,t_c)
		# File for rosetta with the correct naming
		sp_file = open('superimposed.pdb','r')
		rosetta = open(templateFile,'r')
		fileOne = sp_file.readlines()
		fileTwo = rosetta.readlines()
		rosetta.close()
		# Variable to count line number in other file
		# used to insert at the right line
		ct = 0
		for i in fileTwo:
			ln = fileOne[ct].split()
			x = self.set_number_digits(float(ln[6]))
			y = self.set_number_digits(float(ln[7]))
			z = self.set_number_digits(float(ln[8]))
			x = self.set_length_digit(x)
			y = self.set_length_digit(y)
			z = self.set_length_digit(z)		
			i = str(i[0:30])+x+y+z+str(i[55:81])
			outfile.append(i)
			ct = ct +1
		outfile.append(VIZ)
		return outfile
