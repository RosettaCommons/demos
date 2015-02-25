#! /usr/bin/env python

import numpy
import operator

dict_3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

dict_1to3 = dict((v,k) for k, v in dict_3to1.iteritems())

def get_3l_code( letter ):
	return dict_1to3[ letter.upper() ]

def get_1l_code( letters ):
	return dict_3to1[ letters.upper() ]

def parse_fasta( file ):
	        f = open ( file, 'r' )
		sequence=""
		grep_sequence = False
		for line in f:
				if ( grep_sequence and line.startswith(">") ):
						grep_sequence  = False
						break
				if ( grep_sequence and not line.startswith(">") ):
						sequence += line.strip()
				if not grep_sequence:
						grep_sequence = True
		f.close()
		return sequence

#parses cmpfile and returns a numpy array
def parse_cmpfile( file ):
	f = open (file, 'r')
	count = 0
	for line in f:
		split = line.split()
		if (count == 0):
			ma = numpy.zeros((len(split), len(split)), numpy.float)
		for i in range (0, len(split)):
			ma[count, i] = numpy.float(split[i])
		count += 1
	f.close()
	return ma

#parse rigid file and return residues in rigid file
def parse_rigid_file( file ):
	converged_residues = []
	f = open( file, "r" )
	for line in f:
		if len( line.strip() ) > 0 and line.split()[0] == "RIGID":
			converged_residues += range( int(line.split()[1]), int(line.split()[2])+1 )
	f.close()
	return converged_residues
	
	
# extract all restraints from matrix, sorted by value
def extract_cm_restraints( matrix, min_sep, one_indexed=True ):
	restraints = []
	for i in range(matrix[0].size):
		for j in range(i+min_sep, matrix[0].size):
			if one_indexed and matrix[i,j]>0:
				restraints.append((i+1,j+1,matrix[i,j]))
			elif matrix[i,j]>0:
				restraints.append((i, j, matrix[i,j]))
	return sorted(restraints, key=operator.itemgetter(2), reverse = True) # Sort according to Value stored in ContactMap


def print_restraints( constraints, lb, ub, atom, function, weight, seq="", use_weights = False ):
	restraints = ""
	a1 = a2 = atom
	for ele in constraints:
		a1 = a2 = atom
		if seq != "" and atom == "CB":
			if seq[ele[0]-1:ele[0]] == "G":
				a1 = "CA"
			if seq[ele[1]-1:ele[1]] == "G":
				a2 = "CA"
		if use_weights:
			we = ele[2] if ele[2] > 0 else 0.0001
			weights = "SCALARWEIGHTEDFUNC {0:.4f}".format( we )
		else:
			weights = ""
		if ele[2] > 0:
			if function == "BOUNDED":
				restraints += "AtomPair {0:>4} {1:>4} {2:>4} {3:>4} {8} BOUNDED {4:.2f} {5:.2f} {6:.2f} #ContactMap: {7:.2f}\n".format(a1, ele[0], a2, ele[1], lb, ub, weight, ele[2], weights )
			elif function == "SIGMOID":
				restraints += "AtomPair {0:>4} {1:>4} {2:>4} {3:>4} {7} SIGMOID {4:.2f} {5:.2f} #ContactMap: {6:.2f}\n".format(a1, ele[0], a2, ele[1], ub, weight, ele[2], weights )
			else:
				print "Function {} not implemented!".format( function )
				return ""
	return restraints
	
def print_padded_restraints( constraints, padding, atom, function, weight, seq=""  ):

	restraints = ""
	for ele in constraints:
		lb = ele[2] - padding
		ub = ele[2] + padding
		a1 = a2 = atom
		if seq != "" and atom == "CB":
			if seq[ele[0]-1:ele[0]] == "G":
				a1 = "CA"
			if seq[ele[1]-1:ele[1]] == "G":
				a2 = "CA"
		if function == "BOUNDED":
			restraints += "AtomPair {0:>4} {1:>4} {7:>4} {2:>4} BOUNDED {3:.2f} {4:.2f} {5:.2f} #ContactMap: {6:.2f}\n".format(a1, ele[0], ele[1], lb, ub, weight, ele[2], a2 )
		elif function == "SIGMOID":
			restraints += "AtomPair {0:>4} {1:>4} {6:>4} {2:>4} SIGMOID {3:.2f} {4:.2f} #ContactMap: {5:.2f}\n".format(a1, ele[0], ele[1], ub, weight, ele[2], a2 )
		else:
			print "Function {} not implemented!".format( function )
			return ""
	return restraints
	
def get_pdb_sequence( pdb_file, chain ):

	models = []
	seq = ""

	pre_res = -100000


	f = open( pdb_file, "r" )
	for line in f:
		if line.startswith("ATOM"):
			if line[21:22] == chain or chain =="*":
				cur_res = int(line[22:26].strip())
				if pre_res == -100000:
					pre_res = cur_res -1
				if cur_res == pre_res +1:
					aa = dict_3to1[line[17:20].strip()]
					seq += aa
					pre_res = cur_res
		elif line.startswith(("TER","END")):
			if len(seq) != 0:
				models.append( seq )
			pre_res = -100000
			seq = ""

	if len(seq ) != 0:
		models.append( seq )

	f.close()
	return models
	
def get_all_coordinates( pdb_file, chain = "*" ):

	models = []
	residues = []
	pre_res = -1
	atoms = {}

	f = open( pdb_file, "r" )
	for line in f:
		if line.startswith("ATOM"):
			if line[21:22] == chain or chain == "*":
				cur_res = int(line[22:26].strip())
				if pre_res == -1:
					pre_res = cur_res -1
				if cur_res == pre_res +1:
					if line[12] != "H":
						atoms[line[12:16].strip()] = (float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()))
				else:
					residues.append(atoms)
					atoms = {}
					if line[12] != "H":
						atoms[line[12:16].strip()] = (float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()))
					pre_res = cur_res - 1
		elif line.startswith(("TER","END")):
			if len(atoms.keys()) != 0:
				residues.append(atoms)
			if len(residues) !=0:
				models.append( residues )
				residues = []
				atoms={}
			pre_res = 0

	if len(atoms.keys()) != 0:
		residues.append(atoms)
	if len(residues) !=0:
		models.append( residues )
		residues = []


	f.close()
	return models
	
def calculate_dist_matrix_from_coords( coords, atom):

	matrix = numpy.zeros((len(coords), len(coords)), numpy.float)

	for i in range( 0,len(coords) ):
		for j in range( i+1, len(coords) ):

			a1 = a2 = atom

			if not a1 in coords[i].keys():
				print "Using atom CA at position {0}".format(i+1)
				a1 = "CA"
			if not a2 in coords[j].keys():
				print "Using atom CA at position {0}".format(j+1)
				a2 = "CA"

			dist = get_distance( coords[i][a1], coords[j][a2])
			matrix[i,j] = dist
			matrix[j,i] = dist

	return matrix
	
def get_distance( coord1 , coord2 ):
	sum=0
	for i in range(0, len(coord1) ):
		sum += ( coord1[i] - coord2[i] ) * ( coord1[i] - coord2[i] )
	return numpy.sqrt( sum )

def calc_dist_matrix_from_coords( coordinates ):
	matrix = numpy.zeros((len(coordinates), len(coordinates)), numpy.float)

	for i in range( 0,len(coordinates) ):
		for j in range( i+1, len(coordinates) ):
			dist = get_distance( coordinates[i], coordinates[j])
			matrix[i,j] = dist
			matrix[j,i] = dist

	return matrix


