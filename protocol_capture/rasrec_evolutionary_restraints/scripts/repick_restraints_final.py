#! /usr/bin/env python

import argparse
import Utils as utils
import numpy as np
from copy import copy




def main():
	parser = argparse.ArgumentParser(description='extract restraints from contactmap and top scoring decoys (pdb)', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument( '-c', '--cmp', help='contact map', required=True )
	parser.add_argument( '-s', '--structures', help='top scoring structures', required=True, nargs="+" )
	parser.add_argument( '-ch', '--chain', help="protein chain", default = "A")
	parser.add_argument( '-a', '--atom', help=' atoms used for distance restraints [CA/CB]', default = "CA")
	parser.add_argument( '-m', '--min_dist', help='minimum sequence seperation for distance restraints', default = 5)
	parser.add_argument( '-f', '--function', help='function for the restraints [BOUNDED/SIGMOID]', default = "BOUNDED")
	parser.add_argument( '-lb', '--lower_bound', help= 'lower bound used for restraints', default = 1.5, type = float )
	parser.add_argument( '-ub', '--upper_bound', help= 'upper bound/x0 used for restraints', default = 8, type = float )
	parser.add_argument( '-p', '--padding', help="use padding instead of mixed boundaries", default = -1, type = float )
	parser.add_argument( '-sd', '--standard_deviation', help="standard deviation used for getting converged distances", default = 1, type = float)
	parser.add_argument( '-o', '--outfile', help ="file name for restraint files", default = "repicked_restraints")
	parser.add_argument( '-pa', '--parameter', help = 'parameter for restraint function [slope, sd]', default = 1, type = float )
	parser.add_argument( '-nr', '--number_res', help = "number of top scoring residue pairs that gets used for repicking", default = 2000, type=int )
	args = parser.parse_args()

	# parse contact map
	cmp = utils.parse_cmpfile( args.cmp )

	# parse pdb structures
	distance_matrices = []
	for pdb_file in args.structures:
		seq = utils.get_pdb_sequence( pdb_file , args.chain )[0]
		coordinates = utils.get_all_coordinates( pdb_file , args.chain )
		for coords in coordinates:
			distance_ma = utils.calculate_dist_matrix_from_coords( coords , args.atom )
			distance_matrices.append( distance_ma )
	#	for model in range(0, utils.get_model_count( pdb_file )):
	#		distance_ma = utils.calc_dist_matrix(pdb_file, args.chain, args.atom, model)
	#		distance_matrices.append( distance_ma )
	if len(cmp) != len(distance_matrices[0]):
		print "Structures do not match CMP file."
		exit(-1)
	else:
		size_matrix = len(cmp)


	#stack matrices and create mean and std array
	multi_array = np.array( distance_matrices )
	mean_array = np.zeros(( size_matrix, size_matrix ),np.float )
	std_array = np.zeros(( size_matrix, size_matrix ),np.float )
	for i in range( size_matrix ):
		for j in range( size_matrix ):
			eles = list( multi_array[:,i,j] )
			mean_array[i,j] = np.mean( eles )
			std_array[i,j] =  np.std( eles )

	# extract restraints
	# atom distances that are similar in top restraints
	lowstd_sr = []
	lowstd_lr = []
	for i in range(0, size_matrix-1):
		for j in range(i+args.min_dist, size_matrix):
			if std_array[i,j] < args.standard_deviation and mean_array[i,j] > 0:
				if mean_array[i,j] <= 8:
					lowstd_sr.append((i+1,j+1, mean_array[i,j]))
				else:
					lowstd_lr.append((i+1,j+1, mean_array[i,j]))



	cm_restraints = utils.extract_cm_restraints( cmp , 4 )

	filtered_restraints = []

	# remove_restraints that are already on other lists ....
	for restraint in cm_restraints[:args.number_res]: # look at top 200 constraints
		i = restraint[0]-1
		j = restraint[1]-1
		if mean_array[i,j] <= args.upper_bound and std_array[i,j] > args.standard_deviation:	# seems to be a bit unsure
			filtered_restraints.append(restraint)


	# print restraint files
	print_files = []

	if len(lowstd_sr) > 0:
		print_files.append((lowstd_sr, "_converged_distances.cst", args.padding))
	if len(filtered_restraints) > 0:
		print_files.append((filtered_restraints, "_filtered_contactmaps.cst", -1))

	for ele in print_files:
		o = open( args.outfile + ele[1], "w")
		if ele[2] == -1:
			o.write(utils.print_restraints(ele[0], args.lower_bound, args.upper_bound, args.atom, args.function, args.parameter,seq))
		else:
			o.write(utils.print_padded_restraints(ele[0], args.padding, args.atom, args.function, args.parameter, seq))
		o.close()






if __name__ == main():
		main()
