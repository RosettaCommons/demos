#! /usr/bin/env python

import argparse
import Utils as utils

def main():

	parser = argparse.ArgumentParser(description='Extract top scoring restraints from contactmap', formatter_class=argparse.ArgumentDefaultsHelpFormatter )
	parser.add_argument( "cmp", help = "contactmap" )

	parser.add_argument("-ms", "--min_sep", help = "minimum sequence seperation", default = 5, type = int)
	parser.add_argument("-w", "--weighted", help = "use weighted atom pair restraints", action="store_true")

	opt = parser.add_argument_group( "Restraint Output Options" )
	opt.add_argument( '-o', help = "restraint file name" )
	opt.add_argument( '-r_f', '--r_func', help = "Restraint Function (default: %(default)s)", choices=["BOUNDED", "SIGMOID"], default="SIGMOID")
	opt.add_argument( '-r_ub', help = "upper bound (default: %(default)s)", type=float, default = 8 )
	opt.add_argument( '-r_lb', help = "lower bound (default: %(default)s)", type=float, default = 1.5 )
	opt.add_argument( '-r_atom', help = "restraint atoms", default = "CB", choices=["CA", "CB"] )
	opt.add_argument( '-r_num', help = "number of restraints", type = int )
	opt.add_argument( '-r_w', help="restraint function weight", default = 1, type = float )
	opt.add_argument( '-r_fasta', help = "fasta file (necessary if CB is selected)" )
	opt.add_argument( '-r_num_perc', help = "percentage of sequence length converted into restraints", type = float )

	args = parser.parse_args()

	cmp = utils.parse_cmpfile( args.cmp )
	restraints = utils.extract_cm_restraints( cmp , args.min_sep )
	fasta = utils.parse_fasta( args.r_fasta ) if args.r_fasta else ""

	if args.r_num or args.r_num_perc:
		num = (int(( len(fasta) * args.r_num_perc)* 1.0)/10) * 10 if args.r_fasta and args.r_num_perc else args.r_num
		res = restraints[:num]
	else:
		num = "all"
		res = restraints
	protein = args.cmp.split('/')[-1][:-4]

	name = protein  + "_ub" + str(args.r_ub) + "_lb" + str(args.r_lb) + "_w" + str(args.r_w) + "_m" +str(args.min_sep)+ "_"+ args.r_atom + "_" + str(num) + "_" + args.r_func + ".cst" if not args.o else args.o

	o = open( name, "w" )
	o.write( utils.print_restraints(res, args.r_lb, args.r_ub, args.r_atom, args.r_func, args.r_w, fasta, args.weighted) )
	o.close()

if __name__ == main():
		main()

