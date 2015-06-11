#! /usr/bin/env python


import argparse
import Utils as utils
import numpy as np


def parse_csv_store_matrix( fi, output, fasta, value_int, a1_int, a2_int, offset, separator=" " ):
	x_l,y_l, dict = [], [], {}
	f = open (fi, "r")
	for line in f:
		split = line.split( separator )
		try:
			if len(split) > 1:
				x = int(split[a1_int])-offset+1
				y = int(split[a2_int])-offset+1
				x_l.append(x)
				y_l.append(y)
				dict[(x,y)] = round(float(split[value_int]),3)
		except:
			print "Problems with line: {0}".format( line )
	fasta_len = len(utils.parse_fasta( fasta ))
	ma = np.zeros([fasta_len, fasta_len])
	max_val = max(dict.values())
	min_val = min(dict.values())
	for ele in dict:
		ma[ele[0]-1, ele[1]-1] = dict[ele]
		ma[ele[1]-1, ele[0]-1] = dict[ele]

	o=open (output, "w")
	for line in ma:
		o.write(" ".join(map(str,line)))
		o.write("\n")
	o.close()


def main():

	parser = argparse.ArgumentParser(description='Parses EVFold all-by-all residue pairings ({jobname}_{scoringmethod}.txt) and generate lxl contact map',formatter_class=argparse.ArgumentDefaultsHelpFormatter )
	parser.add_argument( '-i', "--input", help="Input residue-residue pairing scores (*_PLM.txt)", required=True )
	parser.add_argument( "-f", "--fasta", help="Fasta file for length comparison", required="True")
	parser.add_argument( "-o", "--output", help="Output Contactmap", default = "contactmap.cmp")
	parser.add_argument( "-c", "--column", help = "index of score column", default = 5, type=int)
	parser.add_argument( "-aa1", help = "index of 1st residue column", default = 0, type=int)
	parser.add_argument( "-aa2", help = "index of 2nd residue column", default = 2, type=int)
	parser.add_argument( "-offset", help = "residue offset", default = 1, type=int)
	parser.add_argument( "-d", "--delimiter", help="pairing score file delimiter", default = " ")

	args = parser.parse_args()

	parse_csv_store_matrix(args.input, args.output, args.fasta, args.column, args.aa1, args.aa2, args.offset, args.delimiter)


if __name__ == main():
	main()
