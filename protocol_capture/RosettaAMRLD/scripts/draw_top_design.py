#!/usr/bin/env python
##
##
## @author Yidan Tang (yidan.tang@vanderbilt.edu)

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys

def main(args):
	col_names = ["pdb", "Total_score", "Interface_delta_X", "LE", "SMILES"]
	designOUT = pd.read_csv(args.input, sep='\t', lineterminator='\n').dropna()
	designOUT.columns = col_names
	query = designOUT['SMILES'].tolist()
	mols = [Chem.MolFromSmiles(mol) for mol in query]
	scores = ["{:.2f}".format(x) for x in designOUT["LE"].tolist()]
	if len(query) < args.topN:
		img = Draw.MolsToGridImage(mols[:len(query)], legends=scores[:len(query)], molsPerRow=5)
	else:
		img = Draw.MolsToGridImage(mols[:args.topN], legends=scores[:args.topN], molsPerRow=5)
	img.save(args.output)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-i', '--input', help="The design file to extract SMILES.")
	parser.add_argument('--output', default='Top10.png', help="Output filename.")
	parser.add_argument('--topN', type=int, default=10, help="plot top N SMILES; default=10")

	args = parser.parse_args()
	main(args)
