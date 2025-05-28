#!/usr/bin/env python3
##
## @author Yidan Tang (yidan.tang@vanderbilt.edu)

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os, sys

def load(filename, c):
	trial = []
	score = []
	with open(filename, 'r') as f:
		for line in f:
			if line.startswith('Stats'):
				break
			splitted = line.split(',')
			trial.append(float(splitted[0]))
			score.append(float(splitted[c]))
	return trial, score

def main(args):
	data = []
	best_pose = []
	for f in args.inputs:
		trial, score = load(f, args.column)
		min_sc = min(score)
		min_i = score.index(min_sc)
		best_pose.append(trial[min_i])
		data.append(trial)
		data.append(score)
	a = np.array([np.array(x) for x in data])

	if args.separate_image == False:
		for i in range(len(args.inputs)):
			plt.plot(a[2*i], a[2*i+1], label=os.path.splitext(args.inputs[i])[0])
		plt.xlabel('Monte Carlo trial number')
		plt.ylabel(args.y)
		plt.legend()
		plt.savefig('interface_MCtrial.png')
	else:
		for i in range(len(args.inputs)):
			fig = plt.figure()
			ax = fig.add_subplot()
			ax.plot(a[2*i], a[2*i+1])
			axes = plt.gca()
			axes.set_xlim([0, args.mc])
			plt.xlabel('Monte Carlo trial number')
			plt.ylabel(args.y)
			plt.savefig(os.path.splitext(args.inputs[i])[0] + '.png')
			plt.close(fig)
	if args.best_trial == True:
		box = plt.figure()
		plt.violinplot(best_pose)
		plt.gca().set_ylim([0, args.mc])
		plt.savefig('best_trial.png')
		plt.close(fig)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('inputs', nargs='*', help="The scores to plot.")
	parser.add_argument('--separate_image', action='store_true', help="Plot scores from files in separate plots.")
	parser.add_argument('--column', type=int, default=2, help="Column to plot; default=1")
	parser.add_argument('--mc', type=int, default=2000, help="max MC cycle; default=1000")
	parser.add_argument('-y', default="LE2(REU)", help="y axis label")
	parser.add_argument('--best_trial', action='store_true', help="Plot voilin plots to show best trial no. distribution.")

	args = parser.parse_args()
	main(args)
