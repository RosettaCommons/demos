#!/blue/meilerlab/apps/Linux2/x86_64/bin/python2.5

from __future__ import division
from optparse import OptionParser
from rosettautil.rosetta import rosettaScore
from rosettautil.util import fileutil
import sys
import math

usage = "%prog [options] --silent <input silent file> --membrane -n <num residues in protein> --outfile <outfile>"
parser = OptionParser(usage)
parser.add_option("--membrane",dest="membrane",help="if want rmsSSE100 or rmsTM100 calculated, specify --membrane",default=False,action="store_true")
parser.add_option("--silent",dest="silent",help="input a single silent file or may work for score file",default="")
parser.add_option("--silent_list",dest="silent_list",help="input a list of silent files or list of score files",default="")
parser.add_option("--pdbs",dest="pdb_list",help="input list of PDB files",default="")
parser.add_option("--outfile",dest="outfile",help="output file",default="")
parser.add_option("-n",dest="n",help="number of residues in protein")
parser.add_option("--rms_tag",dest="rms_tag",help="tag (TM or SSE",default="TM")
(options,args) = parser.parse_args()

if options.silent == "" and options.pdb_list == "" and options.silent_list=="":
	parser.error("you must specify --silent, --silent_file_list, or --pdbs")
if options.n == 0:
	parser.error("you must put the number of residues, which cannot equal 0")

#data = [] #list of tuples in form (tag,x_score,y_score)
rms100 = []
rmsSSE100 = []
rms_list =[]
rmsSSE_list =[]

if options.silent != "":
	scores = rosettaScore.SilentScoreTable()
	scores.add_file(options.silent)
	print "added score file"
	print scores
	rms = scores.score_generator("rms_"+str(options.rms_tag))
	print rms
	score = scores.score_generator("total_score")
	print rms
	# make another generator because can't reuse rms generator once it's at the end (for later)
	rms_list = [x for x in scores.score_generator("rms_"+str(options.rms_tag))]
	score_list = [s for s in scores.score_generator("score")]
	for rms_i in rms:
		print rms_i
		if rms_i[1] != "rms_"+str(options.rms_tag):
			rms100_i = float(rms_i[1])/(1+math.log(math.sqrt(float(options.n)/100)))
		rms100.append(rms100_i)
#	for score_i in score:
#		if score_i[1] != "score":
		
	if options.membrane:
		rmsSSE = scores.score_generator("rms_"+str(options.rms_tag))
		# same goes here
		rmsSSE_list = [x for x in scores.score_generator("rms_"+str(options.rms_tag))]
		score_list = [s for s in scores.score_generator("score")]
		for rmsSSE_i in rmsSSE:
			if rmsSSE_i[1] != "rms_"+str(options.rms_tag):
				rmsSSE100_i = float(rmsSSE_i[1])/(1+math.log(math.sqrt(float(options.n)/100)))
			rmsSSE100.append(rmsSSE100_i)

elif options.silent_list != "":
	silent_list = fileutil.universal_open(options.silent_list,"rU")
	for path in silent_list:
		scores = rosettaScore.SilentScoreTable()
		scores.add_file(path.rstrip())
		rms = scores.score_generator("rms_"+str(options.rms_tag))
		score = scores.score_generator("score")
		for rms_i in rms:
			rms100_i = float(rms_i[1])/(1+math.log(math.sqrt(float(options.n)/100)))
			rms100.append(rms100_i)
		if options.membrane:
			rmsSSE = scores.score_generator("rms_"+str(options.rms_tag))
			for rmsSSE_i in rmsSSE:
				rmsSSE100_i = float(rmsSSE_i[1])/(1+math.log(math.sqrt(float(options.n)/100)))
				rmsSSE100.append(rmsSSE100_i)
	silent_list.close()

if options.pdb_list != "":
	pdb_list = fileutil.universal_open(options.pdb_list,"rU")
	for pdb in pdb_list:
		scores = rosettaScore.ScoreTable(pdb.rstrip())
		rms = scores.score_generator("rms_"+str(options.rms_tag))
		score = scores.score_generator("score")
		for rms_i in rms:
			rms100_i = float(rms_i[1])/(1+math.log(math.sqrt(float(options.n)/100)))
			rms100.append(rms100_i)	
		if options.membrane:
			rmsSSE = scores.score_generator("rms_"+str(options.rms_tag))
			for rmsSSE_i in rmsSSE:
				rmsSSE100_i = float(rmsSSE_i[1])/(1+math.log(math.sqrt(float(options.n)/100)))
				rmsSSE100.append(rmsSSE100_i)
	pdb_list.close()

outfile = open(options.outfile,'w')
if options.membrane:
	outfile.write("tag\tscore\trms\trms100\trms"+str(options.rms_tag)+"\trms100_"+str(options.rms_tag)+"\n")
	for i in range(len(rms100)):
		outfile.write(str(rms_list[i][0])+"\t"+str(score_list[i][1])+"\t"+str(rms_list[i][1])+"\t"+str(rms100[i])+"\t"+str(rmsSSE_list[i][1])+"\t"+str(rmsSSE100[i])+"\n")
else:
	outfile.write("tag\tscore\trms\trms100\n")
	for i in range(len(rms100)):
		outfile.write(str(rms_list[i][0])+"\t"+str(score_list[i][1])+"\t"+str(rms_list[i])+"\t"+str(rms100[i])+"\n")
outfile.close()
