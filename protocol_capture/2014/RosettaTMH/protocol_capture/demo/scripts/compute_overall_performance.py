#!/bin/env python2.7
from optparse import OptionParser
import numpy
import math
import re
from operator import itemgetter

# options stuff
usage = "%prog [options]"
parser = OptionParser(usage)
parser.add_option("--filelist",dest="filelist",help="filelist",default="list")
parser.add_option("--metric",dest="metric",help="header of column wanting to average",default="rms_dssp_bcl100")
parser.add_option("--score_fraction",dest="score_fraction",help="what fraction of models do you consider for TP etc",default="0.10")
parser.add_option("--quality_fraction",dest="quality_fraction",help="what fraction of models do you consider for TP etc",default="0.10")
parser.add_option("--p_ratio",dest="ratio",help="(p+n)/p. this sets your max enrichment",default="10.0")
parser.add_option("--enrichment_output",dest="enrichment_output",help="file to output enrichment values for each pdb and weight",default="enrich.txt")
parser.add_option("--outfile",dest="outfile",help="outfile",default="out.txt")
(opts,args) = parser.parse_args()

filelist = open(opts.filelist,'r')

# parse file list
def parse_filelist(filelist):
	metric = str(opts.metric)
	# loop through filenames in filelist
	pdb = ""
	score = 0.0
	rmsd = 0.0
	info_tuple = ()
	info_list = [] # list of tuples, where tuple[0] = score, tuple[1] = rmsd, tuple[2] = tag
	tag_and_values = {}
	for filename in filelist:
		weight_dict = {}
		filename = filename.rstrip('\n')
		split_filename = re.split('_',filename)
		pdb = split_filename[0]
		weight = split_filename[4]
		# open the file of interest and get column needed
		filename = filename.rstrip('\n')
		file = open(filename,'r')
		column_index = 0
		for index,line in enumerate(file):
			split_line = line.split()
			if index == 0:
				metric_field = split_line.index(opts.metric)
				score_field = split_line.index("score")
				tag_field = split_line.index("tag")
			else:
				rmsd = float(split_line[metric_field])
				score = float(split_line[score_field])
				tag = str(split_line[tag_field])
				info_tuple = (score,rmsd,tag)
				info_list.append(info_tuple)
		file.close()
		#print filename,pdb,method,len(rmsd_list)
		#for this method, you have this list of values
		weight_dict[weight] = info_list
		info_list = []
		try:
			tag_and_values[pdb].update(weight_dict)
		except KeyError:
			tag_and_values[pdb] = weight_dict
		weight_dict = {}
	# tag_and_values is a dictionary of dictionaries. The key of the outer dict is a pdb. The key of the inner dictioanry (value) is a weight, and the value of the inner dict is a list of rmsd values
	return tag_and_values

# take in dictionary of dictionaries and determine how many tp, fp, tn, fn for each pdb and for each weight
def process_data(tag_and_values):
	pdb_dict = {}
	new_weight_dict = {}
	sorted_score_list = []
	sorted_quality_list = []
	num_score_filter_models = 0
	num_quality_filter_models = 0
	top_scoring = []
	top_quality = []
	tp = []
	fp = []
	tn = []
	fn = []
	#loop through each pdb in tag_and_values 
	for pdb,weight_dict in tag_and_values.iteritems():
		# loop through weight_dict, where key is a weight, value is a list of tuples
		for weight,list in weight_dict.iteritems():
			# sort the score and rmsd lists
			sorted_score_list = sorted(list,key=itemgetter(0))
			sorted_quality_list = sorted(list,key=itemgetter(1))
			num_score_filter_models = int(round(float(opts.score_fraction) * len(sorted_score_list)))
			num_quality_filter_models = int(round(float(opts.quality_fraction) * len(sorted_quality_list)))
			# put top scoring models in top_scoring. do the same for top quality models
			for i in range(0,num_score_filter_models):
				top_scoring.append(sorted_score_list[i])
			for i in range(0,num_quality_filter_models):
				top_quality.append(sorted_quality_list[i])
			# set tp, fp, tn, fn
			for model in list:
				# low scoring and low rmsd
				if model in top_scoring and model in top_quality:
					tp.append(model)
				# low scoring and high rmsd
				elif model in top_scoring and model not in top_quality:
					fp.append(model)
				# high scoring and low rmsd
				elif model not in top_scoring and model in top_quality:
					fn.append(model)
				# high scoring and high rmsd
				elif model not in top_scoring and model not in top_quality:
					tn.append(model)
			# dictionary where key = weight and value = tuple of 4 lists
			new_weight_dict[weight] = (tp,fp,tn,fn)
			# reset variables
			tp = []
			tn = []
			fp = []
			fn = []
			sorted_score_list = []
			sorted_quality_list = []
			top_quality = []
			top_scoring = []
		# pdb_dict is a dictionary of dictionaries. key is pdb, value is dictionary. key of dictionary2 = weight, value is tuple of lists
		try:
			pdb_dict[pdb].update(new_weight_dict)
		except KeyError:
			pdb_dict[pdb] = new_weight_dict
		new_weight_dict = {}
	return pdb_dict

# compute enrichment for each pdb and weight
def compute_enrichment(subdivided_data):
	tp = 0
	fp = 0
	tn = 0
	fn = 0
	enrichment = 0.0
	ratio = float(opts.ratio)
	weight_dict = {}
	pdb_dict = {}
	# loop through subidivided_data, where key = pdb,value is dictionary
	for pdb,weight_dict in subdivided_data.iteritems():
		# loop through weight_dict, where key = weight, value = tuple of lists
		for weight,tuple in weight_dict.iteritems():
			tp = float(len(tuple[0]))
			fp = float(len(tuple[1]))
			tn = float(len(tuple[2]))
			fn = float(len(tuple[3]))
			# compute enrichment
			enrichment = float((tp/(tp+fp)) * ratio)
			# enrichment for each weight
			weight_dict[weight] = enrichment
		# each pdb has its own dictionary
		try:
			pdb_dict[pdb].update(weight_dict)
		except KeyError:
			pdb_dict[pdb] = weight_dict
		weight_dict = {}
	return pdb_dict

# get overall performance for each weight as one value
def compute_overall_performance(enrichments):
	sum = 0.0
	performance = 0.0
	dict = {}
	weight_list = []
	enrichment_list = []
	final_result = {}
	# loop through pdb_dict, where key = pdb, value = dict
	for pdb,weight_dict in enrichments.iteritems():
		# loop through weight_dict, where key = weight, value = enrichment
		for weight,enrichment in weight_dict.iteritems():
			# append weight to weight_list unless already found
			if weight in weight_list:
				continue
			else:
				weight_list.append(weight)
	# loop through weights in weight_list
	for weight in weight_list:
		# for every pdb, find the enrichment for this particular weight in the weight_list
		for pdb,weight_dict in enrichments.iteritems():
			for new_weight,enrichment in weight_dict.iteritems():
				if new_weight == weight:
					enrichment_list.append(enrichment)
				else:
					continue
		# now each weight has its list of enrichments (len(list) == # proteins)
		try:
			dict[weight].append(enrichment_list)
		except KeyError:
			dict[weight] = enrichment_list
		enrichment_list = []
	# the performance for a given weight is given by the sqrt of the sum of enrichments for that weight
	for weight,list in dict.iteritems():
		print weight,len(list)
		#sum = math.fsum(list)
		performance = numpy.mean(list)
		sd = numpy.std(list)
		final_result[weight] = (performance,sd)
	# final result is a dictionary, where key = weight, value = (mean,sd)
	return final_result 

# parse file list and get rmsds
all_rmsds = parse_filelist(filelist)
for key,value in all_rmsds.iteritems():
	for index,list in value.iteritems():
		pass
		#print key,index,len(list)
		for i in list:
			pass
			#print key,index,len(list),i

# take in dictionary of dictionaries and determine how many tp, fp, tn, fn for each pdb and for each weight
subdivided_data = process_data(all_rmsds)
for pdb,weight_dict in subdivided_data.iteritems():
	pass
	for weight,tuple in weight_dict.iteritems():
		#pass
		print pdb,weight,len(tuple),len(tuple[0]),len(tuple[1]),len(tuple[2]),len(tuple[3])

# compute enrichment for each pdb and weight
enrichments = compute_enrichment(subdivided_data)

enrich_file = open(opts.enrichment_output,'w')
enrich_file.write("pdb\tweight\tenrichment\n")
for pdb,weight_dict in enrichments.iteritems():
#	pass
	for weight,enrichment in weight_dict.iteritems():
		enrich_file.write(str(pdb)+"\t"+str(weight)+"\t"+str(enrichment)+"\n")
enrich_file.close()

overall_performance = compute_overall_performance(enrichments)

outfile = open(opts.outfile,'w')
outfile.write("weight\tavg_enrichment\tstdev\n")
for weight in overall_performance.keys():
	outfile.write(str(weight)+"\t"+str(overall_performance[weight][0])+"\t"+str(overall_performance[weight][1])+"\n")
outfile.close()
