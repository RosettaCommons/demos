#!/Library/Frameworks/EPD64.framework/Versions/7.3/bin/python

import glob

for file in glob.glob('distance_output/*'):
	infilename = file.split("/")[1][0:5]
	outfile = open('pdb_lists_filtered/'+infilename+'.txt','w')
	
	struct_dict = {}
	f_data = open(file,'r')
	for line in f_data.readlines():
		struct_dict[line.split()[0]]=line.split()[1:]
	
	#Don't use the alignment biased list for filtering in benchmark studies 
	#cutoff_list = [1.155,1.724,1.91,2.927,1.911,1.799,1.121,1.029,1.093,0.773,0.818,0.729,1.624,1.615,1.615,2.236,1.337,0.857,0.511,0.815,1.507,1.168,0.967,1.268,1.243,0.876,0.806,0.876,1.142]
	cutoff_list = [1.155,1.724,1910,2927,1.911,1.799,1.121,1.029,1.093,0.773,0.818,0.729,1624,1615,1.615,2.236,1.337,0.857,0.511,0.815,1.507,1.168,0.967,1.268,1.243,0.876,0.806,0.876,1.142]
	
	
	for key in struct_dict.keys():
		n=0
		m=0
		for distance in struct_dict[key]:
			if float(distance) > float(cutoff_list[n] + 0.5): m+=1
			n += 1
		if m==0: outfile.write(key+'\n')
	outfile.close()
		
