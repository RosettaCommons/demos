#!/usr/bin/env python2.7

# Extract Data and plot figures from *.out file
import argparse

parser = argparse.ArgumentParser(description="Extract Data and plot figures from *.out file",
                                 add_help=True)

parser.add_argument("-file", help="out file");
parser.add_argument("-datax", help="used data x", default='rms');
parser.add_argument("-datay", help="used data y", default='angle_constraint');
parser.add_argument("-title", help="figure title", default='figure x');
args = parser.parse_args()

fileout=open(args.file)
fileList = fileout.readlines()

lengthfile=len(fileList)
used_datax=[]
used_datay=[]
data_numx=0
data_numy=0
SCORE_num=0


#for line in fileList:
#	tags=line.split;
#	if tags[0]!='SCORE': continue
#	if tags[1]=='score':
		#update_column_number
#	else:
#		score_data.append(float(tags[score_col]))

#extract data from files
for i in range (0,lengthfile):
	if fileList[i].split()[0]=='SCORE:' and fileList[i].split()[1]=='score':
		SCORE_num=i
		break
print SCORE_num

lengthline=len(fileList[SCORE_num].split())
print lengthline
for i in range (0, lengthline):
	if fileList[SCORE_num].split()[i]==args.datax:
		data_numx=i
		break
print data_numx
for i in range (0, lengthline):
	if fileList[SCORE_num].split()[i]==args.datay:
		data_numy=i
		break
print data_numy
for i in range (0,lengthfile):
	if fileList[i].split()[0]=='SCORE:' and fileList[i].split()[1]!='score':
		used_datax.append(float(fileList[i].split()[data_numx]))
		used_datay.append(float(fileList[i].split()[data_numy]))

#plot the figure
import matplotlib.pyplot as plt
plt.plot(used_datax, used_datay, 'ro')
plt.xlabel(args.datax)
plt.ylabel(args.datay)
plt.title(args.title)
plt.show()
