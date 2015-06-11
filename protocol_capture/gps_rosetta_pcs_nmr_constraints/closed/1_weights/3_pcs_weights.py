#!/usr/bin/env python

import sys, os, argparse

# run as <python 3_pcs_weights.py -h> or 
# python 3_pcs_weights.py  --csfile dvp_cs.csc --pcsfile dvp_pcs.csc --ntags 4
# Generates two output files pcs.wts and t13_pcs.wts
# 'pcs.wts' is used to patch the centroid scoring function in the low-resolution phase of fragment assembly step.
# 't13_pcs.wts' is used to rescore all-atom silent files with Talaris2013 wts along with pcsTs[1-4] wts for each of the lanthanide tag sites.

ten_percent=3
# ten_percent=100
# replace ten_percent with 100 if 1000 structures are generated in 1_Rosetta_wts.sh else use 10 for 100 structures.

################### Parse Cmdline args ###################

parser = argparse.ArgumentParser(description='calculate PCS weights relative to Rosetta centroid scoring function, run as <python 3_pcs_weights.py  --csfile dvp_cs.csc --pcsfile dvp_pcs.csc --ntags 4> .')
parser.add_argument('--csfile', type=argparse.FileType('r'), help='input CS score file')
parser.add_argument('--pcsfile', type=argparse.FileType('r'), help='input PCS unity score file')
parser.add_argument('--ntags', type=int, help='number of tags used')

args = parser.parse_args()

################### Parse Score files ###################
cs_score= []
cslines=args.csfile.readlines()
for csline in cslines:
    cs_data=csline.split()
    if cs_data[0]=='SCORE:':
        if cs_data[1]!='total_score':
            cs_score.append(float(cs_data[1]))
args.csfile.close()

pcs1w=[]
for i in range(0, args.ntags):
    pcs1w.append([])
pcslines=args.pcsfile.readlines()
for pcsline in pcslines:
    pcs_data=pcsline.split()
    if pcs_data[0]=='SCORE:':
        if pcs_data[1]!='score':
            for i in range (0, args.ntags):
                pcs1w[i].append(float(pcs_data[i+2]))    
args.pcsfile.close()

################### Calculate Weights ###################

output=open("pcs.wts",'w',1)
output1=open("t13_pcs.wts",'w',1)

cs_low, cs_high = 0,0
cs_score.sort()

for i in range (1,ten_percent+1):
    cs_low = cs_low + cs_score[i-1]
    cs_high =cs_high + cs_score[-1*(i)]
a_high = cs_high/(100.00)
a_low = cs_low/(100.00)

pcs_weights = []
for i in range(0,args.ntags):
    pcs_ts=(pcs1w[i])
    pcs_ts.sort()      
    pcs_low, pcs_high = 0, 0
    for j in range (1,ten_percent+1):
        pcs_low = pcs_low + pcs_ts[j-1]
        pcs_high = pcs_high + pcs_ts[-1*(j)]

    c_high = pcs_high/(100.00)
    c_low = pcs_low/(100.00)    
    pcs_weights.append((a_high-a_low)/(c_high-c_low))

################### Write output ###################

#append PCS wights to Talaris2013 all atom weights
output1.write("METHOD_WEIGHTS ref 0.592942 0.354993 -1.28682 -1.55374 0.43057 0.140526 0.357498 0.831803 -0.287374 0.602328 0.158677 -0.94198 -0.219285 -1.17797 -0.14916 0.176583 0.16454 0.744844 0.92933 0.131696 \nfa_atr 0.8 \nfa_rep 0.44 \nfa_sol 0.75 \nfa_intra_rep 0.004 \nfa_elec 0.7  \npro_close 1  \nhbond_sr_bb 1.17 \nhbond_lr_bb 1.17  \nhbond_bb_sc 1.17 \nhbond_sc 1.1  \ndslf_fa13 1.0  \nrama 0.2  \nomega 0.5  \nfa_dun 0.56  \np_aa_pp 0.32  \nref 1  \n")
for i in range(0,len(pcs_weights)):
    tmp2=round((pcs_weights[i]/args.ntags),3)    
    print "Scaled weight PCSTs"+str(i+1)+" "+str(tmp2)
    tmp= "pcsTs"+str(i+1)+" = "+str(tmp2)+"\n"
    output.write(tmp)
    tmp= "pcsTs"+str(i+1)+"  "+str(tmp2)+"\n"
    output1.write(tmp)
output.close()
output1.close()

