#!/usr/bin/python
## make mammoth structure alignments


import string
from glob import glob
from sys import argv,stderr,exit,stdout,stdin
from os import popen,system
from os.path import exists
from operator import add
from math import sqrt

def getMedian(numericValues):
  theValues = sorted(numericValues)
  if len(theValues) % 2 == 1:
    Q2= theValues[(len(theValues)+1)/2-1]
  else:
    lower = theValues[len(theValues)/2-1]
    upper = theValues[len(theValues)/2]
    Q2= (float(lower + upper)) / 2
  L1=int((len(theValues)*25.0/100))
  L2=int((len(theValues)*75.0/100))
  Q1=(theValues[L1-1]+theValues[L1])/2
  Q3=(theValues[L2-1]+theValues[L2])/2
  return (Q2,Q1,Q3)


def mean(nums):
  if len(nums):
    return float( sum(nums) ) / float( len(nums) )
  else:
    return 0.0

def std(nums):
  if len(nums):
    N=len(nums)
    m=mean(nums)
    sqr_sum=0
    for num in nums:
      sqr_sum=sqr_sum+(num-m)*(num-m);
    return sqrt( float( sqr_sum )/(N-1) )

#main program

lines = stdin.readlines()
words = string.split(lines[0])

#figure out the columns with numbers
ct=0;
cols=[]
for word in words:
    try:
       float(word)
       cols.append(ct)
       ct=ct+1
    except ValueError:
       ct=ct+1

nums=[]
for col in cols:
    nums.append([])

for line in lines:
    words=string.split(line)
    ct=0
    for col in cols:
    	num=float(words[col])
   	nums[ct].append(num)
	ct=ct+1

med=[]
q1=[]
q3=[]
hi=[]
lo=[]

ave=[]

ct=0;
for col in cols:
  (Q2,Q1,Q3)=getMedian(nums[ct])
  med.append(Q2)
  q1.append(Q1)
  q3.append(Q3)
  hi.append(max(nums[ct]))
  lo.append(min(nums[ct]))
  average = mean(nums[ct])
  ave.append(average)
  ct=ct+1

res=[ med, q1, q3, hi, lo, ave ]
text=[ "median", "Q1", "Q3", "hi", "lo","mean" ];
ct=0;
keepcols=0
if keepcols:
  for r in res:
    print "%10s"%text[ct],
    ct=ct+1
    ct_col=0
    for col in cols:
      print "%8.3f"%r[ct_col],
      ct_col=ct_col+1
    print
else:
  for r in text:
    print "%8s"%r,
  print
  ct_col=0
  for col in cols:
    for r in res:
      print "%8.3f"%r[ct_col],
    print
    ct_col=ct_col+1
