#!/usr/bin/python

import sys
import os
import string

if len(sys.argv) < 2:
  print "Syntax: %s <silentfile> <columnname1> [  <columnname2> ... ]"%sys.argv[0]
  sys.exit(0)

file = open(sys.argv[1],"r")

#lines = string.split(file.read(),"\n")

showfilename = False;

names = []
for arg in sys.argv[2:]:
  if arg == "-f":
    showfilename = True
    continue
  names.append((arg,-1))

donetitles = False

for (i,(n,j)) in enumerate(names):
  names[i] = (string.lower(n), j)

if len(names)<=0:
  sys.stderr.write("Warning: not columns names supplied: defaulting to rms vs score\n")
  names.append( ("rms", -1) )
  names.append( ("score", -1) )


for l in file:
  token = string.split(l)
  if len(token) <= 0: continue
  ## SCORE LINE HEADERS
  #if not donetitles:
  if (token[0] == "SCORE ") or ((token[0] == "SCORE:") and (token[1] == "score")):
    ## look for names
    for ti,t in enumerate(token):
      for (i,(n,j)) in enumerate(names):
        if string.lower(t) == n:  ## if names match
          names[i] = (n,ti)  ## remember index

    Error = False
    for i,(n,j) in enumerate(names):
      if j < 0:
        sys.stderr.write("Error: Cannot find column named '%s' \n"%n)
        Error = True

    if Error: sys.exit(1)
    donetitles = True
    continue


  ## SCORE LINE
  if token[0] == "SCORE:":
    token = string.split(l)
    if showfilename: sys.stdout.write("%s "%sys.argv[1])
    for (n,i) in names:
      sys.stdout.write("%10s "%token[i])
    sys.stdout.write("\n")
