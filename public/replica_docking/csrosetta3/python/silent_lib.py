##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'
###
###
### This file is part of the CS-Rosetta Toolbox and is made available under
### GNU General Public License
### Copyright (C) 2011-2012 Oliver Lange
### email: oliver.lange@tum.de
### web: www.csrosetta.org
###
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <http://www.gnu.org/licenses/>.
###
###

import string
from library import MissingInput

class ReadSilentData:
    def __init__( self, scores ):
        self.donetitles=False
        self.names = scores
		  if len(scores)<1:
			  raise MissingInput("at least one score-column name has to be supplied")

        for (i,n) in enumerate(scores):
            self.names[i] = (string.lower(n), -1)


    def read_line( self, line ):
        l=line
        token = string.split(l)
        if len(token) <= 0: return
  ## SCORE LINE HEADERS
        if (token[0] == "SCORE ") or ((token[0] == "SCORE:") and (token[1] == "score")):
    ## look for names

            #reset columns
            for (i,(n,j)) in enumerate(self.names):
                self.names[ i ] = (n, -1 )

            self.format_str=""
            #find columns
            for ti,t in enumerate(token):
                for (i,(n,j)) in enumerate(self.names):
                    if string.lower(t) == n:  ## if names match
                        self.names[i] = (n,ti)  ## remember index
                        self.format_str=self.format_str+"%10s "
            #check that all columns have been found
            Error = False
#            for i,(n,j) in enumerate(self.names):
        #        if j < 0:
        #            raise MissingInput("Error: Cannot find column named '%s' \n"%n)
            self.donetitles = True
            return None


  ## SCORE LINE
        if token[0] == "SCORE:":
			  if not self.donetitles:
				  raise MissingInput("Error: Cannot find column position, missing header line")
			  token = string.split(l)
			  str=""
			  if len(self.names)>1:
				  res=[]
				  for (n,i) in self.names:
						if i<0: res.append('nan')
						else: res.append(token[i])
				  return res
					#return [token[i] for (n,i) in self.names]
			  else:
				  i=self.names[0][1]
				  if i<0: return 'nan'
				  return token[i]
#                str=str+"%10s "%token[i]
#            return str
		  return None


def read_lowscore_tags( filename, col, num, sign=-1 ):
	scores=[]
	sfd=ReadSilentData([col,'description'])
	file=open( filename, 'r' )
	min_score = 1000000
	for l in file:
		if l[0:6]=="SCORE:":
			tags=sfd.read_line( l )
			if tags:
				scores.append((float(tags[0]),tags[1]))

#	print scores
	scores=sorted( scores, key=lambda x: x[0] )
	if num>len(scores):
		return [ i[1] for i in scores ]
	if sign<0:
		return [ i[1] for i in scores[0:num] ]
	else:
		return [ i[1] for i in scores[-num:] ]


def extract_decoys( in_file, out_file, tags ):
	file=open(in_file, 'r')
	out=open(out_file,'w')
	remarks=''
	score_header=''
	sequence_header=''
	for l in file:
		if l[0:len('REMAKR')]=="REMARK":
			remarks=remarks+l
			continue
		if l[0:len('SEQUENCE')]=='SEQUENCE':
			sequence_header=l
			continue
		if l[0:6]=="SCORE:" and 'score' in l:
			score_header=l
			continue

		tt=l.split()
		if tt[-1] in tags:
			if sequence_header:
				out.write(sequence_header)
				sequence_header=''
			if score_header:
				out.write(score_header)
				score_header=''
			if remarks:
				out.write(remarks)
				remarks=''
			out.write(l)
		else: remarks=''


