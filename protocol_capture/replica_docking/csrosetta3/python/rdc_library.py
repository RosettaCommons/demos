#!/usr/bin/env python2.7
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


class RDC_Line:
    def __init__( self, res1, atom1, res2, atom2, rdc ):
        self.atom1=atom1
        self.atom2=atom2
        self.res1=res1
        self.res2=res2
        self.rdc=rdc

class RDC_Data:
    def __init__( self ):
		 self.data=[]
		 pass

    def read_file( self, file ):
        self.data=[]
        lines = open( file  ).readlines()
        for line in lines:
            if line[0]=='#': continue
            if len(line) < 5: continue
            cols = line.split()
            self.data.append( RDC_Line( int( cols[0] ), cols[1], int( cols[2] ), cols[3], float( cols[4] ) ) );

    #return list of the RDC values
	 def size(self):
		 return len(self.data)

    def rdcs( self ):
        r =[]
        for l in self.data:
            r.append( l.rdc )
        return r

    def estimate_Da_and_R_hist( self, rdcs=None, binwidth=None ):
        if not rdcs:
            rdcs=self.rdcs()

        mindata = min(rdcs)
        maxdata = max(rdcs)

        if not binwidth:
            binwidth = (maxdata - mindata) / 50.0;

        histogram = []
        bincenter = []
        currentbincenter = mindata + binwidth/2.0;
        numbins =  int((maxdata-mindata)/binwidth)
        for bin in range(numbins):
            histogram.append( 0.0 )
            bincenter.append( currentbincenter )
            currentbincenter += binwidth

        total = 0
        for rdc in rdcs:
            bin = int( (rdc - mindata)/binwidth )
            if bin<0       : bin = 0
            if bin>=numbins: bin = numbins - 1
            histogram[bin] += 1
            total = total + 1


        for bin in range(numbins):
            histogram[bin] = histogram[bin]/total * 100.0

#compute Da and R
        if max(abs(bincenter[0]),abs(bincenter[numbins-1]))==abs(bincenter[0]):
            Dzz=abs(bincenter[0])*abs(bincenter[0])/bincenter[0]
            Dyy=abs(bincenter[numbins-1])*abs(bincenter[numbins-1])/bincenter[numbins-1]
        else:
            Dzz=abs(bincenter[numbins-1])*abs(bincenter[numbins-1])/bincenter[numbins-1]
            Dyy=abs(bincenter[0])*abs(bincenter[0])/bincenter[0]
        Dxx=-Dzz-Dyy
        R=(1+2*Dxx/Dzz)*2/3
        rdc_range=maxdata-mindata
        inv_rdc_range2=1/(rdc_range*rdc_range)
        return Dxx,R,rdc_range,inv_rdc_range2

