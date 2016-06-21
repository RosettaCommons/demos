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

import string
from os.path import basename,dirname

class Job:
	def __init__(self, fn_job ):
		self.job=splitext(basename(fn_job))[0]
		self.jobpath=dirname(fn_job)
		jobbase=self.jobpath+"/"+self.job
            #platform_file = $jobpath/$( echo $jobname | awk -v FS="_" '{print $1}' ).generic
		tags=string.split( self.job, "_" )
		self.platform_file= self.jobpath+"/"+tags[0]+".generic"

	def __str__(self ):
		return "job %s  -- for platform %s\n"%(self.job, self.platform_file )
