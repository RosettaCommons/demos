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
#!/bin/env tcsh -f
#
# CS-ROSETTA: System for Chemical Shifts based protein structure prediction using ROSETTA
#

#path names to find the CS-Rosetta toolbox files
setenv csrosettaDir /work/olange/csrosetta3
setenv csrosettaCom ${csrosettaDir}/com_alias

#provide CS-Rosetta applications in system-path
setenv PATH ${PATH}:${csrosettaCom}:/work/olange/rosetta/rosetta_source/bin

#make sure that python modules can be found
 if ( ! ($?PYTHONPATH) ) then
    setenv PYTHONPATH
 endif
setenv PYTHONPATH ${PYTHONPATH}:${csrosettaDir}/python

#this is useful for benchmarking and used by the scripts setup_target.py display_target.py and setup_run.py
setenv CS3_BENCH_TARGETLIB $HOME/cs_targetlib

#the path to the ROSETTA installation (yields defaults for setup_run options -database and -binaries_prefix)

	setenv ROSETTA3_PATH /work/olange/rosetta
	setenv ROSETTA3_BIN /work/olange/rosetta/rosetta_source/bin
setenv ROSETTA3_DB /work/olange/rosetta/rosetta_database

#Making sure that TALOS+ commands are available, too
setenv TALOSP_DIR /work/olange/csrosetta3/NMRPipe//talosplus

#Sparta+ installation
setenv SPARTAP_DIR /work/olange/csrosetta3/NMRPipe//spartaplus

#Sparta+ installation
setenv PROMEGA_DIR /work/olange/csrosetta3/NMRPipe//NMRPipe/promega
