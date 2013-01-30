#!/bin/bash
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

mkdir -p pick_fragments
cd pick_fragments
#1. prepare input files
wget http://rest.bmrb.wisc.edu/bmrb/NMR-STAR2/15339 -O 2jrm.bmrb
bmrb2talos 2jrm.bmrb > 2jrm.tab

#3.1 figure out flexible regions in original protein sequence
echo running talos+ to figure out flexible regions
mkdir -p untrimmed
cd untrimmed
ln -sf ../2jrm.tab
talos+ -in 2jrm.tab > talos.log
cat pred.tab
cd ..

#3.2 from the pred.tab file we found that we want to trim off residue 1-6 at N-terminal and 53-60 at C-terminal
renumber_talos -s 7 -e 52 2jrm.tab 2jrm_trim.tab
talos2fasta 2jrm_trim.tab > 2jrm_trim.fasta

#3.3 now pick the fragments
echo picking fragments...
mkdir -p fragments
cd fragments
ln -sf ../2jrm_trim* .
pick_fragments -cs 2jrm_trim.tab -hom
cd ..


