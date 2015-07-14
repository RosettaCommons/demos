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
import library
import os
import math
import string
#################
def cst_map_rename_atom( atom, map_atom='CEN' ):
	if atom[0:2]=='QA':
		return atom, 1
	if atom[0:2]=='QQ':
		return map_atom,6
	if atom[0]=='Q':
		return map_atom,3
	if atom[0:3] in ['HE1','HE2','HE3','HZ3','HZ2','HD1','HD2','HD3','HG1','HG2','HG3']:
		return map_atom,1
	if atom[0:2] in ['HG','HE','HB','HH','HZ']:
		return map_atom,1
	return atom, 1

def cst_map_to_CB( file, outfile ):
	cst_map_to_CEN( file, outfile, map_atom='CB', mode='historical' )

def cst_map_to_CEN( file, outfile, map_atom='CEN', mode='simple', fasta=None ):
	print 'cst_map_to_CEN...'
	if not mode in ['simple','simple_short','historical']:
		print 'defer to modern aa-dep method'
		cst_map_to_CEN_by_aa_type( file, outfile, fasta, mode=mode )
		return

	lines = open(file).readlines()
	if isinstance(outfile, str):
      outfile=open(outfile,'w')

	for line in lines:
      if (line[0]=="#"): continue;
      l=string.split(line);
      try:
          resi1 = int( l[2] );
          resi2 = int( l[4] );
      except:
          print '\n\n expected integer in col 3 and 5 of:\n'
          print l[2],l[4],'\n',l
          exit()

      atom1,met1 = cst_map_rename_atom( l[1], map_atom=map_atom );
      atom2,met2 = cst_map_rename_atom( l[3], map_atom=map_atom );
      bound_low = float( l[6] );
      try:
          bound_up = float( l[7] );
          bound_sd = float( l[8] );
      except:
          print '\n\n expected float in col 8 and 9 of:\n'
          print l[7],'\n',l
          exit()

      pad = met1 * met2
      hyd = (met1 > 1) + (met2 > 1 )
      #print line, pad, hyd, pad**(1.0/6)
      bound_sd_ori = bound_sd
      bound_up_ori = bound_up

		if mode=='historical':
			bound_sd = bound_sd + 0.5*((pad-1)**(1.0/6))
			bound_up = (bound_up+0.5*hyd)*(pad**(1.0/6))
			atom1=atom1.replace("CEN","CB")
			atom2=atom2.replace("CEN","CB")
		elif mode=='simple':
			bound_sd = bound_sd+hyd
			bound_up = bound_up+hyd
		elif mode=='simple_short':
			bound_sd = bound_sd+0.15*hyd
			bound_up = bound_up+0.25*hyd

      cst_type='AtomPair'
      if ( len(atom1)>2 or len(atom2)>2 ) or atom1[0]=='Q' or atom2[0]=='Q':  #HA1 HA2 HA3
          cst_type='AmbiguousNMRDistance'
      if ( resi1==1 or resi2==1 ):     #possibly H on resi 1
          cst_type='AmbiguousNMRDistance'

      outfile.write('%20s %5s %5d %5s %5d BOUNDED %5.3f %5.3f %5.3f NOE; rawdata %5.3f %5.3f\n'%(cst_type,atom1,resi1,atom2,resi2,bound_low,bound_up,bound_sd,bound_up_ori, bound_sd_ori))


def cst_map_to_CEN_by_aa_type( file, outfile, seq, mapping_library_file=None, verbose=False, mode='aadep' ):
	def needs_mapping( atom ):
		return ((not atom in [ 'H', '1H', '2H', '3H', 'HA', '1HA', '2HA' ]) and ( atom.find('H') >=0 or atom.find('Q')>=0))

	if not mapping_library_file and 'csrosettaDir' in os.environ:
		mapping_library_file=os.environ['csrosettaDir']+"/database/CEN_to_proton_distances.txt"
	else:
		print mapping_library_file
		raise library.MissingInput("require a library of proton to CEN distances for mapping" )

	if isinstance(outfile, str):
      outfile=open(outfile,'w')

	print 'read CEN to proton distances from %s..'%mapping_library_file
	librarylist=open(mapping_library_file,'r').readlines()
#@zak use smarter data structures for this: dictionary for aa and dictionary for atom-names
	dist_map_library={}
	for m in librarylist:
		cols=m.split()
		if '#'==cols[0][0]: continue
		aa=cols[2]
		dist=(float(cols[5]), float(cols[4]) )
		atom=cols[3]
		if not aa in dist_map_library:
			dist_map_library[aa]={}
		dist_map_library[aa][atom]=dist

	for line in open(file,'r'):
		l=string.split(line)
		if not l[0] in ['AmbiguousNMRDistance','AtomPair' ]:
			outfile.write(line)
			continue
		res1=int(l[2])
		res2=int(l[4])
		atom1=l[1]
		atom2=l[3]
		bound_low = float( l[6] );
      try:
          bound_up = float( l[7] );
          bound_sd = float( l[8] );
      except:
          print '\n\n expected float in col 8 and 9 of:\n'
          print l[7],'\n',l
          exit()

		noe_dist=float(l[7])
		try:
			aa1=seq[res1-1]
			aa2=seq[res2-1]
		except IndexError:
			raise library.InconsistentInput( "fasta sequence is inconsistent with constraint file, cannot find residue %d"%max( res1, res2) )

		#translate into Rosetta atom names
		pdb_atom1=parse_NMR_name(atom1,aa1,res1)
		pdb_atom2=parse_NMR_name(atom2,aa2,res2)
		first_atom1=pdb_atom1[0]
		first_atom2=pdb_atom2[0]

		#do we have to map at all ?
		dist1=(0,0)
		dist2=(0,0)

		try:
			if needs_mapping( first_atom1 ):
				atom1 = 'CEN'
				dist1 = dist_map_library[aa1][first_atom1]

			if needs_mapping( first_atom2 ):
				atom2 = 'CEN'
				dist2 = dist_map_library[aa2][first_atom2]
		except:
			raise library.ProgramError("cannot find one of ( %s, %s ) and ( %s, %s ) in dist_map_library file %s"%( aa1, first_atom1, aa2, first_atom2, mapping_library_file) )
   		#now we have tuples (min, max) for both protons to their own CEN
			#now mode decides what to do with it.
		pad=0
		sd=0
		if needs_mapping( first_atom1 ) or needs_mapping( first_atom2 ):
			if mode=='aadep':
				pad=dist1[0]+dist2[0]
				sd=math.sqrt(dist1[1]-dist1[0]+dist2[1]-dist2[0])
			elif mode=='aadep_padonly':
				pad=dist1[0]+dist2[0]
				sd=0
			elif mode=='aadep_mid':
				pad=0.5*(dist1[0]+dist1[1]+dist2[0]+dist2[1])
				sd=0
			elif mode=='aadep_mid_sd':
				pad=0.5*(dist1[0]+dist1[1]+dist2[0]+dist2[1])
				sd=0.5*math.sqrt(dist1[1]-dist1[0]+dist2[1]-dist2[0])
			elif mode=='aadep_mid_sdfix':
				pad=0.5*(dist1[0]+dist1[1]+dist2[0]+dist2[1])
				sd=min(dist1[0],1)+min(dist2[0],1)  #since max-dist is > 1 this should basically be 1 or 0 for mapable atoms
			else:
				raise library.ProgramError("ma/home/olange/cs_targetlib/ilv_sgr145/rasrec/map_aadep/nmr_data/final_and_manual_noQF.cst.centroidpping mode %s is unknown"%mode )

		#write output
		cst_type='AtomPair'
      if ( len(atom1)>2 and not 'CEN' in atom1 ) or ( len(atom2)>2 and not 'CEN' in atom2 ):
			cst_type='AmbiguousNMRDistance'
		if ( res1==1 or res2==1 ):     #possibly H on resi 1
			cst_type='AmbiguousNMRDistance'

		bound_sd=bound_sd+sd
		bound_up=bound_up+pad
		outfile.write('%15s %5s %5d %5s %5d BOUNDED %5.3f %5.3f %5.3f #remark %5.3f pad %5.3f %5.3f\n'%(cst_type,atom1,res1,atom2,res2,bound_low,bound_up,bound_sd,bound_up-pad,pad,sd))


def parse_NMR_name( name, aa, pos ):
#returns list of atoms that correspond to the input name:
# HG11 ---> [ 1HG1 ]
# QG --> [ 1HG1, 2HG1, 3HG1 ]
# etc.
	atoms=[]
	if ( ( name[0:2] == "HD" ) and aa ==  'H' ):
		atoms.append( "CD2"  )
	if name[0:2] == "HE"  and aa =='H':
		atoms.append( "CE1"  )
	if ( name == "H" and pos == 1) :
		atoms.append( "1H"  )
		atoms.append( "2H"  )
		atoms.append( "3H"  )
	elif ( (name[0:2] == "QA" or name == "HA") and ( aa ==  'G' )) :
		atoms.append( "1HA"  )
		atoms.append( "2HA"  )
	elif ( ( name == "QB" or name == "HB" ) and aa in 'ALSNQPKCDERYFWHM' ) :
		atoms.append( "1HB"  )
		atoms.append( "2HB"  )
		if ( aa ==   'A' ):
			atoms.append( "3HB"  );
	elif ( ( name == "QD1" or  name == "HD1" ) and aa in 'IL' ) :
		atoms.append( "1HD1"  )
		atoms.append( "2HD1"  )
		atoms.append( "3HD1"  )
	elif ( ( name == "QD2" or name == "HD2" ) and aa in 'LN' ) :
		atoms.append( "1HD2"  )
		atoms.append( "2HD2"  )
		if ( aa ==  'L' ):
			atoms.append( "3HD2"  )
	elif ( ( name == "QQD" or name == "QD" ) and aa ==  'L' ) :
		atoms.append( "1HD2"  )
		atoms.append( "2HD2"  )
		atoms.append( "3HD2"  )
		atoms.append( "1HD1"  )
		atoms.append( "2HD1"  )
		atoms.append( "3HD1"  )
	elif ( ( name == "QD" or name == "HD" ) and aa in  'PKRYFH' ) :
		if ( aa ==  'R' or aa ==  'K' or aa ==  'P' ) :
			atoms.append( "1HD"  )
			atoms.append( "1HD"  )
		else : # 'Y', 'F', 'H' ('H' doesn't get down here... )
			atoms.append( "HD1"  )
			atoms.append( "HD2"  )
	elif ( ( name == "QE" or name == "HE" ) and aa in  'YFWMK') :
		if ( aa ==  'F' or aa ==  'Y' ) :
			atoms.append( "HE1"  )
			atoms.append( "HE2"  )
		elif ( aa ==  'W' ) :
			atoms.append( "HE1"  )
			atoms.append( "HE3"  )
		elif ( aa ==  'M' or aa ==  'K' ) : #'M' 'K'
			atoms.append( "1HE"  )
			atoms.append( "2HE"  )
			if ( aa ==  'M' ) :
				atoms.append( "3HE"  )
		 #QE 'M' 'K'
	elif ( ( name == "QG" or name == "HG" ) and aa in  'MQPKER' ) :
		atoms.append( "1HG"  )
		atoms.append( "2HG"  )
        #		atoms.append( "3HE"  )
	elif (( name == "QG2" or name == "HG2" ) and aa in  'ITV'  ) :
		atoms.append( "1HG2"  )
		atoms.append( "2HG2"  )
		atoms.append( "3HG2"  )
	elif ( ( name == "QQG" or name == "HG" ) and aa in  'IV'  ) :
		atoms.append( "1HG2"  )
		atoms.append( "2HG2"  )
		atoms.append( "3HG2"  )
		atoms.append( "1HG1"  )
		atoms.append( "2HG1"  )
		if ( aa !=  'I' ):
			atoms.append( "3HG1"  )
	elif (( name == "QG1" or name == "HG1" ) and aa in 'IV' ) :
		atoms.append( "1HG1"  )
		atoms.append( "2HG1"  )
		if ( aa !=  'I' ):
			atoms.append( "3HG1"  )
	elif (( name == "QE2" or name == "HE2" or name =="HE" ) and aa ==  'Q' ) :
		atoms.append( "1HE2"  )
		atoms.append( "2HE2"  )
	elif (( name == "HZ" or name == "QZ" ) and  aa in 'WK' ) :
		if ( aa ==  'K' ) :
			atoms.append( "1HZ"  )
			atoms.append( "2HZ"  )
			atoms.append( "3HZ"  )
		else : #aa=='W'
			atoms.append( "HZ2"  )
			atoms.append( "HZ3"  )
	elif ( name == "HB1" ) :
		atoms.append( "1HB"  )
	elif ( name == "HB2" ) :
		atoms.append( "2HB"  )
	elif ( name == "HB3" ) :
		if (  aa !=   'A' ) :
			atoms.append( "1HB"  )  #yeah they call it 2HB and 3HB...
		else :
			atoms.append( "3HB"  )
	elif ( name == "HD1" and not aa in 'FYW') :
		atoms.append( "1HD"  )
	elif ( name == "HD2" and not aa in 'FYW') :
		atoms.append( "2HD"  )
	elif ( name == "HD3" ) : #'K', 'P', 'R'  no other has HD3
		atoms.append( "1HD"  )
	elif ( name == "HG1" and aa !=  'T' ) :
		atoms.append( "1HG"  )
	elif ( name == "HG2" ) :
		atoms.append( "2HG"  )
	elif ( name == "HG3" ) : #'E', 'R', 'Q', 'M'
		atoms.append( "1HG"  )
	elif ( name == "HA1" ) :
		atoms.append( "1HA"  )
	elif ( name == "HA2" ) :
		atoms.append( "2HA"  )
	elif ( name == "HA3" ) : #'G'
		atoms.append( "1HA"  )
	elif (( name == "HE1" or name == "HE2" or name=="HE3" ) and not aa in 'FYW' ): #!='F' and aa!='Y' and aa!='W' ):
		if ( name == "HE3" and aa !=  'M' ) :
			atoms.append( "1HE"  ) #e.g. 'K'
		else :
			print ("has a problem here")
                #atoms.append( id::AtomID( pose.residue_type(res).atom_index(name.substr(2,1)+name[0:2])
	elif ( name == "HD11" ) :
		atoms.append( "1HD1"  )
	elif ( name == "HD12" ) :
		atoms.append( "2HD1"  )
	elif ( name == "HD13" ) :
		atoms.append( "3HD1"  )
	elif ( name == "HD21" ) :
		atoms.append( "1HD2"  )
	elif ( name == "HD22" ) :
		atoms.append( "2HD2"  )
	elif ( name == "HD23" ) :
		atoms.append( "3HD2"  )
	elif ( name == "HG11" ) :
		atoms.append( "1HG1"  )
	elif ( name == "HG12" ) :
		atoms.append( "2HG1"  )
	elif ( name == "HG13" ) :
		if ( aa ==  'I' ) :
			atoms.append( "1HG1"  )
		else :
			atoms.append( "3HG1"  )
	elif ( name == "HG21" ) :
		atoms.append( "1HG2"  )
	elif ( name == "HG22" ) :
		atoms.append( "2HG2"  )
	elif ( name == "HG23" ) :
		atoms.append( "3HG2"  )
	elif ( name == "HE11" ) :
		atoms.append( "1HE1"  )
	elif ( name == "HE12" ) :
		atoms.append( "2HE1"  )
	elif ( name == "HE13" ) :
		atoms.append( "3HE1"  )
	elif ( name == "HE21" ) :
		atoms.append( "1HE2"  )
	elif ( name == "HE22" ) :
		atoms.append( "2HE2"  )
	elif ( name == "HH11" ) :
		atoms.append( "1HH1"  )
	elif ( name == "HH12" ) :
		atoms.append( "2HH1"  )
	elif ( name == "HH21" ) :
		atoms.append( "1HH2"  )
	elif ( name == "HH22" ) :
		atoms.append( "2HH2"  )
	elif ( (name == "HH1" or name == "QH1" ) and aa ==  'R' ) :
		atoms.append( "1HH1"  )
		atoms.append( "2HH1"  )
	elif ( (name == "HH2" or name == "QH2" ) and aa ==  'R' ) :
		atoms.append( "1HH2"  )
		atoms.append( "2HH2"  )
	elif ( (name == "HH" or name == "QQH" ) and aa ==  'R' ) :
		atoms.append( "1HH1"  )
		atoms.append( "2HH1"  )
		atoms.append( "1HH2"  )
		atoms.append( "2HH2"  )
	else :
              #tr.Trace << "adding " << id:: ( name  ) << " "
        #				 << pose.residue_type( res ).name3() << std::endl
		atoms.append( name  )
          #tr.Trace << "as atom: " << atoms.back()

	return atoms;

class test():
	def __init__(self):
        #atomfile=open('atoms.txt')
        #atomlist=atomfile.readlines()
		atomlist=self.test_data()
		seq='IDLCLSSEGSEVILATSSDEKHPPENIIDGNPETFWTTTGMFPQEFIICFHKHVRIERLVIQSYFVQTLKIEKSTSKEPVDFEQWIEKDLVHTEGQLQNEEIVAHGSATYLRFIIVSAFDHFASVHSVSAEGTVVS';
		success=True
		for ct,r in enumerate(atomlist):
			tags=r.split()
			resid=int(tags[1])
			atoms=parse_NMR_name(tags[0],seq[resid-1],resid)
#			 print atoms
			success=success and self.check_results( ct, atoms, r )
		if success:
			print "UNIT TEST SUCCESSFUL"
		else:
			print "there were unit test failures!!!"
	def test_data(self):
		str=['H 1','H 98','H 99','HA2 106','HA2 132','HA2 30','HA3 106','HA3 132','HA3 95','HB 1','HB 118','HB2 11','HB2 110','HB3 73','HD21 26','HD21 31','HD21 99','HD22 26','HD22 31','HD22 99','HD3 112','HD3 70','HD3 73','HD3 88','HE1 36','HE1 85','HE21 67','HE21 84','HE21 96','HE22 67','HE22 84','HE22 96','HE3 36','HE3 85','HG 111','HG 14','HG 3','HG 59','HG 69','HG 90','HG 97','HG12 114','HG12 115','HG12 27','HG3 8','HG3 88','HH2 36','HH2 85','HZ 113','HZ 119','HZ 46','HZ 50','HZ 65','HZ 82','HZ2 36','HZ2 85','HZ3 85','QB 101','QB 105','QD 77','QD 82','QD1 111','QD1 14','QD1 3','QD2 69','QD2 90','QE 110','QE 119','QE 122','QE 21','QE 82','QE 88','QG 100','QG 131']
		return str

	def check_results(self, ct, result, input):
		res = [['1H', '2H', '3H'], ['H'], ['H'], ['2HA'], ['2HA'], ['2HA'], ['1HA'], ['1HA'], ['1HA'], ['HB'],  ['1HB', '2HB', '3HB'], ['2HB'], ['2HB'], ['1HB'], ['1HD2'], ['1HD2'], ['1HD2'], ['2HD2'], ['2HD2'], ['2HD2'], ['1HD'], ['1HD'], ['1HD'], ['1HD'], ['HE1'], ['HE1'], ['1HE2'], ['1HE2'], ['1HE2'], ['2HE2'], ['2HE2'], ['2HE2'], ['HE3'], ['HE3'], ['HG'], ['HG'], ['HG'], ['HG'], ['HG'], ['HG'], ['HG'], ['2HG1'], ['2HG1'], ['2HG1'], ['1HG'], ['1HG'], ['HH2'], ['HH2'], ['HZ'], ['HZ'], ['HZ'], ['HZ'], ['HZ'], ['HZ'], ['HZ2'], ['HZ2'], ['HZ3'], ['1HB', '2HB'], ['1HB', '2HB'], ['1HD', '1HD'], ['HD1', 'HD2'], ['1HD1', '2HD1', '3HD1'], ['1HD1', '2HD1', '3HD1'], ['1HD1', '2HD1', '3HD1'], ['1HD2', '2HD2', '3HD2'], ['1HD2', '2HD2', '3HD2'], ['HE1', 'HE2'], ['HE1', 'HE2'], ['HE1', 'HE2'], ['1HE', '2HE'], ['HE1', 'HE2'], ['1HE', '2HE'], ['1HG', '2HG'], ['1HG', '2HG']]
		if res[ct]!=result:
			print "UNIT-TEST FAILED: wrong result: input %s returned %s but should have returned %s"%(input, result, res[ct] )
			return False
		return True
#test()

