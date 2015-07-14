#!/usr/bin/env python2.7
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
## -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-

import argparse
from os.path import basename
import traceback, sys
import os, time
from os import path
import sets
import subprocess

delim_width=100
parser = argparse.ArgumentParser(prog=basename(__file__),description='Installation of CS-Rosetta toolbox')
parser.add_argument("-notalos", action='store_false', dest='install_talos', default=True, help='do not download NMRPipe (TALOS+) package' );
parser.add_argument('-nopicking', default=False, action='store_true', help='suppress installation of frag-picker, this takes a lot of space and is maybe not required (e.g., Laptop/Desktop installation vs. Mainframe installation)' )
parser.add_argument('-update', default=False, action='store_true', help='update genome and vall databases if updates available')
parser.add_argument('-nocheck', default=True, dest='check', action='store_false', help='do not check for availability of updates for genome and vall databases')
parser.add_argument("-talos_dir", default=None, help='where to find an existing TALOS+ installation, use this if autodetection fails' );
parser.add_argument('-rosetta', default=None, help='path to an existing ROSETTA package' )
parser.add_argument('-rosetta_database',default=None, help='path to rosetta_database (only needed if not in same location as rosetta_source')
parser.add_argument('-blast', default=None, help='specify path to pre-installed BLAST executable. If not given BLAST will be downloaded automatically')
parser.add_argument('-vall', default=None, help='specify path to pre-installed CS-Rosetta VALL' )

parser.add_argument('-genomes_db', default=None, help='provide path to genome database, expect to find files nr.00.phd, nr.00.phi, etc. at this location')
parser.add_argument('-blast_platform', default='x64-linux', help='platform string to select appropriate blast executable when downloading')
parser.add_argument("-target_lib", default='cs_targetlib', help='generate target entries at this path, if not absolute path this will be prefixed with $HOME' )
parser.add_argument('-talos_version',default='2012', help='force talos version (give the year: e.g., 2012 )',type=int)
args = parser.parse_args()

wd=os.getcwd()
final_msgs=sets.Set()
com = os.path.dirname(os.path.normpath(os.path.join(os.getcwd(), sys.argv[0])))
root= os.path.normpath(com.replace('com',''))
com = root+'/com'
hroot=root
home=None

class DBNotFound(Exception):
	pass
def check_db_path(db_path):
	if path.exists('%s/chemical/atom_type_sets/fa_standard/atom_properties.txt'%db_path):
		return db_path
	if path.exists('%s/rosetta_database/chemical/atom_type_sets/fa_standard/atom_properties.txt'%db_path):
		return '%s/rosetta_database'%db_path
	raise DBNotFound

def check_rosetta_version( rosetta_source, quiet=False ):
	rosetta_source=rosetta_source.replace('bin/','')
	if path.exists(rosetta_source+'/bin') and path.exists(rosetta_source+'/bin/minirosetta.default.linuxgccrelease'):
		if not quiet: print 'AUTO-DETECTION: check rosetta version in %s via minirosetta.default.linuxgccrelease'%rosetta_source
		cmd='%s/bin/minirosetta.default.linuxgccrelease -help | grep iterative:normalize | wc -l'%rosetta_source
		pipe=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True)
		line = pipe.stdout.readline()
		return int(line)>1
	if path.exists(rosetta_source+'/src/core'):
		if not quiet: print 'AUTO-DETECTION: check rosetta version in %s via source files'%rosetta_source
		return path.exists(rosetta_source+'/src/protocols/jd2/archive/NormalizedEvaluatedArchive.cc')
	return False


def center_text( text, width, fill=' ' ):
	t=text.strip()
	nfill=(width-len(t))/len(fill)
	return fill*(nfill/2-1)+t+fill*(nfill/2)
def silent_wget( link, dir, verbose=False,options="-N -c",progress=True ):
	import math
	def int_size_str(size):
		return int(size.replace('K','000').replace('M','000000'))
	def progressbar(skipped,current,full):
		N=100
		percentage=int(float(current)/full*100)
		s=''
		s+='+'*int(float(skipped)/full*N)
		s+='='*int(float(current-skipped)/full*N)
		s+='>'
		s+='.'*(int(float(full-current)/full*N)-2)
		return s,percentage
	print 'Execute: wget %(options)s %(link)s'%locals()
	pipe=subprocess.Popen('wget %(options)s %(link)s'%locals(), cwd=dir, shell=True, stderr=subprocess.PIPE )
	full=None
	skipped=0
#	skip=5
#	s=skip
	printed_progress=False
	downloaded=False
	header=''
	for line in pipe.stderr:
		try:
#			print '-->'+line[:-1]
			tags=line.split()
			if len(tags)==0: continue
			elif 'Length:' in tags:
				full=int(tags[tags.index('Length:')+1])
				continue
			elif 'skipping' in line:
				skipped=int_size_str(tags[tags.index('skipping')+1])
				continue
			elif 'saved' in line:
				continue
			elif not full: header+=line; continue
			elif 'RETR' in line or 'CWD' in line or 'PASV' in line: continue
			elif '=>' in line: continue
			else:
				try:
					current=int_size_str(tags[0])
#					if s>0: s-=1; continue
#					s=skip
					eta=tags[-1]
					speed=tags[-2]
					progress,percentage=progressbar(skipped,current,full)
					if progress: sys.stdout.write('\r %(percentage)d%% [%(progress)s] %(eta)s       '%locals())
					if progress: sys.stdout.flush()
					printed_progress=True
					downloaded=True
				except ValueError:
					if printed_progress: print; printed_progress=False
					print 'WGET: %s'%line[:-1]
				except KeyboardInterrupt:
					print
		except KeyboardInterrupt:
			raise
		except:
			print 'Exception in line %s'%line
			raise
	if printed_progress: print
	if not downloaded and not 'Remote file no newer' in header and not 'Wrote HTML-ized' in header and not 'nothing to do' in header:
			print header
try:
	home= os.path.realpath(os.environ['HOME'])
	pref=os.path.commonprefix([home,root])
	hroot=root;#.replace(pref,'$HOME')
except KeyError:
	pass

if args.target_lib[0]!='/':
	target_lib='$HOME/'+args.target_lib
else:
	target_lib = args.target_lib

talos_dir=None
if args.talos_dir:
	talos_dir=args.talos_dir
if talos_dir and talos_dir[0]!='/':
	talos_dir=os.path.normpath(os.path.abspath(wd+'/'+talos_dir))

rosetta_dir=None
rosetta_bin=None
if args.rosetta:
	rosetta_dir=args.rosetta
if rosetta_dir and rosetta_dir[0]!='/':
	rosetta_dir=wd+rosetta_dir
if not rosetta_dir:
	paths=os.environ['PATH'].split(':')
	rosetta_paths=sets.Set()
	for p in paths:
		rosetta_paths.add(p.rstrip('/'))
	print 'AUTO-DETECTION: check  $PATH'
	rosetta_paths=[ x for x in rosetta_paths if check_rosetta_version( x.replace('/bin','') ) ]
	if len(rosetta_paths)>1:
		print 'AUTO-DETECTION: found rosetta in $PATH'
		pp=rosetta_paths.pop()
		rosetta_dir=pp.replace('/rosetta_source/bin','')
		rosetta_bin=pp

if not rosetta_dir:
	print 'AUTO-DETECTION: look in $HOME/rosetta'
	if check_rosetta_version( home+'/rosetta/rosetta_source' ):
		rosetta_dir=home+'/rosetta'
else:
	rosetta_dir=rosetta_dir.replace('rosetta_source/bin','').replace('rosetta_source','').rstrip('/')

#final paranoia check: quiet because it is probably double
if rosetta_dir:
	if not check_rosetta_version( rosetta_dir+'/rosetta_source', quiet=True ):
		rosetta_dir=None
		rosetta_bin=None
	else:
		if not rosetta_bin:
			rosetta_bin='%(rosetta_dir)s/rosetta_source/bin'%locals()

if rosetta_bin:
	if not check_rosetta_version( rosetta_bin.replace('/bin',''), quiet=True ):
		rosetta_bin=None

if not args.rosetta and rosetta_bin:
	print 'AUTO_DETECTED: ROSETTA bin directory at '+rosetta_bin

rosetta_db_path=None
db_source=None
if args.rosetta_database:
	try:
		rosetta_db_path=check_db_path(args.rosetta_database)
		db_source='cmd-line'
	except DBNotFound:
		exit('Cannot find the ROSETTA Database directory at the provided path %s'%args.rosetta_database)
elif rosetta_bin:
	top_level_dir=rosetta_bin.replace('/bin','')
	while not rosetta_db_path and len(top_level_dir)>3:
		try:
			rosetta_db_path=check_db_path(top_level_dir)
			db_source='closest top-level directory from %(rosetta_bin)s: %(top_level_dir)s'%locals()
		except:
			top_level_dir=path.dirname(top_level_dir)
			pass
	try:
		rosetta_db_path=check_db_path(os.environ['ROSETTA3_DB'])
		db_source='ROSETTA3_DB environment variable'
	except:
		pass

if rosetta_db_path and not args.rosetta_database:
	print 'AUTO_DETECTED: ROSETTA database directory at %(rosetta_db_path)s\n               via %(db_source)s'%locals()

if not rosetta_bin:
	msg= 'AUTO-DETECTION FAILED: cannot find ROSETTA binaries.\n'\
	     '                     Run with -rosetta to point to binaries if toolbox is intended for automatic setup of CS-ROSETTA calculations'
  print msg
  final_msgs.add(msg)

if not rosetta_db_path:
	msg= 'AUTO-DETECTION FAILED: cannot find ROSETTA database.\n'\
	     '                     Run with -rosetta_database to point to database if toolbox is intended for automatic setup of CS-ROSETTA calculations'
  print msg
  final_msgs.add(msg)


print "Installation of CS-Rosetta Toolbox in path %s ..."%root

print "Generate init script..."
if rosetta_bin:
	path_extension="${csrosettaCom}:%(rosetta_bin)s"%locals()
else:
	path_extension="${csrosettaCom}"%locals()

init_text="""#!/bin/env bash
#
# CS-ROSETTA: System for Chemical Shifts based protein structure prediction using ROSETTA
#

#path names to find the CS-Rosetta toolbox files
export csrosettaDir=%(hroot)s
export csrosettaCom=${csrosettaDir}/com_alias

#provide CS-Rosetta applications in system-path
export PATH=${PATH}:%(path_extension)s

#make sure that python modules can be found
#CSHELL: if ( ! ($?PYTHONPATH) ) then
#CSHELL:    export PYTHONPATH
#CSHELL: endif
export PYTHONPATH=${PYTHONPATH}:${csrosettaDir}/python

#this is useful for benchmarking and used by the scripts setup_target.py display_target.py and setup_run.py
export CS3_BENCH_TARGETLIB=%(target_lib)s

#the path to the ROSETTA installation (yields defaults for setup_run options -database and -binaries_prefix)
"""%locals()
if rosetta_bin and rosetta_dir:
	init_text+="""
	export ROSETTA3_PATH=%(rosetta_dir)s
	export ROSETTA3_BIN=%(rosetta_bin)s
"""%locals()
#{'hroot':hroot,'target_lib':target_lib,'rosetta_bin':rosetta_bin,'rosetta_root':rosetta_dir})

if rosetta_db_path:
	init_text+='export ROSETTA3_DB=%(rosetta_db_path)s\n'%locals()

#make com_aliases
import subprocess
for ext in ['py','com']:
	subprocess.call('cd %s; mkdir -p com_alias; cd com_alias;'%root+
								'for i in $( ls ../com/*.%s ); do ln -sf $i $( basename $i .%s ); done'%(ext,ext), shell=True)

#add TALOS to init-script if possible...
if not talos_dir:
  #try autodection of TALOS+
	keys=['talospDir','TALOS_DIR','TALOSP_DIR','NOT_FOUND']
	for k in keys:
		if k in os.environ:
			break
	if k != 'NOT_FOUND':
		talos_dir=os.path.normpath(os.path.realpath(os.environ[k]))
		if path.exists(talos_dir+'/bin/TALOS+.linux'):
			print 'AutoDetected: TALOS: %s in environment variable %s'%(talos_dir, k)
		else:
			talos_dir=None

sparta_dir=None
promega_dir=None
if not talos_dir:
	if not args.install_talos:
		print "\n","*"*delim_width
		print """       WARNING         WARNING         WARNING         WARNING
TALOS+ is essential for many functions of CS-Rosetta, if not installed
only limited functionality is available.
Please provide TALOS+ root directory with option -talos_dir"""
		print '*'*delim_width
	else:
		print '='*delim_width
		print center_text('Install TALOS+ from NMRPipe (skip with -notalos)',delim_width)
		print '='*delim_width
		nmrPipe_dir=root+'/NMRPipe/'
		talos_link='http://spin.niddk.nih.gov/NMRPipe/install/download/talos.tZ'
		talos_file=nmrPipe_dir+'/'+path.basename(talos_link)
		stat_old=None
		if path.exists(talos_file):
			talos_dir=nmrPipe_dir+'/talosplus'
			sparta_dir=nmrPipe_dir+'/spartaplus'
			promega_dir=nmrPipe_dir+'/NMRPipe/promega'
			stat_old=os.stat(talos_file)
		attempts=2
		while attempts>0 and (args.check or not stat_old):
			attempts-=1
			subprocess.call('mkdir -p %(nmrPipe_dir)s'%locals(),shell=True)
			#cd %(nmrPipe_dir)s;
			print "Downloading from NIH website... "
			silent_wget(talos_link,nmrPipe_dir)
			print 'finished download %(talos_link)s'%locals()
			stat_new=os.stat(talos_file)
			unpack=False
			if stat_old!=stat_new:
				if stat_old: print 'downloaded talos.tZ is different than older file'
				unpack=True
			if not unpack and not path.exists('%(nmrPipe_dir)s/talosplus/bin/TALOS+.linux'%locals()):
				unpack=True
				print '%(nmrPipe_dir)s/talosplus/bin/TALOS+.linux not found, unpacking'%locals()
			if unpack:
				print 'unpacking NMRPipe-TALOS modules...'
				try:
					tarout=subprocess.check_output('cd %(root)s/NMRPipe; tar -xvzf %(talos_file)s'%locals(), shell=True)
					talos_dir=nmrPipe_dir+'/talosplus'
					subprocess.call('ln -sf ../NMRPipe/spartaplus/bin/SPARTA+.linux %(root)s/com_alias/sparta+'%locals(), shell=True)
					subprocess.call('ln -sf ../NMRPipe/promega/bin/PROMEGA.linux %(root)s/com_alias/promega'%locals(), shell=True)
					sparta_dir=root+'/NMRPipe/spartaplus'
					promega_dir=root+'/NMRPipe/promega'
					break
				except subprocess.CalledProcessError as e:
					print e
					print 'something wrong with unpacking that archive... remove file %(talos_file)s and start-over'%locals()
					subprocess.call('rm %s'%nr_file, shell=True)
					continue

talos_version='post_2012'
if talos_dir:
	print 'Checking Version of TALOS+...'
	pipe=subprocess.Popen('%s/bin/TALOS+.linux'%talos_dir,shell=True,stderr=subprocess.PIPE)
	for line in pipe.stderr:
		if 'Version' in line:
			print line
			years=['2009','2012']
			talos_version=args.talos_version
			for year in years:
				if year in line:
					print 'found year %s in TALOS output'%year
					talos_version=int(year)
	if talos_version>=2012:
		talos_version='post_2012'
	else:
		talos_version='pre_2012'

talos_env_name='TALOSP_DIR'
if talos_version=='pre_2012':
	talos_env_name='TALOS_DIR'
if talos_dir:
	init_text=init_text+"""
#Making sure that TALOS+ commands are available, too
export %(talos_env_name)s=%(talos_dir)s
"""%locals()

localized_talos=talos_dir.replace('%(root)s','..')
subprocess.call('ln -sf %(localized_talos)s/bin/TALOS+.linux %(root)s/com_alias/talos+'%locals(), shell=True)

if sparta_dir:
	init_text=init_text+"""
#Sparta+ installation
export SPARTAP_DIR=%s
"""%sparta_dir

if promega_dir:
	init_text=init_text+"""
#Sparta+ installation
export PROMEGA_DIR=%s
"""%promega_dir


#write init script
open(com+'/init.bashrc','w').write(init_text)

#write cshrc version
open(com+'/init.cshrc','w').write(init_text.replace('export','setenv').replace('=',' ').replace('bash','tcsh -f').replace('#CSHELL:',''))

#figure out login shell
try:
	shell= os.path.realpath(os.environ['SHELL'])
	print 'detected shell: ',shell
	wd=os.getcwd()
	os.chdir(com)
	if 'bash' in shell:
		link='init.bashrc'
	else:
		link='init.cshrc'

	try:
		os.symlink(link,'init')
		print 'created symlink %(root)s/com/%(link)s --> %(root)s/com/init'%locals()
	except OSError as exc:
		if exc.errno==17: #file exists (link exists already)
			os.remove('init')
			os.symlink(link,'init')
			print 'created symlink %(root)s/com/%(link)s --> %(root)s/com/init'%locals()
		else:
			raise

	os.chdir(wd)
except:
	print 'cannot detect input shell...'
#######################################
#######################################
#frag_picking
if args.nopicking: exit()
print '\n\n','* * '*(delim_width/4)
print center_text('Installing Fragment Picker',delim_width)
print '='*delim_width
print "For the fragment picker the installer will now download: Vall, BLAST and Genome NR database"
print 'to skip this step run with -nopicking or press ^C now'

print 'If modules are already downloaded the installer will check if updates are available'
print 'skip check with -nocheck and run with -update to download new files'
print '='*delim_width
print '* * '*(delim_width/4)
print
#figure out location of VALL
#is alredy installed?
import glob
import fnmatch
import shutil
import subprocess

def make_file_dict(files):
	dict={}
	for f in files:
		dict[path.basename(f)]=f
	return dict


#
if args.vall:
	vall_path=path.abspath(args.vall)
else:
	print '='*delim_width
	print center_text('               Install Vall',delim_width)
	print '='*delim_width

	installed=None
	installed_year=0
	vall_path=root+'/frag_picker/csrosetta_vall'
	vall_link='http://www.csr-dev.org/sites/default/files/downloads/csrosetta_vall_current.tgz'
	vall_archive=vall_path+'/'+path.basename(vall_link)
	md5_file_base=path.basename(vall_link)+'.md5'
	update_required=False
	test_dir=vall_path+'/test_md5'
	md5_file=vall_path+'/'+md5_file_base
	new_md5_file=test_dir+'/'+md5_file_base

	subprocess.call('mkdir -p %s'%vall_path, shell=True);
	if path.exists(root+'/frag_picker/blosum62.qij'):
		shutil.copy(root+'/frag_picker/blosum62.qij',root+'/frag_picker/csrosetta_vall/')
	else:
		if not path.exists(root+'/frag_picker/csrosetta_vall/blosum62.qij'):
			print 'ERROR: cannot find file %s - which should have been present in tarball of toolbox'%(root+'/frag_picker/blosum62.qij')
			exit()

	if not path.exists(vall_path):
		update_required=True
	else:
		installed=glob.glob(vall_path+'/vall.dat*')
		if not len(installed):
			update_required=True
		else:
			installed=installed[0].replace(root+'/','')
			installed_year=int(path.basename(installed).split('.')[-3])
			print 'installed vall found: %s -- issue year: %d'%(installed, installed_year)

			if not path.exists(md5_file):
				print 'no md5 file -- update required'
				update_required=args.update
				final_msgs.add('Chemical Shift VALL is possibly corrupted, rerun with -update')
			elif args.check:
				print 'Check for update on www.csrosetta.org: get server-md5-file for comparison....'
				subprocess.call('mkdir -p %s'%test_dir, shell=True);
	#	silent_wget(vall_link,test_dir)
	#	subprocess.call('cd %s; wget -N -c %s.md5'%(test_dir,vall_link), shell=True )
				silent_wget(vall_link+'.md5',test_dir,options='-N',progress=False)
				sys.stdout.write('check %(md5_file_base)s for update ...'%locals())
				if open(md5_file,'r').readline().split()[0]!=open(new_md5_file,'r').readline().split()[0]:
					print 'FAIL'
	#		print '\n\n'+'='*delim_width+'\ndetected difference in MD5 checksum between\n%s\nand %s\n'%(new_md5_file,md5_file)
	#		print '='*delim_width+'\n\n'
					update_required=args.update
					if update_required:
						print 'Will update %s from www.csrosetta.org server'%path.splitext(new_md5_file)[0]
					else:
						print 'Run with -update to retrieve newest vall database'
						final_msgs.add('Chemical Shift VALL outdated')
						final_msgs.add('Run with -update to retrieve newest version(s)')
				else:
					print 'Checksum OK: No update required.'
			else:
				print 'not checking for updates (-nocheck)...'

	if update_required:
		attempts=2
		while attempts>0:
			attempts-=1
			#subprocess.call('cd %(vall_path)s; wget -N -c %(vall_link)s'%locals(), shell=True )
			silent_wget(vall_link,vall_path)
			print 'finished download %(vall_link)s'%locals()
			print 'unpacking...'
			try:
				subprocess.check_call('cd %s;tar -xvzf %s'%(path.dirname(vall_path),path.basename(vall_path)+'/'+path.basename(vall_link)), shell=True )
				if path.exists(new_md5_file):
					print 'copy %(md5_file)s to %(new_md5_file)s'%locals()
					shutil.copy(new_md5_file,md5_file)
				else:
					#subprocess.call('cd %(vall_path)s; wget -N -c %(vall_link)s.md5'%locals(), shell=True )
					silent_wget(vall_link+'.md5',vall_path,options='-N')
				print 'file is unpacked, remove tar-archive now to save disc space...'
				subprocess.check_call('rm -f %(vall_archive)s'%locals(), shell=True)
				break
			except subprocess.CalledProcessError as e:
				print e
				print 'something wrong with unpacking that archive... remove file %(talos_file)s and start-over'%locals()
				subprocess.call('rm %(vall_archive)s'%locals(), shell=True)
				continue
	print "done installing VALL\n"

vall_dat=glob.glob(vall_path+'/vall.dat.*vCS')
if not len(vall_dat)==1:
	exit('did not find any or no unique match for vall.dat.*.vCS and vall.blast.* at specified vall location: %s'%vall_path)
vall_dat=vall_dat[0]
vall_blast=vall_dat.replace('vall.dat.','vall.blast.').replace('.vCS','')
if not path.exists( vall_blast ):
	exit('failed to find the corresponding blast vall file: %s to vall %s'%(vall_blast,vall_dat))

print '\n'+'='*delim_width
print center_text('               Install Blast',delim_width)
print '='*delim_width

if args.blast:
	print 'BLAST executable was overridden by commandline.\nUse %s instead of standard executable in csrosetta3-tree. Good luck... \n'%args.blast
	if 'blastpgp' in path.basename(args.blast):
		blast=args.blast
	else:
		blast=args.blast+'/bin/blastpgp'
	if blast[0]=='/':
		blast_path=blast
	else:
		blast_path=os.path.normpath(os.path.abspath(wd+'/'+blast))

else:
	blast_tgz=glob.glob(root+'/frag_picker/blast-*%s*.tar.gz'%args.blast_platform)
	if len(blast_tgz)==0:
		print '-'*delim_width
		print 'cannot find tarball for BLAST executable. Attempt to download...\n '
		print 'will download executable with platform-key "%s" (-blast_platform)'%args.blast_platform
		print 'if this is wrong platform download manually and copy into directory %s/frag_picker, or provide correct platform key with -blast_platform'%root
		print 'check here for list of files: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/'
		print '-'*delim_width+'\n'
	#	subprocess.call('cd %s/frag_picker;'%root+'wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/ -O blast_index.html;'+
		silent_wget('ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/','%s/frag_picker'%root,options=' -O blast_index.html')
		subprocess.call('grep /blast- blast_index.html | grep %s > blast_url'%args.blast_platform,cwd='%s/frag_picker'%root, shell=True)
		silent_wget(' blast_url','%s/frag_picker'%root,options=' -F -i')
		subprocess.call('rm blast_url blast_index.html',
										cwd='%s/frag_picker'%root, shell=True)
		blast_tgz=glob.glob(root+'/frag_picker/blast-*%s*.tar.gz'%args.blast_platform)
		print '\nfinished downloading.\n'+'-'*delim_width
		if len(blast_tgz)==0:
			print 'ERROR: cannot find tarball for legacy BLAST executable and failed to download... '
			print 'you should be able to find the executable here: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/'
			exit()
	if len(blast_tgz)>1:
		blast_tgz.sort()
		print 'Found multiple choices for blast tarball...'
		for i in blast_tgz:
			print i
		print 'will now use last: %s, if this is wrong choice delete all alternatives and run install.py again'%blast_tgz[-1]
	blast_tgz=blast_tgz[-1]
	print 'Install BLAST from %s'%blast_tgz
	os.chdir(root+'/frag_picker')
  if not path.exists( 'blast/bin/blastpgp' ):
		tar_out=subprocess.check_output('tar -xvzf %s'%blast_tgz, shell=True)
#	print tar_out
		dir=tar_out.split('\n')[0]
		print 'installed into dir %s'%dir
		try:
			os.symlink(dir,'blast')
			print 'created symlink frag_picker/blast'
		except OSError as exc:
			if exc.errno==17: #file exists (link exists already)
				print 'symlink existed, overwrite'
				os.remove('blast')
				os.symlink(dir,'blast')
		  else:
				raise
	blast_path='$frag_dir/blast/bin/blastpgp'
	os.chdir(wd)

print '\n'+'='*delim_width
print center_text('               Install Genome Database',delim_width)
print '='*delim_width


print "Installing Genome-Database which is required for frag-picking..."
genomes_path=None
if args.genomes_db:
	if args.genomes_db[0]!='/':
		genomes_path=os.path.normpath(os.path.abspath(wd+'/'+args.genomes_db))+'/'
	else:
		genomes_path=os.path.normpath(args.genomes_db+'/')+'/'
	if not path.exists(genomes_path+'nr.00.phd'):
		print 'ERROR: Cannot find NR database at specified location %s'%genomes_path
		exit()
	else:
		print 'Using existing NR genome database at location %s'%genomes_path
else:
	genomes_path=root+'/frag_picker/genomes/'
	download=False
	if path.exists(genomes_path):
		print 'Found existing NR genome database at location %s'%genomes_path
	#	genome_files=glob.glob(genomes_path+'nr.??.tar.gz')
	#	if len(genome_files)>=6:
	#		print 'Found compressed genome/nr database files...'
	#		for i in genome_files:
	#			print i
	#		print 'will simply use these files, force update with -force_genome_update (takes several hours)'
	#	else:
	#		download=True
	#else:
	if not args.check:
		print 'not checking for updates (-nocheck) ...'
	else:
		print '-'*delim_width
		print 'Check if updated genome database is available at ftp://ftp.ncbi.nlm.nih.gov/blast/db/...'
		print '-'*delim_width
		test_dir='frag_picker/genomes/genomes_test'
		subprocess.call('mkdir -p %s'%test_dir, shell=True);
#check that files are recent
#	new_md5_files=glob.glob('%s/nr.??.tar.gz.md5'%test_dir)
#	now=time.time()
#	for f in new_md5_files:
#		delta=now-path.getmtime(f)
#		if delta>86400: #more than 24h
#			print 'file %s older than %d days... removing it.'%(f,int(delta/86400))
#			subprocess.call('rm -f %s'%f, shell=True);

#	subprocess.call('cd %s; wget -N -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz.md5'%test_dir, shell=True )
		silent_wget('ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz.md5',test_dir,options='-N')

		md5_files=glob.glob(genomes_path+'nr.??.tar.gz.md5')
		new_md5_files=glob.glob('%s/nr.??.tar.gz.md5'%test_dir)

		md5_files=make_file_dict(md5_files)
		new_md5_files=make_file_dict(new_md5_files)
		for i in new_md5_files.keys():
			update_required=False
			if i not in md5_files:
				print "don't have corresponding md5 file for previously installed NR."\
						  "Downloading archive %s"%new_md5_files[i].replace('.md5','').replace('genomes_test/','')
				update_required=True
			else:
				sys.stdout.write('check %s for update...'%path.basename(new_md5_files[i]))
				if open(md5_files[i],'r').readline().split()[0]==open(new_md5_files[i],'r').readline().split()[0]:
					print 'OK'
				else:
					print 'FAIL'
#\n\n'+'='*80+'\ndetected difference in MD5 checksum between\n%s\nand %s\n'%(new_md5_files[i],md5_files[i])
#					print '='*80+'\n\n'
					update_required=args.update
					if update_required:
						print 'Will update %s from ftp-server'%path.splitext(new_md5_files[i])[0]
					else:
						final_msgs.add('Genome NR file: %s outdated'%path.basename(new_md5_files[i]).replace('.md5',''))
						final_msgs.add('Run with -update to retrieve newest version(s)')
			if update_required:
#				os.chdir(genomes_path)
				nr_file=path.splitext(path.basename( new_md5_files[i] ))[0]
				attempts=2
				while attempts>0:
				#subprocess.call('wget -N -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/%s'%nr_file, cwd=genomes_path, shell=True )
					silent_wget('ftp://ftp.ncbi.nlm.nih.gov/blast/db/%s'%nr_file,genomes_path )
#				pipe=subprocess.Popen('wget -N -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/%s'%nr_file, cwd=genomes_path, shell=True, stderr=subprocess.PIPE )
#				for i in pipe.stderr:
#					print '--->'+i[:-1]
					print 'finished download'
					try:
						sys.stdout.write('unpack ...(can take several minutes)');sys.stdout.flush()
						tar_out=subprocess.check_output('tar -xvzf %s'%nr_file, shell=True, cwd=genomes_path )
						subprocess.check_call('cp -f %s/%s/%s .'%(wd,test_dir,nr_file+'.md5'),cwd=genomes_path,shell=True )
						print 'file is unpacked, remove tar-archive now to save disc space...'
						subprocess.check_call('rm -f %s'%nr_file, cwd=genomes_path, shell=True)
						break
					except subprocess.CalledProcessError as e:
						print e
						print 'something wrong with unpacking that archive... remove file %s and start-over'%nr_file
						subprocess.call('rm -fr %s'%nr_file, cwd=genomes_path, shell=True)
						attempts=attempts-1
						continue

#			os.chdir(wd)

genomes_path=genomes_path.replace('%(hroot)s/frag_picker'%locals(),'$frag_dir')
blast_path=blast_path.replace('%(hroot)s/frag_picker'%locals(),'$frag_dir').replace('%(hroot)s'%locals(),'$base_dir')
if not rosetta_bin:
	print 'Cannot install frag-picker without valid ROSETTA installation'
	list=[m for m in final_msgs]
	list.sort()
	for m in list:
		#	if not m in printed:
		print m
	exit(1)
rosetta_bin=rosetta_bin.replace('%(rosetta_dir)s'%locals(),'$rosetta')
rosetta_db_path=rosetta_db_path.replace('%(rosetta_dir)s'%locals(),'$rosetta')

vall_dat=vall_dat.replace('%(hroot)s/frag_picker'%locals(),'$frag_dir')
vall_blast=vall_blast.replace('%(hroot)s/frag_picker'%locals(),'$frag_dir')

text="""
sub add_paths(%%) {

my $home="$ENV{HOME}";
my $base_dir="%(hroot)s";
my $frag_dir="$base_dir/frag_picker";
my $rosetta="%(rosetta_dir)s";

${$_[0]}{talos_version} = "%(talos_version)s";
${$_[0]}{vall}          = "%(vall_dat)s";
${$_[0]}{blastpgp}      = "%(blast_path)s";
${$_[0]}{blast_db}      = "%(genomes_path)s/nr";
${$_[0]}{vall_blast_db} = "%(vall_blast)s";

${$_[0]}{picker}        = "%(rosetta_bin)s/fragment_picker.default.linuxgccrelease";
${$_[0]}{mini_db}       = "%(rosetta_db_path)s";
}

1;
"""%locals()
open(root+'/frag_picker/setup_paths.pl','w').write(text)


#printed=sets.Set()
print "-"*delim_width+'\nInstallation Finished Successfully.'
list=[m for m in final_msgs]
list.sort()
for m in list:
#	if not m in printed:
		print m
#		printed.add(m)


# have now a python wrapper for pick_fragments
#try:
#	os.symlink('../frag_picker/pick_fragments.pl','pick_fragments')
#except OSError as exc:
#	if exc.errno==17: #file exists (link exists already)
#		print 'symlink existed, overwrite'
#	else:
#		raise



#matches=glob.glob('%s/csrosetta_vall*.tgz'%root)
# #matches=glob.glob('../csrosetta_vall*.tgz')
# matches=glob.glob('%s/../csrosetta_vall*.tgz'%root)
# for sroot, dirnames, filenames in os.walk(root):
#   for filename in fnmatch.filter(filenames, 'csrosetta_vall*.tgz'):
# 		matches.append(os.path.join(root,sroot, filename))
# if len(matches)>0:
# 	print "Have found the following vall possible for installation..."
# 	issue_years=[]
# 	best_year=installed_year
# 	for i in matches:
# #		print path.basename(i).split('.')
# 		try:
# 			issue_years.append(int(path.basename(i).split('.')[-3]))
# 			if issue_years[-1]>best_year:
# 				best_year=issue_years[-1]
# 				best_vall=i
# 		except:
# 			issue_years.append(0)

# 		print "%30s -- issue year %d"%(i.replace(root+'/',''),issue_years[-1] )

# print '='*80
# print '               Install Vall'
# print '='*80

# if best_vall:
# 	print "going to install %s..."%best_vall
# 	os.chdir(root+'/frag_picker')
# 	if installed:
# 		print 'backing up old VALL directory...'
# 		shutil.move('csrosetta_vall','csrosetta_vall_OUTDATED')
# 	subprocess.call('tar -xvzf %s'%best_vall, shell=True )
# 	print "done installing VALL\n"
# 	os.chdir(wd)
# elif installed_year>0:
# 	print "up-to-date vall is present will keep installed vall: %s"%installed
# else:
# 	print "ERROR: Cannot install fragment picker - Vall required! \n  - download from http://www.csrosetta.org/downloads"
# 	exit()
