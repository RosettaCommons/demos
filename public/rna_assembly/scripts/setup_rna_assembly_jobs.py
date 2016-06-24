#!/usr/bin/python
# Script for creating RNA helices and inter-helix junctions/loop ('motifs')
# and then assembling them together.
# (C) Rhiju Das 2008-2011

import string
from sys import argv,exit
import os
from os import system, getcwd
from os.path import exists, expanduser, abspath, dirname

if len( argv ) < 2:
    print argv[0], ' <sequence.fasta>  <secstruct file> [<native.pdb> <constraints.cst> <torsions.torsions>]'
    exit( 0 )

#EXE_DIR =  expanduser('~')+'/src/rosetta_TRUNK/rosetta_source/bin/'
scripts_path = dirname( abspath( argv[0] ) ) + '/'

if os.environ.has_key("ROSETTA3"):
	EXE_DIR = os.environ.get("ROSETTA3")+"/bin/"
else:
	EXE_DIR =  scripts_path +'/../../../main/source/bin/'
print EXE_DIR
EXE_DIR = abspath( EXE_DIR )

# for pdbslice.py
if os.environ.has_key("ROSETTA_TOOLS"):
	tools_scripts_path = os.environ.get("ROSETTA_TOOLS")+"/rna_tools/bin"
else:
	tools_scripts_path = abspath( dirname( abspath( argv[0] ) ) + '/../../../../tools/rna_tools/bin' )

if not exists( EXE_DIR ):
    print 'Need to set EXE_DIR in '+argv[0]+' to match an existing directory'
    exit( 0 )

fasta_file = argv[1]
params_file = argv[2]

native_exists = 0
data_exists = 0
cst_exists = 0
torsions_exists = 0
EXE_extension = ""
for i in range( 3, len( argv )):
	if argv[i][-4:] == '.pdb':
		native_pdb_file = argv[i]
		native_exists = 1
	if argv[i][-5] == '.data':
		data_file = argv[i]
		data_exists = 1
	if argv[i][-4:] == '.cst':
		cst_file = argv[i]
		cst_exists = 1
	if argv[i][-9:] == '.torsions':
		torsions_file = argv[i]
		torsions_exists = 1
	if argv[i][:14] == 'exe_extension=':
		EXE_extension = argv[i][14:]
		rna_helix_exe = EXE_DIR + '/rna_helix'+EXE_extension

if EXE_extension == '':
	EXE_extension = '.linuxgccrelease'
	rna_helix_exe = EXE_DIR + '/rna_helix'+EXE_extension
	if not exists( rna_helix_exe ):
		EXE_extension = '.macosgccrelease'
		rna_helix_exe = EXE_DIR + '/rna_helix'+EXE_extension

if not exists( rna_helix_exe ):
	print 'Cannot find '+rna_helix_exe+'.  Different extension than .linuxgccrelease or macosgccrelease?'
	exit( 0 )

# Read in files
lines = open( fasta_file ).readlines()
sequence = lines[1][:-1]
print sequence
numres = len( sequence )

# Read in data information
data_info = []
if data_exists:
    backbone_burial_info = []
    lines = open( data_file ).readlines()
    for line in lines:
        if len( line ) > 6 and line[:6]=='EXPOSE':
            cols = string.split( line )
            for i in range( len(cols)/3 ):
                data_info.append( [int( cols[ 3*i+1] ), cols[3*i+2],cols[3*i+3]] )
        if len( line ) > 15 and line[:15]=='BACKBONE_BURIAL':
            cols = string.split( line )
            for i in range( 1,len(cols) ):
                backbone_burial_info.append( int( cols[i] ) )

pyrimidines = ['c','u']
purines = ['a','g']
cst_info = []
if cst_exists:
    lines = open( cst_file ).readlines()
    for line in lines:
        if len( line ) > 6 and line[0]!='[':
            cols = string.split( line )
            atom_name1 = cols[0]
            res1 = int( cols[1] )
            atom_name2 = cols[2]
            res2 = int( cols[3] )
            if ( sequence[ res1 - 1 ] in pyrimidines and atom_name1=='N1'):
                atom_name1 = 'N3'
                print 'correcting atom name for ', res1
            if ( sequence[ res2 - 1 ] in pyrimidines and atom_name2=='N1'):
                atom_name2 = 'N3'
                print 'correcting atom name for ', res2
            if ( sequence[ res1 - 1 ] in purines and atom_name1=='N3'):
                atom_name1 = 'N1'
                print 'correcting atom name for ', res1
            if ( sequence[ res2 - 1 ] in purines and atom_name2=='N3'):
                atom_name2 = 'N1'
                print 'correcting atom name for ', res2

            cst_info.append( [ atom_name1, res1, atom_name2, res2, string.join( cols[4:] )] )

pair_map = {}
all_pairs = []

complement = {'a':['u'], 'u':['a','g'], 'c':['g'], 'g':['c','u']};

stem_out_files = []
motif_out_files = []

# Parse out stems
lines = open( params_file ).readlines()
cutpoints_original = []
for line in lines:
    if line[:4] == 'STEM':
        cols = string.split( line )
        for i in range( len( cols )):
            if cols[i] == 'PAIR':
                #Offset to get to python numbering (starts with zero)
                res1 = int(cols[i+1])-1
                res2 = int(cols[i+2])-1
                pair_map[ res1 ] = res2
                pair_map[ res2 ] = res1
                all_pairs.append( [res1,res2] )
                assert ( sequence[res1] in complement[ sequence[res2] ] )
    elif line.count( '(' ) > 0:         #maybe dot/bracket notation (((...)))
        print line
        left_brackets = []
        for i in range( len(line) ):
            if line[i] == '(':  left_brackets.append( i )
            if line[i] == ')':
                res1 = left_brackets[-1]
                res2 = i
                del( left_brackets[-1] )
                pair_map[ res1 ] = res2
                pair_map[ res2 ] = res1
                all_pairs.append( [res1,res2] )
                assert ( sequence[res1] in complement[ sequence[res2] ] )
        assert( len (left_brackets) == 0 )
    else:
        try:
            cols = string.split( line[:-1] )
            res1 = int( cols[ 0 ] ) - 1
            res2 = int( cols[ 1 ] ) - 1
            pair_map[ res1 ] = res2
            pair_map[ res2 ] = res1
            all_pairs.append( [res1,res2] )
            assert ( sequence[res1] in complement[ sequence[res2] ] )
        except:
            continue

    if line[:13] == 'CUTPOINT_OPEN':
        cutpoints_original = map( lambda x:int(x), string.split( line[14:] ) )

#print pair_map

def make_tag_with_dashes( int_vector ):
    tag = ''

    start_res = int_vector[0]
    for i in range( 1, len(int_vector)+1 ):
        if i==len( int_vector)  or  int_vector[i] != int_vector[i-1]+1:

            stop_res = int_vector[i-1]
            if stop_res > start_res:
                tag += ' %d-%d' % (start_res, stop_res )
            else:
                tag += ' %d' % (stop_res )

            if ( i < len( int_vector) ): start_res = int_vector[i]

    return tag

# Parse out stems
already_in_stem = {}
for i in range( numres ): already_in_stem[ i ] = 0

stems = []
stem_count = 0
for i in range( numres ):
    if pair_map.has_key( i ) and not already_in_stem[ i ]:  # In a base pair
        k = i
        stem_count += 1
        stem_res = []

        stem_res.append( [k, pair_map[k]] )
        already_in_stem[ k ] = 1
        already_in_stem[ pair_map[k] ] = 1

        # Can we extend in one direction?
        while( pair_map.has_key( k + 1 ) and  pair_map[ k+1 ] == pair_map[ k ] - 1  and  not already_in_stem[k+1] ):
            k += 1
            stem_res.append( [k, pair_map[k]] )
            already_in_stem[ k ] = 1
            already_in_stem[ pair_map[k] ] = 1

        # Do not allow single WC base pairs.
        if ( len( stem_res ) <2 ):
            print 'All stems must have length > 1 bp '
            print stem_res
            exit()
        stems.append( stem_res )

# Parse out motifs
already_in_motif = {}
for i in range( numres ): already_in_motif[ i ] = 0

motif_count = 0
motif_cutpoints = []
motif_res_maps = []
motif_stem_sets = []
motifs = []
for i in range( numres ):

    if not already_in_motif[i] and ( not already_in_stem[ i ] or
                                     ( i > 0 and already_in_stem[i-1] and \
                                           pair_map[i]+1 != pair_map[i-1] ) ):

        motif_count += 1
        motif_res = []
        motif_stem_set = []
        cutpoints = []

        if ( i > 1 ):

            # back up to beginning of stem.
            k = i - 1
            motif_stem = []

            # first base pair.
            motif_stem.append( [ k, pair_map[k] ] )
            motif_res.append( k )
            motif_res.append( pair_map[k] )

            k-= 1
            while k >= 0 and already_in_stem[k] and \
                    (pair_map[k]-1 == pair_map[k+1]):
                motif_stem.append( [ k, pair_map[k] ] )
                motif_res.append( k )
                motif_res.append( pair_map[k] )
                k -= 1
            motif_stem_set.append( motif_stem )
            k += 1
            cutpoints.append( pair_map[k] )

        #print 'AFTER FIRST HELIX: ', motif_res

        k = i
        while ( k not in motif_res and k < numres):
            # Move forward to next stem:
            while ( k < numres and not already_in_stem[ k ] ):
                if (already_in_motif[ k ] ):
                    print 'Hey cant deal with pseudoknots!'
                    exit()
                motif_res.append( k )
                k += 1

            stem_start = k

            if k >= numres:
                cutpoints.append( k-1 )
                break

            if k in motif_res : break

            motif_stem = []
            motif_stem.append( [ k, pair_map[k] ] )
            motif_res.append( k )
            motif_res.append( pair_map[ k ] )
            k+=1

            while ( k < numres and already_in_stem[ k ] and \
                        pair_map[k-1] == pair_map[k]+1 and not k in motif_res):
                motif_stem.append( [ k, pair_map[k] ] )
                motif_res.append( k )
                motif_res.append( pair_map[ k ] )
                k += 1
            motif_stem_set.append( motif_stem )
            cutpoints.append( k-1 )

            # Next non-helical part..
            k = pair_map[ stem_start ] + 1

            #print 'AFTER NEXT HELIX: ', motif_res

        motif_res.sort()

        motif_res_map = {}
        for k in range( len( motif_res ) ):
            motif_res_map[ motif_res[k] ] = k
            already_in_motif[ motif_res[k] ] = 1

        motifs.append( motif_res )
        motif_stem_sets.append( motif_stem_set )
        motif_cutpoints.append( cutpoints )
        motif_res_maps.append( motif_res_map )
        #print 'CUTPOINTS ', cutpoints

#print motifs
#print motif_stem_sets
#print motif_cutpoints

# Output stem definition jobs
readme_stems_file = 'README_STEMS'
fid_README_STEMS = open( readme_stems_file,'w')
#print 'INPUT_RES',
for i in range( stem_count ):

    # Fasta
    tag = 'stem%d_%s' % (i+1, fasta_file)
    fid = open( tag , 'w' )
    fid.write( '>'+tag+'\n')

    stem_res = stems[i]
    stem_length = len( stem_res )

    for k in range( stem_length ):
        fid.write( sequence[stem_res[k][0]] )
        #print stem_res[k][0]+1,
    for k in range( stem_length ):
        fid.write( sequence[stem_res[stem_length-k-1][1]] )
        #print stem_res[stem_length-k-1][1]+1,

    #print
    fid.write('\n')
    fid.close()
    print 'Created: ', tag

    # pdb_file
    if native_exists:
        command = 'python %s/pdbslice.py  %s -segments %d %d %d %d stem%d_' %( \
            tools_scripts_path,
            native_pdb_file,
            stem_res[0][0]+1,
            stem_res[-1][0]+1,
            stem_res[-1][-1]+1,
            stem_res[0][-1]+1,
            i+1 )
        print command
        system( command )

    outfile = 'stem%d_%s.out' % (i+1, fasta_file.replace('.fasta',''))
    stem_out_files.append( outfile )
    command = '%s/rna_helix%s  -fasta %s -out:file:silent %s' % (EXE_DIR, EXE_extension, tag, outfile)
    fid_README_STEMS.write(command+'\n')

#print

fid_README_STEMS.close()
print 'Created: ', readme_stems_file
print ' This has the command lines that you need to make helical stems'
print



def make_tag_with_dashes( int_vector ):
    tag = ''

    start_res = int_vector[0]
    for i in range( 1, len(int_vector)+1 ):
        if i==len( int_vector)  or  int_vector[i] != int_vector[i-1]+1:

            stop_res = int_vector[i-1]
            if stop_res > start_res:
                tag += ' %d-%d' % (start_res, stop_res )
            else:
                tag += ' %d' % (stop_res )

            if ( i < len( int_vector) ): start_res = int_vector[i]

    return tag

# Output motif jobs
short_commands = []
readme_motifs_file = 'README_MOTIFS'
fid_README_MOTIFS = open( readme_motifs_file,'w')
for i in range( motif_count ):

    # Fasta
    motif_fasta_file = 'motif%d_%s' % (i+1, fasta_file)
    fid = open( motif_fasta_file , 'w' )
    fid.write( '>'+motif_fasta_file+'\n')

    motif_res = motifs[i]
    motif_length = len( motif_res )

    for k in range( motif_length ):
        fid.write( sequence[motif_res[k]] )
    fid.write('\n')
    fid.close()
    print 'Created: ', motif_fasta_file

    # params file
    motif_stem_set = motif_stem_sets[ i ]
    motif_res_map = motif_res_maps[ i ]
    motif_cutpoint = motif_cutpoints[ i ]

    motif_params_file = 'motif%d_%s.params' % (i+1, fasta_file.replace('.fasta',''))
    fid = open( motif_params_file , 'w' )

    which_stems = []
    stem_chunk_res = []
    for k in range( len(motif_stem_set) ):
        motif_stem = motif_stem_set[ k ]
        fid.write( 'STEM   ' )
        for n in range( len( motif_stem ) ):
            fid.write( '  PAIR %d %d W W A'  % \
                       ( motif_res_map[ motif_stem[n][0] ]+1,
                         motif_res_map[ motif_stem[n][1] ]+1 ) )
        fid.write( '\n' )

        fid.write( 'OBLIGATE PAIR %d %d W W A \n\n'  % \
                       ( motif_res_map[ motif_stem[-1][0] ]+1,
                         motif_res_map[ motif_stem[-1][1] ]+1 ) )

        #need to find in stems
        for n in range( len( stems )  ):
            stem = stems[n]
            found_match = 0
            for q in range( len( stem ) ) :
                if motif_stem[0][0] in stem[q]:
                    found_match = 1
                    break
            if found_match: break
        which_stems.append( n )

        for q in range( len( stem ) ):
            stem_chunk_res.append( motif_res_map[ stem[q][0] ]+1 )
        for q in range( len( stem ) ):
            stem_chunk_res.append( motif_res_map[ stem[ len(stem) - q - 1][1] ]+1 )

    motif_cutpoint.sort()

    #print motif_res
    #print motif_cutpoint

    if ( len( motif_cutpoint ) > 1 ):
        fid.write( 'CUTPOINT_OPEN ' )
        for k in range( len( motif_cutpoint ) ):
            if motif_res_map[ motif_cutpoint[k] ] < (len( motif_res )-1):
                fid.write( ' %d' % (motif_res_map[ motif_cutpoint[k] ]+1) )
    fid.write('\n')
    fid.close()
    print 'Created: ', motif_params_file

    # pdb_file
    native_tag = ''
    if native_exists:
        command = 'python %s/pdbslice.py  %s -subset ' %  ( tools_scripts_path, native_pdb_file )
        for k in range( motif_length ): command += ' %d' % (motif_res[k]+1)
        command += ' motif%d_' % (i+1)
        print command
        system( command )
        native_pdb_file_subset =  'motif%d_%s' % (i+1, native_pdb_file )
        native_tag = '-native %s ' % native_pdb_file_subset
        print 'Created: ', native_pdb_file_subset

    if data_exists:
        motif_data_file = 'motif%d_%s' % ( i+1, data_file )
        fid_data = open( motif_data_file, 'w' )
        fid_data.write( 'EXPOSE' )
        for data in data_info:
            if data[0]-1 in motif_res_map.keys():
                fid_data.write( '   %d %s %s ' % (motif_res_map[data[0]-1]+1,data[1],data[2]) )
        fid_data.write('\n')

        if len( backbone_burial_info ) > 0:
            fid_data.write( 'BACKBONE_BURIAL ' )
            for k in backbone_burial_info:
                if k-1 in motif_res_map.keys():
                    fid_data.write( ' %d' % (motif_res_map[ k-1 ] + 1) )
            fid_data.write( '\n' )
        fid_data.close()
        print 'Created: ', motif_data_file

    cst_found = 0;

    if cst_exists:
        motif_cst_file = 'motif%d_%s' % ( i+1, cst_file )
        fid_cst = open( motif_cst_file, 'w' )
        fid_cst.write( '[ atompairs ]\n' )
        for cst in cst_info:
            if cst[1]-1 in motif_res_map.keys() and cst[3]-1 in motif_res_map.keys():
                fid_cst.write( '%s %d %s %d %s\n' % (cst[0], motif_res_map[cst[1]-1]+1,cst[2],motif_res_map[cst[3]-1]+1,cst[4]) )
                cst_found = 1
        fid_cst.close()
        print 'Created: ', motif_cst_file

    motif_out_file = motif_params_file.replace( '.params','.out')
    motif_out_files.append( motif_out_file )
    NSTRUCT = 100
    command = '%s/rna_denovo%s %s -fasta %s -params_file %s -nstruct %d -out:file:silent %s -cycles 5000 -mute all -close_loops -close_loops_after_each_move -minimize_rna' % \
        ( EXE_DIR, EXE_extension, native_tag, motif_fasta_file, motif_params_file, NSTRUCT, motif_out_file )

    if data_exists: command += ' -data_file %s ' % motif_data_file
    if cst_exists and cst_found: command += ' -cst_file %s ' % motif_cst_file
    if torsions_exists: command += ' -vall_torsions %s ' % torsions_file

    # new -- forcing exactly the same stems in all motifs.
    command += ' -in:file:silent_struct_type rna -in:file:silent '
    for n in which_stems: command += ' stem%d_%s.out' % (n+1, fasta_file.replace('.fasta',''))


    command += ' -input_res '
    command += ' '+make_tag_with_dashes( stem_chunk_res )

    short_command = command.replace('-nstruct 100', '-nstruct 1')
    short_command = short_command.replace('-cycles 5000', '-cycles 1')
    short_commands.append(short_command)

    fid_README_MOTIFS.write( command+'\n' )


fid_README_MOTIFS.close()


with open("README_MOTIFS.short", 'w') as shortfile:
	for short in short_commands:
		shortfile.write(short+'\n')


print 'Created: ', readme_motifs_file
print ' This has the command lines that you need to make interhelical motifs, like hairpin loops and internal junctions.'
print

# Output assembly job
#Where are the jumps and chainbreaks?
jumps = []

if len( cutpoints_original ) > 0:
    fid.write( 'CUTPOINT_OPEN ' )
    for cutpoint in cutpoints_original:
        fid.write( ' %d'  % (cutpoint+1) )
    fid.write( '\n' )


cutpoints  = []
for i in range( motif_count ):

    motif_stem_set = motif_stem_sets[ i ]

    motif_stem = motif_stem_set[ 0 ]

    possible_cutpoints =  [ motif_stem[ 0 ][ 0 ], motif_stem[ 1 ][ 1 ] ]
    possible_cutpoints.sort()
    #print possible_cutpoints
    if ( possible_cutpoints[0] not in cutpoints):
        cutpoints.append( possible_cutpoints[ 0 ] )


params_file = fasta_file.replace('.fasta','_assemble.params' )
fid = open( params_file, 'w')

if len( cutpoints ) > 0:
    fid.write( 'CUTPOINT_CLOSED ' )
    #cutpoints.sort()
    for cutpoint in cutpoints:
        fid.write( ' %d'  % (cutpoint+1) )
    fid.write( '\n' )

#for cutpoint in cutpoints:
#    fid.write( 'OBLIGATE   PAIR %d %d W W A\n' % (cutpoint+1, pair_map[cutpoint]+1) )

for i in range( stem_count ):
    stem_res = stems[i]
    fid.write( 'STEM ')
    for k in range( len( stem_res )):
        fid.write( ' PAIR %d %d W W A ' % \
                       ( stem_res[k][ 0 ]+1, stem_res[k][ 1 ]+1 ) )
    fid.write('\n')
    fid.write( 'OBLIGATE PAIR %d %d W W A \n\n'  % \
                   ( stem_res[-1][0] + 1,\
                     stem_res[-1][1] + 1 ) )

fid.close()


print 'Created: ', params_file

########
assemble_cst_file = params_file.replace('.params','.cst')
if cst_exists:
    assemble_cst_file = cst_file.replace('.cst','_assemble.cst')
fid = open( assemble_cst_file,'w')
fid.write('[ atompairs ]\n')
for cutpoint in cutpoints:
    fid.write( "O3'  %d  P     %d  HARMONIC  1.619  2.0\n" % \
                   ( cutpoint+1, cutpoint+2 ) )
if cst_exists:
    for cst in cst_info:
            fid.write( '%s %d %s %d %s \n' % (cst[0], cst[1], cst[2], cst[3], cst[4]) )

fid.close()
print 'Created: ', assemble_cst_file


#########
readme_assemble_file = 'README_ASSEMBLE'
fid = open( readme_assemble_file, 'w' )

native_tag = ''
if native_exists: native_tag = '-native '+native_pdb_file

outfile = params_file.replace( '.params','.out' )
command = '%s/rna_denovo%s %s -fasta %s -in:file:silent_struct_type binary_rna  -cycles 10000 -nstruct 200 -out:file:silent %s -params_file %s -cst_file %s -close_loops  -in:file:silent ' % \
( EXE_DIR, EXE_extension, native_tag, fasta_file, outfile, params_file, assemble_cst_file )

for stem_out_file in stem_out_files:
    command += ' '+stem_out_file
for motif_out_file in motif_out_files:
    command += ' '+motif_out_file

chunk_res = []
command += ' -input_res '
for n in range( len( stems )  ):
    stem = stems[n]
    for q in range( len( stem ) ):        chunk_res.append(stem[q][0] + 1)
    for q in range( len( stem ) ):        chunk_res.append(stem[ len(stem) - q - 1][1] + 1)

for n in range( motif_count ):
    motif_res = motifs[n]
    for m in motif_res: chunk_res.append(m+1)

command += make_tag_with_dashes( chunk_res )


if torsions_exists: command += ' -vall_torsions %s ' % torsions_file
if data_exists:
    command += ' -data_file '+data_file

fid.write( command+'\n')
fid.close()

with open('README_ASSEMBLE.short','w') as assemble_shortfile:
	command = command.replace('-cycles 10000', '-cycles 1')
	command = command.replace('-nstruct 200', '-nstruct 1')
	assemble_shortfile.write( command )

print 'Created: ', readme_assemble_file
print ' This has the command lines that you need to make full-length structures after making the individual stems and interhelical motifs'
print

os.system('chmod a+x README_*')
