#!/usr/bin/env python2.7
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#
from heapq import *
import argparse, os, os.path, re, shutil
from sys import exit, stderr, stdout, path
from os import popen, system
from os.path import basename, exists
from multiprocessing import Pool
import pprint
from glob import glob

import denovo_model_building_util


def extract_silent( args, name_map ):
    '''Writes the low-scoring entries from the input silent file to the output silent file'''
    tags = name_map.keys()
    hash = set(tags)

    with open(args.silent_extract, 'w') as fout:
        with open(args.silent_in) as fin:
            fout.write(fin.readline())  # SEQUENCE
            fout.write(fin.readline())  # SCORE header

            for line in fin:
                cols = line.split()
                tag = cols[-1].strip()

                if tag in hash:
                    fout.write(line)


def run_extract_systemcall( silentfile, name_map="" ):
    config = denovo_model_building_util.read_config_file( args.config_file )
    PATH = config.get("path", "demo_dir")
    bin = PATH + "/rosetta/extract_pdbs.static.linuxgccrelease"
    db  = PATH + "/rosetta/rosetta_trunk_database"
    assert exists( bin ), 'Failed to locate extract_pdbs binary, not in %s' %bin

    extract = [
                bin,
                ' -in:file:silent ' + silentfile,
                ' -in:file:silent_struct_type ' + args.silent_struct_type,
                ' -database ' + db,
                ' -chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation SpecialRotamer protein_cutpoint_upper protein_cutpoint_lower VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm ',
                ' -mute all '
               ]

    if args.fa_param:
        param_file = ' -extra_res_fa %s ' % args.fa_param
        extract.append( param_file )

    if args.cen_param:
        param_file = ' -extra_res_cen %s ' % args.cen_param
        extract.append( param_file )


    if args.xtal_refine:
        xtal_refine = ' -crystal_refine '
        extract.append( xtal_refine )

    if args.symm_def:
        symm_def = ' -symmetry_definition %s ' % args.symm_def
        print symm_def
        extract.append( symm_def )

    if name_map:
        tags = ' -tags ' + ' '.join( name_map.keys() )
        extract.append( tags )

    cmd = " ".join( extract )
    print cmd

    os.system( cmd )



def extract_dir( args, name_map ):
    '''Writes the low-scoring entries from the input silent file to PDBs in the output directory'''

    stdout.write("Extracting %s... \n" % args.silent_extract)

    # Execute the command, then move extracted files to the output directory
    run_extract_systemcall( args.silent_extract, name_map )


    # Adding tags and move to output_dir if specified
    for tag in name_map.keys():
        src = tag + '.pdb'
        if args.label:
            link = name_map[tag][1]+'.pdb'
            if exists( link ): os.unlink( link )
            os.symlink(src, link)

        if args.output_dir:
            dst = args.output_dir + '/' + src
            if exists(src):
                if args.label:
                    #print link, dst
                    #os.system("mv %s %s"%(link, args.output_dir))
                    shutil.move(link, args.output_dir)
                shutil.move(src, dst)

    if args.output_dir:
        shutil.move(args.silent_extract, args.output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config_file', default="denovo_model_building_scripts.cfg", help='')

    # Input args
    parser.add_argument('silent_in', help='Input silent file')
    parser.add_argument('-ist', '--silent_struct_type', type=str, default="binary", choices=["protein", "binary"], help='Input silent struct type ')

    # Output args
    parser.add_argument('-so', '--silent_only', action="store_true", help='Output only silent file without extractint to pdbs')
    parser.add_argument('-o', '--silent_extract', default='tmp_SELECT_silent.out', help='Output silent file')
    parser.add_argument('-d', '--output_dir', help='Directory to store extracted models')
    parser.add_argument('-ow', '--overwrite', action="store_true", help='Delete the output if it exists')
    parser.add_argument('-n', '--num_models', type=int, default=10, help='Number of models to extract')
    parser.add_argument('-p', '--percentage', type=float, help='percentage of models to extract - this would overwrite num_models option')
    parser.add_argument('-l', '--label', action="store_true", help='label the output with the ranking by making symlink')
    parser.add_argument('--fa_param', help="ligand param")
    parser.add_argument('--cen_param', help="ligand param")

    parser.add_argument('--muteall', action="store_true", default=False, help='')

    # Filtering args
    parser.add_argument('--taglist', nargs="+", help='If you want to get the best GDTMM decoys, use this option!')
    parser.add_argument('--taglist_col', default=1, type=int, help='If you want to get the best GDTMM decoys, use this option!')
    parser.add_argument('-r', '--reverse_sort', action="store_true", help='If you want to get the best GDTMM decoys, use this option!')
    parser.add_argument('-c', '--column_name', default='score', help='Energy terms to sort in heap.')
    parser.add_argument('--exclude_terms', default='', nargs='+', help='Energy terms to subtract from score. Supports regular expressions.')
    parser.add_argument('--symm_def', help='symm_def file')
    parser.add_argument('-x', '--xtal_refine', action="store_true", help='allow chain break')

    args = parser.parse_args()

    #assert (args.silent_extract or args.output_dir), 'Must provide either --silent_extract or --output_dir'
    #assert args.silent_extract, 'Must provide either --silent_extract or --output_dir'
    assert exists( args.silent_in ), 'Failed to open input silent file'
    assert args.num_models > 0


    # Removes the output directory if it already exists
    if args.output_dir:
        if exists(args.output_dir):
            if not args.overwrite:
                stderr.write("ERROR: folder %s exists. Did you forget pass -ow/--overwrite ?\n" % args.output_dir); exit()
            else:
                shutil.rmtree(args.output_dir)
        os.mkdir(args.output_dir)

    if not args.taglist:
        with open(args.silent_in) as file:
            file.readline()

            # Columns containing score, tag
            header = file.readline().split()
            try:
                idx_score = header.index(args.column_name)
            except ValueError:
                #pprint.pprint( header )
                stdout.write("%s\n" % "\n".join(header))
                stderr.write("ERROR: couldn't find the score term [ %s ].\n" % args.column_name);exit()

            idx_descr = header.index('description')

            # Columns containing score terms to exclude
            excluded_idxs = []
            for e in args.exclude_terms:
                pattern = re.compile(e)

                for (i, term) in enumerate(header):
                    if pattern.search(term):
                        print 'Excluding %s, index %d' % (term, i)
                        excluded_idxs.append(i)

            heap = []

            n_total_models = int( popen("grep ^SCORE: %s | wc -l" % args.silent_in ).readline().strip() ) - 1
            if n_total_models > args.num_models:
                count = args.num_models
            else:
                count = n_total_models

            if args.percentage:
                count = int(n_total_models*args.percentage/100)

            for line in file:
                line = line.strip()
                if not line.startswith('SCORE:'):
                    continue

                cols = line.split()
                score = float(cols[idx_score])
                try:
                    descr = cols[idx_descr]
                except IndexError:
                    stderr.write("WARNING: read through 1 bad decoy %s \n" % line.split()[-1] )

                # Subtract all excluded terms from the score
                additional = [float(cols[i]) for i in excluded_idxs]
                score -= sum(additional)

                # If want lowest-score: nagate score; vice versa - because we're inserting into a min heap
                reverse = 1 if args.reverse_sort else -1
                #print "reverse:", reverse
                entry = ( reverse*score, descr )

                if len(heap) < count:
                    heappush(heap, entry)
                else:
                    heappushpop(heap, entry)

        #tags = []

        sort_ = "highest" if reverse == 1 else "lowest"
        if not args.muteall:
            if not args.silent_only:
                stdout.write("You are going to extract %s %s %s PDBs...\n" %( count, sort_, args.column_name ))
            else:
                stdout.write("You are going to extract %s %s %s as a silent file...\n" %( count, sort_, args.column_name ))


        name_map = {}
        while heap:
            score, descr = heappop(heap)
            #tags.append(descr)
            score = score*-1
            link_name = "%s_%s_model_%s" %( sort_, args.column_name, count )
            name_map[descr] = ( score, link_name )
            if not args.muteall:
                stdout.write( "%.3f \t %s \t %s_%s_model_%s.pdb\n" %( score, descr, sort_, args.column_name, count ))
            count -= 1

        #args.output_dir =  "%s_%s_models" %( sort_, args.column_name )
        #args.silent_extract =  "%s_%s_models_SELECT_silent.out" %( sort_, args.column_name )

        # Write the low-score decoys to the appropriate sink
    else:
        name_map = {}
        for taglist in args.taglist:
            if exists( taglist ):
                with open( taglist, "r" ) as f:
                    for l in f:
                        if not l.strip(): continue
                        if l.startswith("#"): continue
                        tag = l.split()[args.taglist_col-1]
                        name_map[ tag ] = 100000
            else:
                name_map[ taglist ] = 100000

        print name_map.keys()

    if exists( args.silent_extract ):
        if not args.overwrite:
            stderr.write("ERROR: silent_extract %s exists. Did you forget pass -ow/--overwrite ?\n" % args.silent_extract); exit()
        else:
            os.remove( args.silent_extract )

    extract_silent( args, name_map ) #if args.silent_extract else extract_dir

    if not args.silent_only:
        extract_dir( args, name_map )

