#!/usr/bin/env python
# python scorefile_to_files scorefile #
# the # represents the number of sets of _0001_0001 to ignore at the end for sorting
import sys
import string
import os

def description_from_line( line ) :
    return line.split(" ")[ -1 ]

def batch_name_from_description( desc, n_underscores ) :
    return desc.rsplit('_', n_underscores)[0]

if __name__ == '__main__':
    fname = sys.argv[ 1 ]
    n_underscores = 1
    if len(sys.argv) >= 3:
        n_underscores = int( sys.argv[ 2 ] )
    flines = open( fname ).readlines()
    score_lines_for_batches = {}
    first = True
    if flines[0].rstrip() == 'SEQUENCE:':
        del flines[0]
    for line in flines :
        if first :
            first = False
            continue
        batchname = batch_name_from_description( description_from_line( line ), n_underscores )
        if not batchname in score_lines_for_batches :
            score_lines_for_batches[ batchname ] = []
            score_lines_for_batches[ batchname ].append( flines[ 0 ] )
        score_lines_for_batches[ batchname ].append( line )

    for batch in score_lines_for_batches :
        lines = score_lines_for_batches[ batch ]
        fout_name = "score." + batch + ".sc"
        open( fout_name, "w" ).writelines( lines )
