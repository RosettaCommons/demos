#!/usr/bin/env python

from __future__ import print_function, division
import argparse,sys

def read_options():
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     description="Extract columns from Rosetta scorefile",
                                     usage="./extract_scores.py <scorefile> --columns score rms nmr_pcs [nmr_rdc ...]")
    parser.add_argument('--columns', type=str, nargs='+', dest='columns', required=True,
                        help='Vector of column labels to be extracted.')
    parser.add_argument('scorefile', type=argparse.FileType('r'), metavar='scorefile',
                        help='Name of the Rosetta scorefile.')
    args = parser.parse_args()

    return args.columns, args.scorefile

def read_scores(scorefile,columnlabels):
    lines=scorefile.readlines()
    nline=0
    foundheader=False
    colidxs=[]
    colla=[]

    while not foundheader and nline<len(lines):
        if "description" in lines[nline]:
            foundheader=True
            line=lines[nline].strip().split()
            ncols=len(line)
            for la in columnlabels:
                try:
                    idx=line.index(la)
                except ValueError as e:
                    print("Warning: Column %s could not be found in scorefile header and will be skipped."%la)
                else:
                    colidxs.append(idx)
                    colla.append(la)
        nline+=1

    if not foundheader:
        print("Error: Could not find scorefile header.")
        sys.exit(1)

    if len(colidxs)==0:
        print("Error: Zero column labels found in scorefile header.")
        sys.exit(2)

    nline=0
    scores=[[] for i in range(len(colidxs))]
    for line in lines:
        nline+=1
        line=line.strip().split()
        if len(line)!=ncols:
            print("Warning: Number of fields in line does not match score file header.")
            print("Warning: Line %d will be skipped."%nline)
            continue
        if "description" in line:
            continue
        for i, idx in enumerate(colidxs):
#            scores[i].append(float(line[idx]))
            scores[i].append(line[idx])

    return colla, scores

def write_scores(outfile,scores,columnlabels):
    with open(outfile,'w') as f:
        f.write("%s\n"%(','.join(columnlabels)))
        for row in map(list, zip(*scores)):
            f.write("%s\n"%(','.join([str(s) for s in row])))

def main():
    # Read input options
    colla,scfile=read_options()

    # Read scores from scorefile
    print("Reading scores from file %s ...\n"%scfile.name)
    collanew,scores=read_scores(scfile,colla)
    print("Finished reading scores for: %s\n"%(', '.join(collanew)))

    # Write extracted scores to comma-separated file
    print("Writing scores ...\n")
    write_scores("scores.csv",scores,collanew)
    print(" * * * DONE * * * ")

main()
