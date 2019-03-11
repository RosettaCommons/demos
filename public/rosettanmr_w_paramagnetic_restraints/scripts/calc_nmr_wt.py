#!/usr/bin/env python

from __future__ import print_function, division
import argparse,sys

def read_options():
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     description="Calculate NMR weights from scores in Rosetta scorefile",
                                     usage="./calc_nmr_wt.py <scorefile> --nmr_scoretypes pcs [ rdc pre ]")
    parser.add_argument('scorefile', type=argparse.FileType('r'), metavar='scorefile',
                        help='Name of the Rosetta scorefile.')
    parser.add_argument('--nmr_scoretypes', type=str, dest='scoretypes', nargs='+',
                        required=True, choices=['pcs','rdc','pre'], help='List of NMR scoretypes.')
    args = parser.parse_args()

    return args.scorefile, args.scoretypes

def read_scores(scorefile,scoretypes):
    lines=scorefile.readlines()
    nline=0
    foundheader=False
    colidxs={}
    scoretypes.append('total_score')
    scoretypes.append('score')
    while not foundheader and nline<len(lines):
        if "description" in lines[nline]:
            foundheader=True
            line=lines[nline].strip().split()
            ncols=len(line)
            for st in scoretypes:
                try:
                    idx=line.index(st)
                except ValueError as e:
                    print("Warning: NMR scoretype %s was not found in scorefile header and will be skipped."%st)
                else:
                    colidxs[st]=idx
        nline+=1

    if not foundheader:
        print("Error: Could not find scorefile header.")
        sys.exit(1)

    if len(colidxs)==0:
        print("Error: No NMR scoretypes found in scorefile header.")
        sys.exit(2)

    if 'total_score' not in colidxs.keys() and 'score' not in colidxs.keys():
        print("Error: No total_score or score column found in scorefile.")
        sys.exit(3)
    else:
        if 'total_score' in colidxs.keys() and 'score' in colidxs.keys():
            del colidxs['score']
        elif 'total_score' not in colidxs.keys() and 'score' in colidxs.keys():
            colidxs['total_score']=colidxs.pop('score')
        else: # Nothing to be done here
            pass

    nline=0
    scores=dict([(k,[]) for k in colidxs.keys()])
    for line in lines:
        nline+=1
        line=line.strip().split()
        if len(line)!=ncols:
            print("Warning: Number of fields in line does not match score file header.")
            print("Warning: Line %d will be skipped."%nline)
            continue
        if "description" in line:
            continue
        for st, idx in colidxs.items():
            scores[st].append(float(line[idx]))

    return scores

def main():
    # Read input options
    scfile,sctypes=read_options()

    # Read score values
    print("Reading scores from file %s ...\n"%scfile.name)
    rosetta_sctype_names={'pcs':'nmr_pcs','rdc':'nmr_rdc','pre':'nmr_pre'}
    scores=read_scores(scfile,[rosetta_sctype_names[st] for st in sctypes])
    
    # Calculate score weights
    print("Calculating NMR score weights ...\n")

    # First, subtract NMR scores from total score
    for st, vals in scores.items():
        if st=='total_score':
            continue
        assert len(vals)==len(scores['total_score']),\
               "Number of total_score values does not match "\
               "number of % score values. Check if scorefileis correct."%st
        for i,v in enumerate(vals):
            scores['total_score'][i]-=v

    # Second, sort scores and calculate averages
    # of lower and upper 10 percent
    avgs={}
    tenpercent=int(len(scores['total_score'])/10)
    sclimit=1.0e5 # score limit to avoid large outliers
    for st, vals in scores.items():
        vals.sort()
        assert tenpercent<len(vals)

        avglo=0
        avghi=0
        i=0
        for v in vals[:tenpercent]:
            if v>sclimit:
                continue
            avglo+=v
            i+=1
        avglo/=float(i)

        i=0
        for v in vals[-tenpercent:]:
            if v>sclimit:
                continue
            avghi+=v
            i+=1
        avghi/=float(i)

        avgs[st]=(avglo,avghi)

    # Third, calculate NMR weights
    ros_sclo=avgs['total_score'][0]
    ros_schi=avgs['total_score'][1]
    for st, avg in avgs.items():
        if st=='total_score':
            continue
        else:
            nmr_sclo=avg[0]
            nmr_schi=avg[1]
            wt=(ros_schi-ros_sclo)/(nmr_schi-nmr_sclo)
            print("%s = %.5f\n"%(st,wt))

    print(" * * * DONE * * * ")

main() 
