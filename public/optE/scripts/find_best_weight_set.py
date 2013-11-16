import sys

# Python script to read through the log file from an optE job and
# determine which weight set was the best of those examined.
# The first (and only) argument should be the name of the log file.


def check_int(s):
    if s[0] in ('-', '+'):
        return s[1:].isdigit()
    return s.isdigit()

def find_best_weight_set_from_optElog( logfilename ) :
    lines = open( logfilename ).readlines()
    inner=0
    outer=0
    lastaccept_inner=0
    lastaccept_outer=0
    lastaccept_seqrec=0
    lastaccept_crossent=0
    gomaster_string = "protocols.optimize_weights.IterativeOptEDriver: (0) go(): master node"
    accept_string = "protocols.optimize_weights.IterativeOptEDriver: (0) decide_if_sequence_recovery_improved(): accepting "
    for line in lines :
        gostart = line.find( gomaster_string )
        if gostart != -1 :
            remainder = line[ gostart + len(gomaster_string) : ]
            #print "remainder:", remainder
            outerinnerstring = remainder.partition("(")[2].partition(")")[0]
            #print "outerinnerstring", outerinnerstring
            if outerinnerstring :
                parts = outerinnerstring.partition(",")
                if check_int(parts[0]) and check_int( parts[2]) :
                    outer = int(parts[0])
                    inner = int(parts[2])
            continue
        acceptstart = line.find( accept_string )
        if acceptstart != -1 :
            lastaccept_inner = inner
            lastaccept_outer = outer
            remainder = line[ acceptstart + len( accept_string ):]
            # new weight set: inner loop recovery rate: 0.374774, outer loop recovery rate: 0.371293, inner loop weighted cross entropy: -0.288942, outer loop weighted cross entropy: -0.289134
            cols = remainder.split()
            lastaccept_seqrec = float(cols[7][:-1])
            if len(cols) > 18 :
                lastaccept_crossent = float(cols[18][:-1])*-10
            else :
                lastaccept_crossent = 0
    return lastaccept_outer, lastaccept_inner, lastaccept_seqrec, lastaccept_crossent

if __name__ == "__main__" :
    outer, inner, seqrec, crossent = find_best_weight_set_from_optElog( sys.argv[1] )
    print "Round that produced the last accepted weight set:", outer, inner, seqrec, crossent
