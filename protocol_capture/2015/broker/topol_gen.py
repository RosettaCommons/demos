import sys
import argparse
import top_file

def process_command_line(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("--strands", required=True, type=str, metavar="(xx-yy),(zz-aa)", 
                        help="The regions (one indexed) that are strands.")
    parser.add_argument("--topology", required=True, type=str, metavar="(w-x:P),(y-z:A)",
                        help="The ways in which the given strands are fit together." )
    parser.add_argument("--len_shift", default=0, type=int, metavar="DELTA",
                        help="Allow strands to shift up to DELTA against each other." )

    args = parser.parse_args(argv[1:])

    strand_defs = []
    for pair in args.strands.split(','):
        (i_str, j_str) = pair.split('-')
        (i_str, j_str) = ( i_str.lstrip(), j_str.rstrip() )
        (i, j) = ( int( i_str[1:] ), int( j_str[0:-1] ) )
        
        strand_defs.append( ( i, j ) )

    strand_pairs = []
    for triplet in [s.rstrip().lstrip()[1:-1] for s in args.topology.split(",")]:
        (p1, substr) = triplet.split("-")
        (p2, para) = substr.split(":")
        (p1, p2) = (int(p1), int(p2))
        assert( para == "P" or para == "A" )
        strand_pairs.append( (p1, p2, para) )

    args.strands = strand_defs
    args.topology = strand_pairs

    if( args.len_shift > 0 ):
        reduced = set()
        for ( s1, s2, para ) in strand_pairs:
            for s in [ s1, s2 ]:
                lindiv = all( [ s1 not in reduced and 
                                s2 not in reduced for ( s1, s2, para ) in strand_pairs
                                if s == s1 or s == s2] )
                if( lindiv ):
                    args.strands[s] = ( args.strands[s][0], 
                                         args.strands[s][1] - args.len_shift )
                    reduced.add( s )
    elif( args.len_shift < 0 ):
        sys.stderr.write( "--len_shift flag values must be >= 0." )
        sys.exit(1)

    return args

def main(argv=None):
    
    args = process_command_line(argv)

    tfile = top_file.TopologyFile( args.strands,
                                   [args.topology] )

    print tfile

    return 1

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)
