#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#
from sys import stderr

def extract_fragid( fragfile ):
    assert fragfile.startswith("after"), fragfile
    assert isinstance( fragfile, str )

    l = fragfile.split(".")
    mer         = int(float( l[1] ))
    pos         = int(float( l[2] ))
    picker_rank = int(float( l[3] ))
    SHD_rank    = int(float( l[4] ))

    return ( mer, pos, picker_rank, SHD_rank )


def fragid_to_fragfile( fragid ):
    assert isFragID( fragid )

    if isNullFrag( fragid ):
        frag_fn = "null"
    else:
        frag_fn = "after_rotation_frags.%s.%s.%s.%s.????.pdb" % fragid()

    return frag_fn



################################### all those assertions ###################################
def isPose( pose ):
    #if pose.__class__.__name__ == "Pose" and pose._initialized: # this is probably not necessary
    if pose.__class__.__name__ == "Pose": # this is probably not necessary
        return True
    else:
        return False


def isFragIdxPose( pose ):
    #if pose.__class__.__name__ == "FragIdxPose" and pose._initialized: # this is probably not necessary
    if pose.__class__.__name__ == "FragIdxPose": # this is probably not necessary
        return True
    else:
        return False


def isFragIdx( frag_idx ):
    ''' what is the best way to determine a given input is frag_idx or not -
        maybe ask ScoreTable to do that?  '''
    if isinstance( frag_idx, int ):
        return True
    else:
        return False


def isFragID( frag_id ):
    if frag_id.__class__.__name__ == "FragID":
        return True
    else:
        return False


def isFragTuple( frag_tuple ):
    if isinstance( frag_tuple, tuple ) and len( frag_tuple )==4:
        return True
    else:
        return False


def isScoreTable( scoretable ):
    if scoretable.__class__.__name__ == "ScoreTable":
        return True
    else:
        return False


def isScoreFxn( scorefxn ):
    if scorefxn.__class__.__name__ == "ScoreFunction":
        return True
    else:
        return False


def isResidue( residue ):
    if residue.__class__.__name__ == "Residue":
        return True
    else:
        return False


def isWeights( weights ):
    if weights.__class__.__name__ == "Weights" :
        return True
    else:
        return False


def isBoltzmann( boltzmann ):
    if boltzmann.__class__.__name__ == "Boltzmann":
        return True
    else:
        return False


def isTracker( tracker ):
    if tracker.__class__.__name__ == "TrajectoryTracker":
        return True
    else:
        return False


def isNullFrag( frag_id ):
    assert isFragID( frag_id ) or isFragIdx( frag_id ) or isFragTuple( frag_id )
    if isFragID( frag_id ):
        if ( frag_id._picker_rank, frag_id._SHD_rank ) == ( 0, 0 ):
            return True
        else:
            return False

    elif isFragIdx( frag_id ):
        if frag_id < 0:
            return True
        else:
            return False

    elif isFragTuple( frag_id ):
        if ( frag_id[2], frag_id[3] ) == ( 0, 0 ):
            return True
        else:
            return False
    else:
        #print "isNullFragThis shouldn't happen", frag_id
        stderr.write("isNullFrag.ERROR: this function takes either FragID(tuple) or FragIdx: %s\n" % frag_id )
        return False

