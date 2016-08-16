#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#
from Util import *

class FragID:
    ''' try to make it immutable by using __getattr__ to intercept changing '''
    ''' the original implemented tuple is safer actually '''
    def __init__( self, frag_id ):
        assert len( frag_id )==4
        assert isinstance( frag_id, tuple )

        self._mer         = frag_id[0]
        self._pos         = frag_id[1]
        self._picker_rank = frag_id[2]
        self._SHD_rank    = frag_id[3]


    def __call__( self ):
        return ( self._mer, self._pos, self._picker_rank, self._SHD_rank )


class Residue:
    def __init__( self, frag_id, score=0.0, density_score=0.0, overlap_score=0.0, closab_score=0.0, clash_score=0.0, rmsd=0.0, boltzmann=None ): # missing Boltzmann
        ''' '''
        assert isFragID( frag_id )

        self._frag_id       = frag_id
        self._score         = score
        self._density_score = density_score
        self._overlap_score = overlap_score
        self._closab_score  = closab_score
        self._clash_score   = clash_score
        self._rmsd          = rmsd
        self._boltzmann     = boltzmann


    def update_boltzmann( self, boltzmann ):
        ''' this is being used right after query a selected_frag_id from scorefxn_object '''
        assert isBoltzmann( boltzmann )
        self._boltzmann = boltzmann


