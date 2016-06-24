Running the fixbb protocol
==========================

KEYWORDS: STRUCTURE_PREDICTION GENERAL

(c) Copyright Rosetta Commons Member Institutions.
(c) This file is part of the Rosetta software suite and is made available under license.
(c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
(c) For more information, see http://www.rosettacommons.org. Questions about this can be
(c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

Design in Centroid Mode:

../../bin/fixbb.macosgccrelease -database ~/minirosetta_database/ -l list.txt -ignore_unrecognized_res -centroid_input -score:weights score3 -mute core.io core.conformation -nstruct 2

This will produce 2 output sequences for each of the pdb files listed in list.txt.

Design in Fullatom Mode:

../../bin/fixbb.macosgccrelease -database ~/minirosetta_database/ -l list.txt -ignore_unrecognized_res  -mute core.io core.conformation -nstruct 2

By default, uses standard weights with score12 patch.

Fixbb is not smart enough, yet, to reuse precomputed pair energies for multiple
designs on the same backbone.

