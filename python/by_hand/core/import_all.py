# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
import sys
sys.path.append("/Users/sheffler/svn/branches/mini/python/by_hand/core/lib")

libs = ("_containers","_chemical","_conformation","_graph","_kinematics","_mm",
        "_optimization","_options","_pose","_scoring","_util")

tplt = "print '%(lib)s'\nimport %(lib)s\nprint dir(%(lib)s)\nprint"
for lib in libs:
  exec tplt%vars()

