#!/bin/sh
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#this script is for "deflating" large files in your demo directory.  It will zip up your big file, create a symlink from the old filename to the gzipped version, and make SVN happy with the changes.

gzip $1
svn rm $1
ln -s $1.gz $1
svn add $1
svn add $1.gz

