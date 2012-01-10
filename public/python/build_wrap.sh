# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
libRos="/Users/sheffler/svn/branches/mini"

CC=g++ #-dp-4.2
#OPTS="  -O3 -c -ftemplate-depth-256  -fno-inline -fPIC"   #" -Wall"
OPTS="-c -pipe -ffor-scope -O0 -g -ggdb -ffloat-store"
INCL="-I. -I/opt/local/include  -I/usr/include -I$libRos/src -I$libRos/src/platform/macos -I/Library/Frameworks/Python.framework/Versions/2.4/include/python2.4 -I/opt/local/boost"
FRMWK="" #"-F/System/Library/Frameworks"
BSTSRC="/opt/local/boost/libs/python/src"
BOOSTFLAGS="-DBOOST_PYTHON_DYNAMIC_LIB -DBOOST_PYTHON_SOURCE"
LFLAGS="-L.  -L$libRos/build/src/debug/macos/10.4/32/x86/gcc"

src=relax

cmd1="$CC $OPTS -DBOOST_PYTHON_DYNAMIC_LIB $FRMWK $INCL -o build/wrap_$src.o $src.cc"
echo $cmd1
$cmd1

cmd2="$CC -Wl,-stack_size,4000000,-stack_addr,0xc0000000 -bundle -o _test.so build/wrap_$src.o $FRMWK -framework Python $LFLAGS -lcore -lprotocols -lnumeric -lObjexxFCL -lutility -lboost_python -lz"
echo $cmd2
$cmd2

