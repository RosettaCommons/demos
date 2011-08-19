# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
libRos="/Users/sheffler/svn/branches/mini"

CC=g++ #-dp-4.2
#OPTS="  -O3 -c -ftemplate-depth-256  -fno-inline -fPIC"   #" -Wall"
OPTS="-c -pipe -ffor-scope -Wno-long-double -O0 -g -ggdb -ffloat-store"
INCL="-I. -I/opt/local/include  -I/usr/include -I$libRos/src -I$libRos/src/platform/macos -I/Library/Frameworks/Python.framework/Versions/2.4/include/python2.4 -I/opt/local/boost"
FRMWK="" #"-F/System/Library/Frameworks"
BSTSRC="/opt/local/boost/libs/python/src"
BOOSTFLAGS="-DBOOST_PYTHON_DYNAMIC_LIB -DBOOST_PYTHON_SOURCE"
LFLAGS="-L.  -L$libRos/build/src/debug/macos/10.4/32/x86/gcc "
lFLAGS="-lcore -lnumeric -lObjexxFCL -lutility -lboost_python -lz -lprotocols"

mkdir -p lib
for i in test_relax; do #containers chemical conformation graph kinematics mm optimization pose scoring util; do # no graph pack io
  if [ ! -e lib/_$i.so ]; then 
    cmd1="$CC $OPTS -DBOOST_PYTHON_DYNAMIC_LIB $FRMWK $INCL -o lib/wrap_$i.o wrap/wrap_$i.cc"
    echo $cmd1
    $cmd1 &> log/build_$i.log
    if [ -e lib/wrap_$i.o ]; then 
      cmd2="$CC -bundle -o lib/_$i.so lib/wrap_$i.o $FRMWK -framework Python $LFLAGS $lFLAGS"
      echo $cmd2
      $cmd2 &> log/build_$i.log      
    fi
  fi
  if [ ! -s log/build_$i.log ]; then 
    rm -f log/build_$i.log ; 
  fi
done


