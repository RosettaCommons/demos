# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
mini=/Users/sheffler/svn/branches/mini/

g++ -o build/relax.o -c -pipe -ffor-scope -W -Wall -pedantic -fno-exceptions -Wno-long-double -O0 -g -ggdb -ffloat-store -I$mini/src -I$mini/external/include -I$mini/src/platform/macos/32/gcc -I$mini/src/platform/macos/32 -I$mini/src/platform/macos -I/usr/local/include -I/usr/include $mini/demo/python/relax.cc

g++ -o python_relax -Wl,-stack_size,4000000,-stack_addr,0xc0000000 build/relax.o -L$mini/lib -L$mini/external/lib -L$mini/build/src/debug/macos/10.4/32/x86/gcc -L/usr/local/lib -L/usr/lib -lnumeric -lutility -lObjexxFCL -lcore -lprotocols -lz -lnumeric -lutility -lObjexxFCL -lcore -lprotocols -lz
