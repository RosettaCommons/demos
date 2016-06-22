AbInitio Protocol Capture
=========================

This protocol runs abinitio with a given beta strand topology and a helix fixed 
at atomic resolution.

Assuming Rosetta is in your home directory, you run it as follows(10000 
structures should be more than enough to see what's going on, but will take a VERY long time to run without a cluster):

    $> rosetta_scripts.default.linuxgccrelease @flags -nstruct 10

where platform is "linux" or "mac" and compiler is "gcc" or "clang". 
