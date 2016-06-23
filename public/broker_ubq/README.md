Topology Broker - AbInitio with a given beta strand topology
=========================

KEYWORDS: STRUCTURE_PREDICTION GENERAL

This protocol runs abinitio with a given beta strand topology and a helix fixed 
at atomic resolution using the Topology Broker.

Assuming Rosetta is in your home directory, you run it as follows. 
(For production runs, 10000 structures should be more than enough to see what's going on, but will take a VERY long time to run without a cluster, so we'll reduce it to 10 for an example):

    $ rosetta_scripts.default.linuxgccrelease @flags -nstruct 10

where platform is "linux" or "mac" and compiler is "gcc" or "clang". 

