AbInitio Protocol Capture
=========================
KEYWORDS: DESIGN GENERAL
This protocol runs abinitio with a given beta strand topology and a helix fixed 
at atomic resolution.

Assuming Rosetta is in your home directory, you run it as follows:

    $ ~/Rosetta/main/source/bin/rosetta_scripts.default.[platform][compiler]release @flags -nstruct [number of structures]

where platform is "linux" or "mac" and compiler is "gcc" or "clang". 10000 
structures should be more than enough to see what's going on.
