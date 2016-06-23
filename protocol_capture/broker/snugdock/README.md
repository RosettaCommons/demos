SnugDock Protocol Capture
=========================
KEYWORDS: DESIGN DOCKING
This protocol runs our snugdock-inspired antibody docking protocol. It uses 
fragment insertion instead of normal loop modeling, but is otherwise very 
similar to the orginal.

Assuming Rosetta is in your home directory, you run it as follows:

    $ ~/Rosetta/main/source/bin/rosetta_scripts.default.[platform][compiler]release @flags -nstruct [number of structures]

where platform is "linux" or "mac" and compiler is "gcc" or "clang".
