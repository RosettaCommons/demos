Domain Insertion Protocol Capture
=================================

This protocol runs ab inito for a domain insertion protein per the Broker 
paper. Assuming Rosetta is in your home directory, you run it as follows:

    $ ~/Rosetta/main/source/bin/rosetta_scripts.default.[platform][compiler]release @flags -nstruct [number of structures]

where platform is "linux" or "mac" and compiler is "gcc" or "clang".

Note a difference from the paper is the "cheat_region" residue selector and 
"cheat" RigidChunkCM which are there to make the loop closure problem easier. 
Otherwise, it takes quite a long time, because we don't have constraints active 
during relaxation, and because scoring still scores the fixed regions, and the 
loop closure fails upwards of 80% of the time.
