Domain Insertion Protocol Capture
=================================
KEYWORDS: DESIGN GENERAL
This protocol runs ab inito for a domain insertion protein per the Broker 
paper. Because the loop closure for this protocol is particularly difficult,
and because the loop closure algorithm in this example relies on fragments
to perform the closure, the loop closure succeeds rather infrequently. The
problem is that, for many non-native loops, fragment coverage of the
required region of conformational space does not exist.

As a result, loop closure failed for me in >90% of cases. Thus, to test
this protocol, you must run with -n 100 or greater. Alternatively, the
closure is always able to close the loops in the native structure, so if
the sampling protocol starts near a correct (or probably any closable)
structure, so you can run with -s native.pdb to produce a closed structure
nearly every time.

An additional tip for producing output more quickly is to disable (or
turn down) the relaxation step at the end. Because relaxation is done
system-wide and not exclusively for the inserted domain, it can have
a very long runtime.

Assuming Rosetta is in your home directory, you run it as follows:

    $ ~/Rosetta/main/source/bin/rosetta_scripts.default.[platform][compiler]release @flags -nstruct [number of structures]

where platform is "linux" or "mac" and compiler is "gcc" or "clang".
