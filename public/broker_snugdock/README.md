Toplogy Broker - SnugDock Protocol Capture
=========================

KEYWORDS: DOCKING GENERAL ANTIBODIES

This protocol runs our snugdock-inspired antibody docking protocol. It uses 
fragment insertion instead of normal loop modeling, but is otherwise very 
similar to the orginal.

Assuming Rosetta is in your home directory, you run it as follows (passing a much higher nstruct than this demo is running):

    $ ~/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease @flags -nstruct 1

where platform is "linux" or "mac" and compiler is "gcc" or "clang".
