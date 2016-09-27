Flexible backbone design with RosettaScripts
============================================

KEYWORDS: SCRIPTING_INTERFACES DESIGN

This is a short demo showing how you can design your protein with flexible backbone using Rosetta Scripts. The flxbb mover can be given a blueprint file to extrapolate into a movemap; does design, can do some relaxation; converts your pose to alanine before designing.

To run the demo simple type:

(`$ROSETTA3`= path-to-Rosetta/main/source)

```
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @flags
```

Blueprint format:

        resnum  residue  (ss_struct)(abego) rebuild
        resnum = consecutive (starting from 1) or 0 (to indicate a new residue not in the input.pdb)
        residue = one letter code amino acid (e.g. V for Valine)
        ss_struct = secondary structure, E,L or H. ss_struct and abego are single-letter and have no space between them.
        abego = abego type (ABEGO), use X if any is allowed
        rebuild = R (rebuild this position) or "." (leave as is)
Examples

        1   V  LE  R   (position 1, Val, loop, abego type E, rebuild)
        0   V  EX  R   (insert a residue, Val, sheet, any abego, rebuild)
        2   V  EB  .   (position 2, Val, sheet, abego type B, do not rebuild)
