# Rhiju Iron Chef Recipe
KEYWORDS: UTILITIES GENERAL
## Author
Rhiju Das, rhiju@stanford.edu

## Protocol Name
Rhiju's live demo on hacking code [creating a minimized helix for arbitrary sequence]

## Brief Description

  How not to code [renamed by Charlie & Rich to: Iron Chef Rosetta 1, Live Demo with Rotating Knives]
  Rhiju Das
  Tuesday Aug. 3, 2010

  A 20-minute attempt to demystify writing C++ protocols in mini. 
  

## Running
### Example Rosetta Command Line
```
~/src/mini/bin/protein_helix_assemble.macosgccrelease  -database ~/minirosetta_database/ -seq MRGSHHHHHHGMASIEGRGSLRDLQYALQEKIEELRQRDALIDELELELDQKDELIQMLQNELDKYRSVIRP
```

Above is a CASP target that is mostly helical (obviously the His-tag in the front is not, but hey I only had 20 minutes).


## Versions

The right Rosetta version:

https://svn.rosettacommons.org/source/branches/das_lab/mini
Revision: 36561


## Other Comments: 
The above is the Das lab branch, but I think the only thing I used from it is a function clear_conformation_viewers(), which you can comment out.

I think you can compile this in trunk if you just copy the file:
```
src/apps/pilot/rhiju/protein_helix_assemble.cc
```
to your pilot directory, and then make sure scons knows about it, in the file:
```
src/pilot_apps.src.settings.all
```
