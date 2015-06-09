# Template Demo Directory

The entire workflow for this demo should be described in a file named README.md, formatted in gollum markdown.  
It should describe an entire work flow, with command lines, tested if possible.

The contents of each demo directory should be:

```
starting_files/
  -- directory in which the raw input files are given - these
     are provided for you and serve as the base for your
     tutorial
rosetta_inputs/
  -- directory in which the modified starting files should
     be placed which will be used as inputs to Rosetta.
     You may need to make modifications like stripping
     extra chains from the input PDB; store the modified
     PDB here and leave the unaltered one in starting_files 
scripts/
  -- python scripts, shell scripts, used in the workflow
  -- awk, grep, sed, etc. command lines

README.dox
  -- A prose or list description of how to perform the protocol

FOR_AUTHORS.txt
  -- A description for the demo creators of what their demo
     should achieve.
  -- Most of what starts in this file should end up in the
     README file as well.
```
