## Working With Rosetta

There are three main interfaces for using Rosetta: the executables themselves, the XML-based RosettaScripts, and the Python-based PyRosetta.

The Rosetta executables themselves are run from the command line; see [here](http://structbio.vanderbilt.edu/comp/unix/ ) for information on how to run programs via the bash shell.

RosettaScripts are run via the rosetta_scripts executable, with the script itself passed via the "-parser_protocol" flag.

PyRosetta is accessed through python, via "from rosetta import \*" or similar.
