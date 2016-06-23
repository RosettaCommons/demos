Protein protein interface design using rosetta scripts
======================================================

KEYWORDS: SCRIPTING_INTERFACES INTERFACES DESIGN

This demo has example script to show how a simple protien protein interface design works. The script is placed in the scripts directory and the input odb is in the inputs directory. Example output files are provided in the example_outputs directory.

To run the script, use:

(`$ROSETTA3`= path-to-Rosetta/main/source)
```
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @flags
```

Please note that in your system, depending on the compilation mode, you may need to change the extensions in the executable.
For more information about what each part of the script does, check [rosetta scripts documentation](https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts).

**NOTE** the values in the script are shorten to enable a fast demo run. You can increase them in multiple places (some mentioned in the xml script) in real application. In particular, you need to increase the MonteCarlo trial numbers.
