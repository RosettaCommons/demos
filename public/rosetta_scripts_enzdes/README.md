# Rosetta scripts: enzdes 
---------------------------------

KEYWORDS: SCRIPTING_INTERFACES DESIGN ENZYMES


The necessary files to run this demo is provided for you in the inputs directory. You need to make sure you change the files to your case.

To run the test:

(`$ROSETTA3`= path-to-Rosetta/main/source)

```
$> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @enzflags -parser:protocol scripts/enzdes.xml -s inputs/enzdes_test.pdb
```

An example out put is provided for you in the example_output directory.
 Please note that for a real run, you need to increase the cycles in each design module. The current demo is designed for fast performance.
Also, please refer to ligand tutorial, constraint tutorial and the [enzyme design documentation](https://www.rosettacommons.org/docs/latest/application_documentation/design/enzyme-design).
