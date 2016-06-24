The script will analyze the zinc site in the PDB file and print type
of metal site. 
KEYWORDS: METALS DESIGN
Input parameter are PDB file and type of metal ion - it has mainly
been tested with Zinc ions. 

Example with 1A4L. 

````bash
python $ROSETTA_TOOLS/zinc_site_redesign/analyze_zinc_site.py -f 1A4L.pdb
```

Output from analysis

```
Number of chains in pdb file:  4
Number of interactions between metal and protein: 4
Number of interactions between hetero atom : 1
PDB chain D contains 5 coordinated metal site
Number of interactions between metal and protein: 4
Number of interactions between hetero atom : 1
PDB chain C contains 5 coordinated metal site
Number of interactions between metal and protein: 4
Number of interactions between hetero atom : 1
PDB chain A contains 5 coordinated metal site
Number of interactions between metal and protein: 4
Number of interactions between hetero atom : 1
PDB chain B contains 5 coordinated metal site
```

To summarize: 4 chains with a triangular bipyramidal zinc site.
