The resfile contains reversions back to the native sequence of the protein. In the Rosetta tools repository, there is a script that generates a template residue file where one can change residue back to native or test certain residues at different position.
KEYWORDS: METALS UTILITIES
```bash
python $ROSETTA_TOOLS/zinc_site_redesign/generate_residuefile.py UM_1_H15H17H214D295Q58_1A4L_clean_A_r_1A4L_clean_A_1__DE_1.pdb
```

The above returns a resfile with type NATAA, NATRO (keep rotamer) - one should
add PIKAA for picking new residue. For example:

```
    1 PIKAA WFY 
```
