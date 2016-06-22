The resfile contains reversions back to the native sequence of the protein. In the Rosetta tools repository, there is a script that generates a template residue file where one can change residue back to native or test certain residues at different position.

```bash
$> python $ROSETTA_TOOLS/zinc_site_redesign/generate_residuefile.py rosetta_cst.pdb
```

The above returns a resfile with type NATAA, NATRO (keep rotamer) - one should
add PIKAA for picking new residue. For example:

```
    1 PIKAA WFY 
```
