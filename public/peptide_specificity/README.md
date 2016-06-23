# Peptide Specificity

KEYWORDS: PEPTIDES DOCKING

## What is this?
This is the demo for the pepspec and pepspec_anchor_dock applications.

## What does it do?
First, a single peptide "anchor" residue will be docked onto the surface of a peptide binding protein (c-CRK SH3 domain) in which the original peptide has been deleted. This docking will be repeated to generate an ensemble of peptide anchor-docked structures. To do this, the pepspec_anchor_dock application uses the relative orientation of peptide anchor residues in three homologous peptide complex structures (1CKB, 1N5Z, 1OEB).
Next, these anchor residue docked structures are used by the pepspec application to design putative binding peptides on the surface of the SH3 domain. Peptide residues are added to the anchor residue while peptide sequences and structures are simultaneously explored. Structures and sequences are saved for later processing.

# Demo Summary 
1. generating an ensemble of anchor prolines docked to an SH3 scaffold, then

2. exploring the sequence specificity of peptides designed from those docked proline anchors

3. generate a specificity PWM from the designed peptides

# Demo Detail
1. run the anchor docking protocol to generate 10 structures of proline docked to 1CKA.align.nopep.pdb
    ```
    ~/mini/bin/pepspec_anchor_dock.linuxgccrelease /path/to/minirosetta_database @dock.args
    ```

    This will generate 10 pdbs, a pdblist file, and a cst (contraint) file. These files are used in the next stage.

    - 1cka.docked_[1-10].pdb are the anchor-docked structures of the input structure 1CKA.align.nopep.pdb
    - 1cka.docked.pdblist is simply a list of the pdb files generated
    - 1cka.docked.cst is a peptide constraint file

2. run the peptide design protocol:
    ```
    ~/mini/bin/pepspec.linuxgccrelease -database /path/to/minirosetta_database @spec.args
    ```

    This will generate a folder "1cka_spec.pdbs" with some designed peptide - protein complexes and
a 1cka.spec file with all the peptide sequences, energies, and some other values.

    This *.spec file can be used for deriving a PWM using the script

3. Genrate a normalized PWM (position-weight-matrix)
    ```
    ~/mini/analysis/apps/gen_pepspec_pwm.py 1cka_spec.spec 3 0.1 binding-prot_score /path/to/minirosetta_database/pepspec_background.binding-prot-0.1.pwm
    ```

    This will sort the peptide sequences in 1cka_spec.spec, in which there are 3 residues n-term to the anchor residue, by the "binding-prot_score" score term (one of many scores calculated in 1cka_spec.spec), then filter out the lowest-scoring 10% peptide sequences, and use these sequences to construct a peptide PWM based on the position-specific amino acid frequencies. This resultant PWM is then normalized by a background PWM (pepspec_background.binding-prot-0.1.pwm) in order to eliminate some residue frequency biases caused by artifacts of the rosetta score function.

    The script creates a PWM file (1cka_spec.binding-prot_score.0.1.pwm), a list of the sequences from which that PWM was constructed (1cka_spec.binding-prot_score.0.1.seq), and a normalized PWM file (1cka_spec.binding-prot_score.0.1.norm.pwm).

    For a real, production-level run, I would reset the values of these flags:
    ```
    <dock.args>
    -pepspec::n_peptides 100
    
    <spec.args>
    -pepspec::n_peptides 1000
    -pepspec::diversify_lvl 5
    ```

**Note:**
	In a real production-level run, the peptide backbone coordinate constraints (e.g. 1cka.docked.cst) may be generated from only one or two homologue complexes. If this is the case, or you have reason to believe that peptides may bind with backbone conformations not represented in the homologue set, then use of the constraints may prevent your simulation from ever generating the correct peptide backbone structures.
	To avoid use of constraints, simply do not include a reference to the constraint file in you command line arguments. (e.g. delete the line "-pepspec::homol_csts 1cka.docked.cst" from spec.args).
	If you're not using constraints, you need to do much, MUCH more sampling (i.e. an order of magnitude more).
The unconstrained sequence+structure space gets big fast.

