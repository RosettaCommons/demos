### Standard fragment picking

Do the best picker option for the starter demo: 

`./BestFragmentsProtocol/` 

Submit the following sequence (in the FASTA format) to 
[[psipred|http://bioinf.cs.ucl.ac.uk/psipred]]:

    > 2JSV X
    MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE

After psipred finishes get the "download the machine learning scores results in 
plain text format" and give it a meaningful name, e.g. `2jsv.psipred.ss2`. This 
should have a header with the string "VFORMAT"

* Put this vall in the flags:

 `$ <path-to-Rosetta>/main/database/sampling/filtered.vall.dat.2006-05-05.gz`

  Update the flags with your correct database location

* Run the picker like this:

 `$> <path-to-Rosetta>/main/source/bin/fragment_picker.linuxgccrelease @BestFragmentsProtocol/flags -in::file::vall $ROSETTA3/main/database/sampling/small.vall.gz ` 

  Output:  
  
  * `frags.200.3mers`: for use with `-in:file:frag3`
  * `frags.200.9mers`: for use with `-in:file:frag9`
  * fsc files for fragment information.

* Now get the quality of these 3-mer fragments:

        $> <path-to-Rosetta>/main/source/bin/r_frag_quality.linuxgccrelease -in:file:native BestFragmentsProtocol/input_files/2jsvX.pdb -f BestFragmentsProtocol/output_files/frags.200.3mers

  Output:
  * frag_qual.dat

  Use gnuplot to visualize the quality data:

        plot "frag_qual.dat" u 2:4 w p

  This command simply makes a scatter plot with column 2 on the x-axis and 
  column 4 on the y-axis. 

More detailed information in:  
[Rosetta Documentation](https://www.rosettacommons.org/docs/wiki/application_documentation/utilities/app-fragment-picker)

### Tweak the fragment picker

There many ways to influence how fragments are being picked.
Many options are described in the [Rosetta Documentation](https://www.rosettacommons.org/docs/wiki/application_documentation/utilities/app-fragment-picker).

The various subdirectories of this demo serve as repository for examples:

```
./NullFragments  
./fragment_picking_for_flexible_loop_design
./fragment_picking_torsion_class_score
./fragment_picking_with_psi_jufo_sam_L1_quota
./fragment_picking_with_quota
./fragment_picking_with_quota_and_restraints  
```

You can test them:

```
$> <path-to-Rosetta>/main/source/bin/fragment_picker.linuxgccrelease @NullFragments/flags -in::file::vall $ROSETTA3/main/database/sampling/small.vall.gz
$ <path-to-Rosetta>/main/source/bin/fragment_picker.linuxgccrelease @fragment_picking_for_flexible_loop_design/flags -in::file::vall $ROSETTA3/main/database/sampling/small.vall.gz
$> <path-to-Rosetta>/main/source/bin/fragment_picker.linuxgccrelease @fragment_picking_torsion_class_score/flags -in::file::vall $ROSETTA3/main/database/sampling/small.vall.gz
$> <path-to-Rosetta>/main/source/bin/fragment_picker.linuxgccrelease @fragment_picking_with_psi_jufo_sam_L1_quota/flags -in::file::vall $ROSETTA3/main/database/sampling/small.vall.gz
$> <path-to-Rosetta>/main/source/bin/fragment_picker.linuxgccrelease @fragment_picking_with_quota/flags -in::file::vall $ROSETTA3/main/database/sampling/small.vall.gz
$> <path-to-Rosetta>/main/source/bin/fragment_picker.linuxgccrelease @fragment_picking_with_quota_and_restraints/flags -in::file::vall $ROSETTA3/main/database/sampling/small.vall.gz
```` 

The fragment files can be found in <whatever version you ran>/output_files. There are example outputs already. The newly generated ones will contains ".200". So, you can check whether it worked.
