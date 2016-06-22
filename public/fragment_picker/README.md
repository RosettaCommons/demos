### Standart fragment picking

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

 `$> <path-to-Rosetta>/main/source/bin/fragment_picker.linuxgccrelease @BestFragmentsProtocol/best-frags.flags -in::file::vall $ROSETTA3/main/database/sampling/small.vall.gz ` 

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
For other options as described in the [Rosetta Documentation](https://www.rosettacommons.org/docs/wiki/application_documentation/utilities/app-fragment-picker) have a look at the sub-directories of this demo.

```
./NullFragments  
./QuotaProtocol          
./Quota_with_restraints  
```

` 