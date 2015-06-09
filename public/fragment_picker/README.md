Do the best picker option for the starter demo: BestFragmentsProtocol 

Submit the following sequence (in the FASTA format) to 
[[psipred|http://bioinf.cs.ucl.ac.uk/psipred]]:

    > 2JSV X
    MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE

After psipred finishes get the "download the machine learning scores results in 
plain text format" and give it a meaningful name, e.g. `2jsv.psipred.ss2`. This 
should have a header with the string "VFORMAT"

* Put this vall in the flags:

        [path]/rosetta/rosetta_database/sampling/filtered.vall.dat.2006-05-05.gz

  Update the flags with your correct database location

* Run the picker like this:

        [path]/rosetta/rosetta_source/bin/fragment_picker.linuxgccrelease @best-frags.flags 

  Output:
  * `frags.200.3mers`: for use with `-in:file:frag3`
  * `frags.200.9mers`: for use with `-in:file:frag9`
  * fsc files for fragment information.

* Now get the quality of these fragments:

        [path]/rosetta/rosetta_source/bin/r_frag_quality.linuxgccrelease -database [path]/rosetta/rosetta_database -in:file:native input_files/2jsvX.pdb -f output_files/frags.200.9mers

  Output:
  * frag_qual.dat

  Use gnuplot to visualize the quality data:

        plot "frag_qual.dat" u 2:4 w p

  To plot this in plotting software (e.g. R or excel, or whatever you prefer) 
  do a scatter plot with column 2 as x-axis and col 4 as y-axis. 

More detailed information in:  
http://www.rosettacommons.org/manuals/archive/rosetta3.4_user_guide/dc/d10/app_fragment_picker.html

