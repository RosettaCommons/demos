# Ubiquitin Fasta Fragment Generation

The entire workflow for this demo should be described in a file named README.md.
It should describe an entire work flow, with command lines, tested if possible.

The contents of each demo directory should be:

```
starting_files/
  -- directory in which the raw input files are given
rosetta_inputs/
  -- directory in which the modified starting files should
     be placed which will be used as inputs to Rosetta.
scripts/
  -- python scripts, shell scripts, used in the workflow
  -- awk, grep, sed, etc. command lines

README.md
  -- A prose

FOR_AUTHORS.txt
  -- A description for the demo creators of what their demo
     should achieve.
  -- Most of what starts in this file should end up in the
     README file as well.
```
     
## Introduction and Preliminary Installation Steps

The following tutorial will allow you to generate protein fragment files using the `tools/fragment_tools/make_fragments.pl` script.
The directory provides also other necessary files, such as fragment database (vall.jul19.2011) or a script useful to convert between secondary structure prediction file formats (ss_pred_converter.py, also in `tools/fragment_tools/` dir).

**Note, that:**
- All commands listed below assume use of the bash shell.
- It requires the use of third party software that is not included as part of the Rosetta distribution. In particular, you will need to use secondary structure prediction software such as PsiPred, SAM, Jufo, or Porter.

### Instructions for Installation of PsiPred

- Download PsiPred into using the following command:
	`curl –O http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/psipred3.5.tar.gz`

- PsiPred should then be extracted with the following command:
	`tar –xzf psipred3.5.tar.gz`

- You then need to compile PsiPred using the following commands, entered sequentially:
	```
	cd psipred/src/

	make

	make install
	```


If you would like to install other secondary structure prediction software, you will need to obtain and
install the standalone versions of the software yourself. This tutorial does not include information on
the use of these secondary structure prediction programs, but information about them can be found at
the following locations:

**Porter:** To obtain this software, you should send an email to gianluca.pollastri@ucd.ie inquiring about
the standalone version.

**Sam:** Sam can be obtained at: compbio.soe.ucsc.edu/sam2src/

**Jufo:** Information about Jufo can be found here: http://www.meilerlab.org/index.php/bclcommons/
show/b_apps_id/1

**Note:** these may be out of date.

### PsiBlast Installation

You will then need to obtain and install PsiBlast. This can be done using the following commands:

`curl –O ftp://ftp.ncbi.nlm.gov/blast/executables/release/LATEST/blast-*.tar.gz`

It’s important that you download an older version of Blast called “legacy Blast” (not Blast+) because the
current version doesn’t write files in a format Rosetta can understand.

You should then also obtain and install the NCBI non-redundant (nr) database with the following
command:

`curl –O ftp://ftp.ncbi.nih.gov/blast/db/nr0?.tar.gz`

If prompted for a username and password, the username is anonymous and password is any email address.

Within the db dir (which will contain the nr.0x.tar.gz files) run the following command to uncompress the archives: 
`for i in *.tar; do tar –xzf $1; done`

## Preparations to the `make_fragments.pl` script

Once the above has been done, you should open the script in order to modify its contents. The section
that must be modified is labeled “USER CONFIGURATION”, and points the script to the pertinent paths
to the database and programs you have just installed. You should change the paths in this section to
reflect the directories in which you have installed all pertinent programs / databases.

## Actual Fragment Picking
All the steps listed above have to be done only once for a given machine, which now is ready to 
pick fragments. This demo will find fragments from the protein ubiquitin (PDB: 1ubq). The commandline entry for the
script is as follows, and should be run from the directory containing the script.

```
./make_fragments.pl -id 1ubqA -nosam -nojufo -noprof -nonnmake -picker -ref_struc ../starting_files/
1ubq.pdb.gz -nohoms ../starting_files/1ubq.fasta
```

The flags `-nosam`, and `-nojufo`, ensure that the script does not attempt to run SAM or Jufo. The -
picker flag signifies use of the fragment picker. The `-ref_struc` flag signifies the input .pdb file to be
used for analysis, and should include a correct path to the file. You will also need as an input a .fasta
corresponding to the sequence of your protein of interest. The `-nohoms` flag ensures that homologues
to the query sequence are not considered.


