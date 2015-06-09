The entire workflow for this demo should be described in a file
named README.dox.   It should describe an entire work flow, with
command lines, tested if possible.

The contents of each demo directory should be:

starting_files/
  -- directory in which the raw input files are given
rosetta_inputs/
  -- directory in which the modified starting files should
     be placed which will be used as inputs to Rosetta.
scripts/
  -- python scripts, shell scripts, used in the workflow
  -- awk, grep, sed, etc. command lines

README.dox
  -- A prose

FOR_AUTHORS.txt
  -- A description for the demo creators of what their demo
     should achieve.
  -- Most of what starts in this file should end up in the
     README file as well.
     
@section Introduction and preliminary installation steps

The following tutorial will allow you to generate protein fragment files using the make_fragments.pl
script which can be found within the following directory:

       <code>https://svn.rosettacommons.org/source/branches/releases/rosetta-3.3/rosetta_fragments</code>
The directory provides also other necessary files, such as fragment database (vall.apr24.2008.extended.gz)
or a script usefull to convert between secondary structure prediction file formats (ss_pred_converter.py)

Note, that:
    -- All commands listed below assume use of the bash shell.
    -- It requires the use of third party software that is not included as part of the Rosetta distribution. In particular, 
    you will need to use secondary structure prediction software such as PsiPred, SAM, Jufo, or Porter.

@subsection Instructions for installation of PsiPred follow:

    -- Download PsiPred into using the following command:
	<code> curl –O http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/psipred321.tar.gz</code>

    -- PsiPred should then be extracted with the following command:
	<code>tar –xzf psipred321.tar.gz</code>

    -- You then need to compile PsiPred using the following commands, entered sequentially:
	<code>
	cd psipred/src/

	make

	make install
	</code>


If you would like to install other secondary structure prediction software, you will need to obtain and
install the standalone versions of the software yourself. This tutorial does not include information on
the use of these secondary structure prediction programs, but information about them can be found at
the following locations:

<bf>Porter:</bf> To obtain this software, you should send an email to gianluca.pollastri@ucd.ie inquiring about
the standalone version.

<bf>Sam:</bf> Sam can be obtained at: compbio.soe.ucsc.edu/sam2src/

<bf>Jufo:</bf> Information about Jufo can be found here: http://www.meilerlab.org/index.php/bclcommons/
show/b_apps_id/1

@subsection PsiBlast installation

You will then need to obtain and install PsiBlast. This can be done using the following commands:

<code>
curl –O ftp://ftp.ncbi.nlm.gov/blast/executables/release/LATEST/blast-*.tar.gz
</code>

It’s important that you download an older version of Blast called “legacy Blast” (not Blast+) because the
current version doesn’t write files in a format Rosetta can understand.

You should then also obtain and install the NCBI non-redundant (nr) database with the following
command:
<code>
curl –O ftp://ftp.ncbi.nih.gov/blast/db/nr0?.tar.gz
</code>

If prompted for a username and password, the username is anonymous and password is any email address.

Within the db dir (which will contain the nr.0x.tar.gz files) run the following command to uncompress the archives: 
<code> for i in *.tar; do tar –xzf $1; done </code>

@section Preparations to the <code>make_fragments.pl</code> script

Once the above has been done, you should open the script in order to modify its contents. The section
that must be modified is labeled “USER CONFIGURATION”, and points the script to the pertinent paths
to the database and programs you have just installed. You should change the paths in this section to
reflect the directories in which you have installed all pertinent programs / databases.

@section Actual fragment picking
All the steps listed above have to be done only once for a given machine, which now is ready to 
pick fragments. This demo will find fragments from the protein ubiquitin (PDB: 1ubq). The commandline entry for the
script is as follows, and should be run from the directory containing the script.

./make_fragments.pl -id 1ubqA -nosam -nojufo -noprof -nonnmake -picker -ref_struc ../starting_files/
1ubq.pdb.gz -nohoms ../starting_files/1ubq.fasta

The flags -nosam, and -nojufo, ensure that the script does not attempt to run SAM or Jufo. The -
picker flag signifies use of the fragment picker. The -ref_struc flag signifies the input .pdb file to be
used for analysis, and should include a correct path to the file. You will also need as an input a .fasta
corresponding to the sequence of your protein of interest. The -nohoms flag ensures that homologues
to the query sequence are not considered.


