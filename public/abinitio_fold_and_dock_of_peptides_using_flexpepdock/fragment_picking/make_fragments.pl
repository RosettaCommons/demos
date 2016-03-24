#!/usr/bin/perl

#  CVS information:
#  $Revision: 19024 $
#  $Date: 2007-12-11 01:01:36 +0200 (Tue, 11 Dec 2007) $
#  $Author: tex $

###############################################################################
# USER CONFIGURATION  -- ONLY CHANGE THIS SECTION.
#                        YOU MUST EDIT THESE PATHS BEFORE USING
###############################################################################

# change the following paths to point to the locations of your copies of these
# files, databases or directories.

# tail of fragment files
#
# v1_4 : new tail for 4ss quota: .6, .4, .4, .2 jones, rdb, prof, juf

my $TAIL = "_v1_3";

# shareware_dir
#
# NOTE: set $BAKER_HOME in your environment to be your home directory,
#   or where ever you installed the "src" directory to, so then you have
#   various shareware apps in $BAKER_HOME/src/shareware/ and nnamke at
#   $BAKER_HOME/src/nnmake/
#
#  e.g. in .cshrc:  setenv BAKER_HOME /users/dylan
#       in .bashrc: export BAKER_HOME=/users/dylan
#
#if (! defined $ENV{'BAKER_HOME'}) {
#    print STDERR ("must define \$BAKER_HOME in environment\n");
#    exit -1;
#}

#make sure all paths to local software are updated (psipred,blast,etc...)
#note that the "NNMAKE" script is found in the demo dir.

#$shareware_dir = $src_dir . '/shareware/';
my $fragsDir   = "/vol/ek/share/rosetta/rosetta_fragments_latest";
my $blastDir   = "/vol/ek/share/bin/blast-2.2.18/bin";
my $nrDir      = "/vol/ek/share/databases/genomes/nr_2006_12";
my $psiPredDir = "/vol/ek/share/src/psipred";
my $samDir     = "/vol/ek/share/src/sam3.5.i686-linux/bin";

# psi-blast
my $PSIBLAST      = "$blastDir/blastpgp";                              # PSI-BLAST (duh.)
my $NR            = "$nrDir/nr";                                       # nr blast database filename
my $VALL_BLAST_DB = "$fragsDir/nnmake_database/vall.blast.2006-05-05"; # vall blast database filename (cvs respository 'nnmake_database')
my $BLOSUM        = "$fragsDir/nnmake_database/";                      # BLOSUM score matrices directory (cvs repository 'nnmake_database')

# psipred
my $FILTNR       = "$nrDir/filtnr";            # filtnr blast database filename
my $MAKEMAT      = "$blastDir/makemat";        # makemat utility (part of NCBI tools)
my $PSIPRED      = "$psiPredDir/bin/psipred";  # psipred
my $PSIPASS2     = "$psiPredDir/bin/psipass2"; # psipass2 (part of psipred pkg)
my $PSIPRED_DATA = "$psiPredDir/data";         # dir containing psipred data files.

# nnmake
my $VALL      = "$fragsDir/nnmake_database";    # dir containing vall database (cvs repository 'nnmake_database')
my $VALL2     = "$fragsDir/nnmake_database";    # alt dir containing vall database (cvs repository 'nnmake_database')
my $VALL_NAME = "vall.dat.2006-05-05";          # filename of vall (vall.dat.<id> and vall_cst_coord.dat.<id> must exist)
my $NNMAKE    = "/vol/ek/liorz06/ProtocolCapture/FlexPepDock_AbInitio/scripts/frags/make.sh";                         # nnmake binary  (cvs respository 'nnmake')
my $TRIMLOOPS = "$fragsDir/nnmake/trimLoopLibrary.pl";                  # trimLoopLibrary.pl (cvs respository 'nnmake')

# chemshift
my $CHEMSHIFT = "$fragsDir/chemshift/pCHEMSHIFT.gnu"; # chemshift binary (cvs repository 'chemshift')
my $TALOS_DB = "$fragsDir/chemshift_database";        # TALOS databases directory (cvs respository 'chemshift_database')

# jufo (secondary structure prediction software)
my $JUFO = "/vol/ek/share/src/jufo/molecule.exe";                              # jufo executable

# sam (secondary structure prediction software)
my $SAM_target99		= "$samDir/target99";	# sam target99 executable
my $SAM_uniqueseq		= "$samDir/uniqueseq";	# sam uniqueseq executable
my $SAM_predict_2nd_dir = "/vol/ek/share/src/sam.predict-2nd/";			# sam predict-2nd directory
my $SAM_predict_2nd		= "$SAM_predict_2nd_dir/predict-2nd";	# sam predict-2nd executable

###############################################################################
#
# MAKE_FRAGMENTS.PL 1.00 -- THE (PEN)ULTIMATE IN HOME FRAGMENT-PICKING SOFTWARE!
#
# CAUTION:  NO USER SERVICEABLE PARTS BELOW!
#
#           TO REDUCE RISK OF ELECTRIC SHOCK, DO NOT REMOVE THE COVER!
#           DO NOT ATTEMPT REPAIRS!  REFER SERVICING TO YOUR AUTHORIZED DEALER!
#           AVOID PROLONGED EXPOSURE TO HEAT OR SUNLIGHT!
#           TO REDUCE THE RISK OF FIRE OR ELECTRIC SHOCK, DO NOT EXPOSE THE
#            PRODUCT TO RAIN AND/OR MOISTURE!
#           DO NOT MOVE THE PRODUCT WHILE IN USE!
#           DO NOT LOOK AT THE PRODUCT WHILE IN USE!
#           DO NOT COMPLAIN ABOUT THE PRODUCT WHILE IN USE!
#           DO NOT DISCUSS THE PRODUCT WHILE IN USE!
#           DO NOT THINK ABOUT THE PRODUCT WHILE IN USE!
#           CLEAN ONLY WITH MILD DETERGENTS AND A SOFT CLOTH!
#           USE ONLY IN WELL-VENTILATED AREAS!
#
#           FOR EXTERNAL USE ONLY!  DO NOT TAKE INTERNALLY!
#           MAY PRODUCE STRONG MAGNETIC FIELDS!
#
#           DO NOT REMOVE THIS TAG UNDER PENALTY OF LAW.
#
#           THIS ARTICLE CONTAINS NEW MATERIAL ONLY.
#
#           THIS LABEL IS AFFIXED IN COMPLAINCE WITH THE UPHOLSTERED AND
#            STUFFED ARTICLES ACT.
#
#
# (IN OTHER WORDS:  DON'T EVEN *THINK* ABOUT CHANGING THINGS BELOW THIS POINT!)
#
###############################################################################

use Cwd;
use bytes;

$DEBUG = 0;
$VERSION = 1.00;

$| = 1;                                              # disable stdout buffering

###############################################################################
# init
###############################################################################

# argv
my %opts = &getCommandLineOptions ();
my $file = $opts{f};
my $run_dir =  cwd();  # get the full path (needed for sam stuff)
my $no_homs = 0;
my $current_file = "";
my $psipred_file = "";
my $sam_file = "";
my $phd_file = "";
my $jufo_file = "";
my $psipred = 1;
my $psipred_iter = 1;
my $psipred_hbias = 1;
my $psipred_sbias = 1;
my $jufo = 1;
my $sam = 1;
my $prof = 1;
my $id = "temp";;
my $chain = "_";
my $xx = "aa";
my @checkpoint_matrix;
my $sequence;
my $num_ss_preds = 0;
my $min_num_ss_preds = 1;
my $cleanup = 1;
my @cleanup_files = ();
my $template = "";
my $barcodefile = "";
my $make_loops = 0;

chop($run_dir) if (substr($run_dir, -1, 1) eq '/');

# wait a random amount of time before proceeding, to avoid disk-choke problems
srand();
sleep(int(rand(10))+1);


if (defined($opts{verbose})) {
  if ($opts{verbose} == 1) {
    $DEBUG = 1;
    print "Run options:\n";
    print "be verbose.\n";
  }
}

if (defined($opts{rundir})) {
  $run_dir = $opts{rundir};
  &checkExist('d',$run_dir);
  chop($run_dir) if (substr($run_dir, -1, 1) eq '/');

  (!$DEBUG) || print "run directory: $run_dir\n";
}

if (defined($opts{homs})) {
  if ($opts{homs} == 0) {
    $no_homs = 1;
    (!$DEBUG) || print "exclude homologs.\n";
  }
}

if (defined($opts{psipred})) {
  if ($opts{psipred} == 0) {
    $psipred = 0;
    (!$DEBUG) || print "don't run psipred.\n";
  }
}

if (defined($opts{prof})) {
  if ($opts{prof} == 0) {
    $prof = 0;
    (!$DEBUG) || print "don't run prof.\n";
  }
}

if (defined($opts{psipred_iter})) {
  $psipred_iter = $opts{psipred_iter};
  (!$DEBUG) || print "psipred_iter = $psipred_iter\n";
}

if (defined($opts{psipred_hbias})) {
  $psipred_hbias = $opts{psipred_hbias};
  (!$DEBUG) || print "psipred_hbias = $psipred_hbias\n";
}

if (defined($opts{psipred_sbias})) {
  $psipred_sbias = $opts{psipred_sbias};
  (!$DEBUG) || print "psipred_sbias = $psipred_sbias\n";
}

if (defined($opts{jufo})) {
  if ($opts{jufo} == 0) {
    $jufo = 0;
    (!$DEBUG) || print "don't run jufo.\n";
  }
}

if (defined($opts{sam})) {
  if ($opts{sam} == 0) {
    $sam = 0;
    (!$DEBUG) || print "don't run sam.\n";
  }
}

if (defined($opts{xx})) {
  $xx = substr($opts{xx}, 0, 2);
  (!$DEBUG) || print "nnmake xx code: $xx\n";
}

if (defined($opts{cleanup})) {
  if ($opts{cleanup} == 0) {
    $cleanup = 0;
  }
}

if (defined($opts{template})) {
     $template = $opts{template};
     $make_loops = 1;
}

if (defined($opts{barcode_file})) {
     $barcodefile = $opts{barcode_file};
}

(!$DEBUG) || print "FILENAME: $file\n";

if (!defined($opts{id})) {
  (!$DEBUG) || print "no id specified. parse filename instead.\n";

  ($id) = $file =~ /(\w+\.\w+)$/;
  (!$DEBUG) || print "INTERMEDIATE: $id\n";
  ($id) = $id =~ /^(\w+)/;

  if (length($id) != 5) {
    die("DANGER WILL ROBINSON! DANGER! Your fasta filename is more than/less than five letters!\n".
	"You should either 1) rename your file to something of the form *****.fasta -or-\n".
	"You should explicitly specify the id and chain with the -id option.\n");
  }


  $chain = substr($id, 4, 1);
  $id = substr($id, 0, 4);
  (!$DEBUG) || print "ID: $id CHAIN: $chain\n";

} else {
  chomp $opts{id};

  (!$DEBUG) || print "id specified by user: $opts{id}\n";

  if (length($opts{id}) != 5) {
    die("The id you specify must be 5 characters long.\n");
  }

  if ($opts{id} =~ /\W+/) {
    die("Only alphanumeric characters and _ area allowed in the id.\n");
  }

  $id = substr($opts{id}, 0, 4);
  $chain = substr($opts{id}, 4, 1);

  (!$DEBUG) || print "ID: $id CHAIN: $chain\n";
}

if (defined($opts{minss})) {
  $min_num_ss_preds = $opts{minss};
  (!$DEBUG) || print "minimum # required ss predictions: $min_num_ss_preds\n";
}

#################################################
# Verify existence of programs and other files. #
#################################################

if ( $DEBUG ) {
    my @programs = (
        $PSIBLAST,
        $PSIPRED,
        $PSIPASS2,
        $NNMAKE,
        $TRIMLOOPS,
        $CHEMSHIFT,
        $JUFO,
        $SAM_target99,
        $SAM_uniqueseq,
        $SAM_predict_2nd_dir,
        $SAM_predict_2nd
    );

    my @files = ($NR, $VALL_BLAST_DB);

    foreach my $prog (@programs) {
        if ( ! -x $prog ) { warn "Error: can't find executable file $prog!\n"; }
    }
    foreach my $file (@files) {
        if ( ! -f $file ) { warn "Error: can't find file $file!\n"; }
    }
}


#########
# determine what ss predictions to run
#########

if (defined($opts{psipredfile})) {
    $psipred_file = $opts{psipredfile};

    if (&fileExist($psipred_file)) {
	if ($psipred) {
	    print "Whoa there, pokey...you've specified -psipredfile without specifying -nopsipred...\n";
	    print "I found your precious psipred file, so I'm NOT going to run psipred...setting -nopsipred...\n";
	    $psipred = 0;
	}
	($psipred_file_base)= $psipred_file =~ /\/*([\w\.]+)$/;
	if ( "$run_dir/$psipred_file_base" ne $psipred_file) {
	    `cp $psipred_file $run_dir/$psipred_file_base`;
	    push(@cleanup_files, $psipred_file_base);
	}
    } else {
	print "Specified psipred file ($psipred_file) not found!  Running psipred instead.\n";
	$psipred = 1;
    }
}

if ( &fileExist("$run_dir/$id$chain.psipred_ss2")) {
    $psipred=0;
    print "Assuming $run_dir/$id$chain.psipred_ss2 is a psipred_ss2 file.\n";

}

#########
if (defined($opts{jufofile})) {
    $jufo_file = $opts{jufofile};

    if (&fileExist($jufo_file)) {
	if ($jufopred) {
	    print "Whoa there, pokey...you've specified -jufofile without specifying -nojufo...\n";
	    print "I found your precious jufo file, so I'm NOT going to run jufo...setting -nojufo...\n";
	    $jufo = 0;
	}
	($jufo_file_base)= $jufo_file =~ /\/*([\w\.]+)$/;
	if ( "$run_dir/$jufo_file_base" ne $jufo_file) {
	    `cp $jufo_file $run_dir/$jufo_file_base`;
	    push(@cleanup_files, $jufo_file_base);
	}
    } else {
	print "Specified jufo file ($jufo_file) not found!  Running jufo instead.\n";
	$jufo = 1;
    }
}

if ( &fileExist("$run_dir/$id$chain.jufo_ss")) {
    $jufo=0;
    print "Assuming $run_dir/$id$chain.jufo_ss is a jufo file.\n";

}
###############

if (defined($opts{samfile})) {
    $sam_file = $opts{samfile};

    if (&fileExist($sam_file)) {
	if ($sampred) {
	    print "Whoa there, pokey...you've specified -samfile without specifying -nosam...\n";
	    print "I found your precious sam file, so I'm NOT going to run sam...setting -nosam...\n";
	    $sam = 0;
	}
	($sam_file_base)= $sam_file =~ /\/*([\w\.]+)$/;
	if ( "$run_dir/$sam_file_base" ne $sam_file) {
	    `cp $sam_file $run_dir/$sam_file_base`;
	    push(@cleanup_files, $sam_file_base);
	}
    } else {
	print "Specified sam file ($sam_file) not found!  Running sam instead.\n";
	$sam = 1;
    }
}

if ( &fileExist("$run_dir/$id$chain.rdb")) {
    $sam=0;
    print "Assuming $run_dir/$id$chain.rdb is a sam file.\n";
}

##############
if (defined($opts{phdfile})) {
    $phd_file = $opts{phdfile};

    if (&fileExist($phd_file)) {
	($phd_file_base)= $phd_file =~ /\/*([\w\.]+)$/;
	if ( "$run_dir/$phd_file_base" ne $psipred_file) {
	    `cp $psipred_file $run_dir/$phd_file_base`;
	    push(@cleanup_files, $phd_file_base);
	}
    } else {
	print "Specified phd file ($phd_file) not found!\n";
    }
}
if ( &fileExist("$run_dir/$id$chain.phd")) {
    $sam=0;
    print "Assuming $rundir/$id$chain.phd is a phd file.\n";
}

###############################################################################
# main
###############################################################################

chdir($run_dir);

# get the sequence from the fasta file
open(SEQFILE, $file);

my $has_comment = 0;
my $has_eof	= 0;
my $eof;
while (<SEQFILE>) {
  $eof = $_;
  s/\s//g;
  if (/^>/) { $has_comment = 1; next; }
  chomp;
  $sequence .= $_;
}
$has_eof = 1 if ($eof =~ /\n$/);
($has_comment && $has_eof) or die "fasta file must have a comment and end with a new line!\n";

(!$DEBUG) || print "Sequence: $sequence\n";

# run blast
unless (&fileExist("$id$chain.check" ))  {
  if (!&try_try_again("$PSIBLAST -t T -i $file -F F -j2 -o $id$chain.blast -d $NR -v10000 -b10000 -K1000 -h0.0009 -e0.0009 -C $id$chain.check -Q $id$chain.pssm",
		      2, ["$id$chain.check"], ["$id$chain.check","$id$chain.blast","$id$chain.pssm","error.log"])) {
    die("checkpoint psi-blast failed!\n");
  }
}

unless (&fileExist("$id$chain.checkpoint")) {
  # parse & fortran-ify the checkpoint matrix.
  @checkpoint_matrix = &parse_checkpoint_file("$id$chain.check");
  @checkpoint_matrix = &finish_checkpoint_matrix($sequence, @checkpoint_matrix);
  &write_checkpoint_file("$id$chain.checkpoint", $sequence, @checkpoint_matrix);
}

push(@cleanup_files,("$id$chain.blast","$id$chain.pssm","error.log"));

#############################################
# Secondary Structure Prediction methods
#############################################

#
# preliminaries -- run psi-blast for jufo and/or psipred
#
if ($psipred || $jufo) {
  unless (&fileExist("sstmp.chk") && &fileExist("sstmp.ascii") && &fileExist("ss_blast")) {
      ## PB - add -b -v args since ss_blast can be used by PROF
    if (!&try_try_again("$PSIBLAST -t T -b10000 -v10000 -j3 -h0.001 -d $FILTNR -i $file -C sstmp.chk -Q sstmp.ascii -o ss_blast",
			2, ["sstmp.chk","sstmp.ascii"], ["sstmp.chk","sstmp.ascii","ss_blast"])) {
      die("jufo/psipred psi-blast failed!\n");
    }
    push(@cleanup_files,("ss_blast","sstmp.chk","sstmp.ascii"));
  }
}

##########
#
# psipred
#
if ($psipred) {
  (!$DEBUG) || print "running psipred.\n";

  unless (&fileExist("psitmp.sn")) {
      &run("echo $file > psitmp.sn",("psitmp.sn"));
  }
  unless (&fileExist("psitmp.pn")) {
      &run("echo sstmp.chk > psitmp.pn", ("psitmp.pn"));
  }

  unless (&fileExist("sstmp.mtx")) {
      if (!&try_try_again("$MAKEMAT -P psitmp", 2, ["sstmp.mtx"], ["sstmp.mtx"])) {
	  die("psipred: makemat failed.");
      }
  }

  unless (&fileExist("psipred_ss")) {
      if (!&try_try_again("$PSIPRED sstmp.mtx $PSIPRED_DATA/weights.dat $PSIPRED_DATA/weights.dat2 $PSIPRED_DATA/weights.dat3 $PSIPRED_DATA/weights.dat4 > psipred_ss",
			  2, ["psipred_ss"], ["psipred_ss"])) {
	  die("psipred failed.");
      }
  }

  unless (&fileExist("psipred_horiz")) {
      print "DEBUG\n";
      if (!&try_try_again("$PSIPASS2 $PSIPRED_DATA/weights_p2.dat $psipred_iter $psipred_hbias $psipred_sbias psipred_ss2 psipred_ss > psipred_horiz",
			  2, ["psipred_ss2","psipred_horiz"], ["psipred_ss2","psipred_horiz"])) {
	  die("psipred/psipass2 failed.");
      }

      rename("psipred_horiz", "$id$chain.psipred") || die ("couldn't move psipred_horiz to $id$chain.psipred: $!\n");
#      rename("psipred_ss2", "$id$chain.psipred_ss2") || die ("couldn't move psipred_ss2 to $id$chain.psipred_ss2: $!\n");
      &run("echo '# PSIPRED VFORMAT (PSIPRED V2.6 by David Jones)' > tmp", ("tmp") );
      &run("cat tmp psipred_ss2 > $id$chain.psipred_ss2", ("$id$chain.psipred_ss2") );
  }

  if (&fileExist("$id$chain.psipred_ss2")) {
      (!$DEBUG) || print "psipred file ok.\n";
  } else {
      print "psipred run failed!\n";
  }

  push(@cleanup_files,(<psitmp*>,"psipred_ss","sstmp.mtx"));

}

#
# jufo  -- jens meiler's ss prediction software
#
if ($jufo) {
    (!$DEBUG) || print "running jufo.\n";

    open(JUFOCMD, ">jufo.input") || die("couldn't create jufo command file: $!\n");
    print JUFOCMD "readfasta $file\n";
    print JUFOCMD "readblast sstmp.ascii\n";
    print JUFOCMD "calcsecondary -n -b\n";
    print JUFOCMD "writesecondary $id$chain.jufo_ss -p\n";
    print JUFOCMD "quit\n";
    close(JUFOCMD) || die("couldn't create jufo command file: $!\n");

    if (!&try_try_again("$JUFO < jufo.input", 2, ["$id$chain.jufo_ss"], ["$id$chain.jufo_ss"])) {
	die("jufo failed!\n");
    }

    push(@cleanup_files,("jufo.input","error.log"));
}

#
# sam -- target99 alignment and predict2nd with 6 state neural net - condensed output to 3 state
#
if ($sam) {
    (!$DEBUG) || print "running sam.\n";

    my $target99_out = $id.$chain.".target99";
    my $target99_a2m_file = $target99_out.".a2m";

    ## run target99 for a2m alignment
    if (!&try_try_again("$SAM_target99 -seed $file -out ".$target99_out, 2, [$target99_a2m_file], [$target99_out_a2m_file])) {
	die "sam target99 failed!\n";
    }

    ## run uniqueseq
    my $uniqueseq_a2m_id	= $id.$chain.".uniqueseq";
    my $uniqueseq_a2m_file	= $uniqueseq_a2m_id.".a2m";
    if (!&try_try_again("$SAM_uniqueseq $uniqueseq_a2m_id -percent_id 0.9 -alignfile ".$target99_a2m_file, 2, [$uniqueseq_a2m_file], [$uniqueseq_a2m_file])) {
	die "sam uniqueseq failed!\n";
    }

    ## run predict-2nd
    chop($SAM_predict_2nd_dir) if (substr($SAM_predict_2nd_dir, -1, 1) eq '/');

    # create samscript
    my $sam_6state = "$run_dir/$id$chain.sam_6state";
    my $sam_ebghtl = "$run_dir/$id$chain.sam_ebghtl";
    my $sam_log    = "$run_dir/$id$chain.sam_log";
    my $samscript_txt_buf =
	qq{ReadAlphabet $SAM_predict_2nd_dir/std.alphabet
	       ReadAlphabet $SAM_predict_2nd_dir/DSSP.alphabet
		   ReadNeuralNet $SAM_predict_2nd_dir/overrep-3617-IDaa13-7-10-11-10-11-7-7-ebghtl-seeded3-stride-trained.net
		       ReadA2M $run_dir/$uniqueseq_a2m_file
			   PrintRDB $sam_6state
			       PrintPredictionFasta $sam_ebghtl
			       };
    my $samscript_txt = "$run_dir/$id$chain.samscript.txt";
    open  (SAMSCRIPT, '>'.$samscript_txt);
    print  SAMSCRIPT $samscript_txt_buf;
    close (SAMSCRIPT);

    # get into $SAM_predict_2nd_dir so sam can read file recode3.20comp
    chdir $SAM_predict_2nd_dir;

    if (!&try_try_again("$SAM_predict_2nd -noalph < $samscript_txt >& $sam_log", 2, [$sam_6state], [$sam_6state])) {
	die "sam predict-2nd failed!\n";
    }

    chdir $run_dir;

    # condense to 3 state prediction
    # formerly a separate perl script: condense_6state.pl
    # integrated here 8/1/03   C. Rohl

    my $rows;
    open IN, "<$sam_6state" || die "Cannot open file $sam_6state: $!\n";
    @{$rows} = <IN>;
    close (IN);
    $rows = &condense_rdb6_rdb($rows);
    open OUT, ">$id$chain.rdb" || die "Cannot open file $id$chain.rdb: $!\n";
    print OUT @{$rows};
    close (OUT);


    (!$DEBUG) || print "sam file ok.\n";

    push(@cleanup_files,($samscript_txt, $sam_log, $sam_6state,$sam_ebghtl,$target99_a2m_file,
                         $uniqueseq_a2m_file,"$target99_out.cst"));

}

#############################################
# Chem shift data stuff
#############################################

# run pCHEMSHIFT if there's a .chsft_in file in the run directory..
if (&fileExist("$id$chain.chsft_in") && !&fileExist("$id$chain.chsft")) {
  if (!&try_try_again("$CHEMSHIFT $id$chain $TALOS_DB/", 2, ["$id$chain.chsft"], ["$id$chain.chsft"])) {
    die(".chsft_in file found in $run_dir, and $CHEMSHIFT failed!\n");
  }
}

#############################################
# Vall and homologue searches
#############################################

if ($no_homs) {
  unless (&fileExist("$id$chain.outn")) {
    if (!&try_try_again("$PSIBLAST -t T -i $file -j 1 -R $id$chain.check -o $id$chain.outn -e 0.05 -d $VALL_BLAST_DB", 2, ["$id$chain.outn"], ["$id$chain.outn"])) {
      die("homolog blast failed!\n");
    }
  }

  &exclude_blast($id,$chain);
  &exclude_outn($id,$chain);
}

#############################################
# nnmake
#############################################

open(PATH_DEFS, ">path_defs.txt");

# make the path_defs file for this run.
print PATH_DEFS "./\n";
print PATH_DEFS "./\n";
print PATH_DEFS "./\n";
print PATH_DEFS "$VALL2/\n";
print PATH_DEFS "./\n";
print PATH_DEFS "./\n";
print PATH_DEFS "./\n";
print PATH_DEFS "./\n";
print PATH_DEFS "./\n";
print PATH_DEFS "./\n";
print PATH_DEFS "./\n";
print PATH_DEFS "./\n";
print PATH_DEFS "$VALL/\n";
print PATH_DEFS "./\n";
print PATH_DEFS "./\n";
print PATH_DEFS "./\n";
print PATH_DEFS "./\n";
print PATH_DEFS "$VALL_NAME\n";
close(PATH_DEFS);

push(@cleanup_files,"path_defs.txt");

unless (-f <$xx$id$chain*$TAIL>) {
  if ($barcodefile) {
    &run("$NNMAKE $xx $id $chain -barcode_file $barcodefile",(<$xx$id$chain*$TAIL>));
  } else {
    &run("$NNMAKE",(<$xx$id$chain*$TAIL>));
  }
# TODO: add new checkups to make sure the picker succeeded
#  ((-f <$xx$id$chain*$TAIL>) && (-s <$xx$id$chain*$TAIL>)) || die ("nnmake failed!\n");
#   push(@cleanup_files,(<zscore_*_$TAIL_$id$chain>,<status.*$TAIL_$xx$id>,<names.*$TAIL$xx$id$chain>));
}
if ($make_loops) {
    unless (&fileExist("$template.loops" ))  {
	if (!&try_try_again("$NNMAKE $xx $id $chain $template",
			    2, ["$template.loops_all"], ["status.200$TAIL_$xx$id"])) {
	    die("NNMAKE loops failed!\n");
	}
	if (!&try_try_again("$TRIMLOOPS $template.loops_all > $template.loops",
			    2, ["$template.loops"], ["$template.loops_all"])) {
	    die("TrimLoops loops failed!\n");
	}
    }
}

#############################################
# cleanup
#############################################

if ($cleanup) {
  my $file;

  foreach $file (@cleanup_files) {
    unlink($file);
  }
}

# done
exit 0;

###############################################################################
# util
###############################################################################

# getCommandLineOptions()
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;

    my $usage = qq{usage: $0  [-rundir <full path to run directory (default = ./)>
		   \t-verbose  specify for chatty output
		   \t-id  <5 character pdb code/chain id, e.g. 1tum_>
		   \t-nopsipred  don\'t run psipred.
		   \t-psipred_iter  \# of psipred iterations
		   \t-psipred_hbias  psipred helix bias
		   \t-psipred_sbias  psipred strand bias
		   \t-nosam   don\'t run sam
		   \t-nojufo  don\'t run jufo
                   \t-noprof  don\'t run prof
		   \t-nohoms  specify to omit homologs from search
		   \t-psipredfile  <path to file containing psipred ss prediction>
		   \t-samfile  <path to file containing sam ss prediction>
		   \t-phdfile  <path to file containing phd ss prediction>
		   \t-jufofile  <path to file containing jufo ss prediction>
		   \t-minss  <minimum \# ss predictions needed (default=1)>
		   \t-xx  <silly little 2-letter code, if you care>
		   \t-nocleanup  specify to keep all the temporary files
                   \t-template <template base name>  make loop library in addition to fragments]
                   \t-barcode_file <barcode file> restrict bb_big_bin and/or ss
		   \t<fasta file>
		   \n\nVersion: $VERSION\n};

    # Get args
    my %opts = ();
    &GetOptions (\%opts, "jufo!", "prof!", "psipred!", "psipred_iter=f",
                 "psipred_hbias=f", "psipred_sbias:f", "sam!", "homs!", "psipredfile=s",
                 "samfile=s", "phdfile=s", "jufofile=s", "xx=s", "verbose!", "rundir=s",
                 "id=s", "minss=i", "cleanup!", "template=s","barcode_file=s");


    if (scalar(@ARGV) != 1) {
      die "$usage\n";
    }

    $opts{f} = $ARGV[0];

    if ($opts{f} =~ /[^\w\.\/]/) {
      die("Only alphanumeric characters, . and _ are allowed in the filename.\n");
    }

    &checkExist("f",$opts{f});

    return %opts;
}

# checkExist()
#
sub checkExist {
    my ($type, $path) = @_;
    if ($type eq 'd') {
	if (! -d $path) {
            print STDERR "$0: dirnotfound: $path\n";
            exit -3;
	}
    }
    elsif ($type eq 'f') {
	if (! -f $path) {
            print STDERR "$0: filenotfound: $path\n";
            exit -3;
	}
	elsif (! -s $path) {
            print STDERR "$0: emptyfile: $path\n";
            exit -3;
	}
    }
}

sub fileExist {
  my $path = shift;

  return 0 if (! -f $path);
  return 0 if (-z $path);

  return 1;
}

sub run {
  my ($cmd, @files) = @_;
  my $pid;

  if ($DEBUG) {
    print "cmd is: $cmd\n";
  }

 FORK: {
    if ($pid = fork) {
      # parent
      local $SIG{TERM} = sub {
	kill 9, $pid;
	`rm -f @files`;
	exit;
      };

      local $SIG{INT} = sub {
	kill 9, $pid;
	`rm -f @files`;
	exit;
      };

      wait;

      return $?;

    } elsif (defined $pid) {
      #child
      exec($cmd);
    } elsif ($! =~ /No more process/) {
      # recoverable error
      sleep 5;
      redo FORK;

    } else {
      # unrecoverable error
      die ("couldn't fork: $!\n");
    }
  }
}

sub try_try_again {
  my ($cmd, $max_tries, $success_files, $cleanup_files) = @_;
  my $try_count = 0;
  my $missing_files = 1;
  my $f;

  # if at first you don't succeed in running $cmd, try, try again.
  while (($try_count < $max_tries) && ($missing_files > 0)) {
    sleep 10;
    &run($cmd,@$cleanup_files);
    ++$try_count;

    $missing_files = 0;

    foreach $f (@$success_files) {
      if (!&fileExist($f)) {
	++$missing_files;
      }
    }
  }

  # but if you've tried $max_tries times, give up.
  # there's no use being a damn fool about things.
  if ($missing_files > 0) {
    return 0;
  }

  return 1;
}

sub exclude_outn {
  my ($pdb, $chain) = @_;
  my @hits;
  my $hit;
  my $hit_pdb;
  my $hit_chain;
  my %uniq_hits;

  @currentFiles = ("$pdb$chain.homolog_vall");

  open(OUTN, "$pdb$chain.outn");
  open(EXCL, ">$pdb$chain.homolog_vall");

  @hits = grep(/^>/,<OUTN>);

  foreach $hit (@hits) {
    $uniq_hits{$hit} = 1;
  }


  foreach $hit (sort keys %uniq_hits) {
    ($hit_pdb) = $hit =~ /^>\s*(\w+)/;
    $hit_chain = substr($hit_pdb, 4, 1);
    $hit_pdb = substr($hit_pdb, 0, 4);

    $hit_pdb =~ tr/[A-Z]/[a-z]/;
    $hit_chain = '_' if ($hit_chain eq ' ');
    $hit_chain = '_' if ($hit_chain eq '0');

    print EXCL "$pdb$chain $hit_pdb$hit_chain\n";
  }

  close(EXCL);

  @currentFiles = ();
}

sub exclude_blast {
  my ($pdb, $chain) = @_;
  my @hits;
  my $hit;
  my $hit_pdb;
  my $hit_chain;
  my %uniq_hits;

  my $blast_file   = "$pdb$chain.blast";
  my $exclude_file = "$pdb$chain.homolog_nr";
  @currentFiles = ($exclude_file);

  # don't try to do this if the $blast_file doesn't exist
  if ( ! -f $blast_file ) {
    warn "Can't open file $blast_file, so not excluding homologs from it!\n";
    return;
  }
  open(BLAST, $blast_file);
  open(EXCL, ">$exclude_file");

  @hits = grep(/pdb\|/, <BLAST>);

  foreach $hit (@hits) {
    $uniq_hits{$hit} = 1;
  }

  foreach $hit (sort keys %uniq_hits) {
      while ($hit =~ s/pdb\|(\w{4})\|(.?)//) {
	  $hit_pdb   = $1;
	  $hit_chain = $2;

	  $hit_pdb   =~ tr/[A-Z]/[a-z]/;
	  $hit_chain = '_' if ($hit_chain =~ /^\s*$/);

	  print EXCL "$pdb$chain $hit_pdb$hit_chain\n";
      }
  }

  close(EXCL);

  @currentFiles = ();
}

## parse_checkpoint_file -- parses a PSI-BLAST binary checkpoint file.
#
# args:  filename of checkpoint file.
# rets:  N x 20 array containing checkpoint weight values, where N
#        is the size of the protein that BLAST thought it saw...

sub parse_checkpoint_file {
  my $filename = shift;
  my $buf;
  my $seqlen;
  my $seqstr;
  my $i;
  my $j;
  my @aa_order = split(//,'ACDEFGHIKLMNPQRSTVWY');
  my @altschul_mapping = (0,4,3,6,13,7,8,9,11,10,12,2,14,5,1,15,16,19,17,18);
  my @w;
  my @output;

  open(INPUT, $filename) || die ("Couldn't open $filename for reading.\n");

  read(INPUT, $buf, 4) || die ("Couldn't read $filename!\n");
  $seqlen = unpack("i", $buf);

  (!$DEBUG) || print "Sequence length: $seqlen\n";

  read(INPUT, $buf, $seqlen) || die ("Premature end: $filename.\n");
  $seqstr = unpack("a$seqlen", $buf);

  for ($i = 0; $i < $seqlen; ++$i) {
    read(INPUT, $buf, 160) || die("Premature end: $filename, line: $i\n");
    @w = unpack("d20", $buf);

    for ($j = 0; $j < 20; ++$j) {
      $output[$i][$j] = $w[$altschul_mapping[$j]];
    }
  }

  return @output;
}

## finish_checkpoint_matrix -- "finishes" the parsed PSI-BLAST checkpoint matrix,
##                             by adding pseudo-counts to any empty columns.
#
# args:  1) sequence string  2) array returned by parse_checkpoint_file
# rets:  "finished" array.  suitable for printing, framing, etc.

sub finish_checkpoint_matrix {
  my ($s, @matrix) = @_;
  my @sequence = split(//,$s);
  my $i;
  my $j;
  my $line;
  my $sum;
  my @words;
  my @b62;
  my @blos_aa = (0,14,11,2,1,13,3,5,6,7,9,8,10,4,12,15,16,18,19,17);
  my %aaNum = ( A => 0,
		C => 1,
		D => 2,
		E => 3,
		F => 4,
		G => 5,
		H => 6,
		I => 7,
		K => 8,
		L => 9,
		M => 10,
		N => 11,
		P => 12,
		Q => 13,
		R => 14,
		S => 15,
		T => 16,
		V => 17,
		W => 18,
		Y => 19,
		X => 0   ### cheep fix for now ##
	      );

  (length($s) == scalar(@matrix)) || die ("Length mismatch between sequence and checkpoint file!\n");

  open(B62,"$BLOSUM/blosum62.qij") || die "couldnt find blosum62 table\n";
  $i = 0;

  # read/build the blosum matrix
  while (<B62>) {
    next if ($_ !~ /^\d/);
    chomp;
    @words = split(/\s/);

    for ($j = 0; $j <= $#words; ++$j) {
      $b62[$blos_aa[$i]][$blos_aa[$j]] = $words[$j];
      $b62[$blos_aa[$j]][$blos_aa[$i]] = $words[$j];
    }

    ++$i;
  }

  # normalize the blosum matrix so that each row sums to 1.0
  for ($i = 0; $i < 20; ++$i) {
    $sum = 0.0;

    for ($j = 0; $j < 20; ++$j) {
      $sum += $b62[$i][$j];
    }

    for ($j = 0; $j < 20; ++$j) {
      $b62[$i][$j] = ($b62[$i][$j] / $sum);
    }
  }

  # substitute appropriate blosum row for 0 rows
  for ($i = 0; $i <= $#matrix; ++$i) {
    $sum = 0;

    for ($j = 0; $j < 20; ++$j) {
      $sum += $matrix[$i][$j];
    }

    if ($sum == 0) {
      for ($j = 0; $j < 20; ++$j) {
	$matrix[$i][$j] = $b62[$aaNum{$sequence[$i]}][$j];
      }
    }
  }

  return @matrix;
}

sub write_checkpoint_file {
  my ($filename, $sequence, @matrix) = @_;
  my $row;
  my $col;
  my @seq = split(//,$sequence);

  open(OUTPUT, ">$filename");

  die ("Length mismatch between sequence and checkpoint matrix!\n") unless (length($sequence) == scalar(@matrix));

  print OUTPUT scalar(@matrix),"\n";

  for ($row = 0; $row <= $#matrix; ++$row) {
    print OUTPUT "$seq[$row] ";
    for ($col = 0; $col < 20; ++$col) {
      printf OUTPUT "%6.4f ", $matrix[$row][$col];
    }
    print OUTPUT "\n";
  }

  print OUTPUT "END";

  close(OUTPUT);
}

## verify_existing_checkpoint
#
# Does what it says...it verifies a checkpoint file to make sure it's kosher.
#
# arg: filename of checkpoint file.
# rets: 1 if good, 0 if bad.

sub verify_existing_checkpoint {
  my $filename = shift;
  my @file;

  open(CPT, $filename) || return 0;

  (!$DEBUG) || print "opened checkpoint file.\n";

  @file = <CPT>;

  ($file[0] =~ /^(\d+)$/) || return 0;

  (!$DEBUG) || print "starting line ok.\n";

  ((scalar(@file) - 2) == $1) || return 0;

  (!$DEBUG) || print "length ok.\n";

  ($file[$#file] =~ /^END/) || return 0;

  (!$DEBUG) || print "end ok.\n";

  close(CPT);
  return 1;
}

## convert 6 state rdb file to 3 state
sub condense_rdb6_rdb {
	my($rows)=@_;
	my($outrows,@parts,$line);
	for (@{$rows}) {
		if ( /^#/) {
			push (@{$outrows}, $_);
		} elsif (/^Pos/) {
			push (@{$outrows}, "Pos\tAA\tE\tH\tL\n");
		} elsif (/^10N/) {
			push (@{$outrows}, "10N\t1S\t5N\t5N\t5N\n");
		} elsif (/^[0-9]/) {
			@parts = split;
			die "Wrong number of columns in results row: $_\n" if $#parts != 7;
			$line = sprintf "%8s\t%1s\t%8s\t%8s\t%8s\n",$parts[0]-1,$parts[1],$parts[2]+$parts[3],$parts[4]+$parts[5],$parts[6]+$parts[7];
			push (@{$outrows}, $line);
		} else {
			die "Unrecognized row: $_\n";
		}
	}
	return $outrows;
}

