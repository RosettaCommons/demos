package hhsearch::Util;

@ISA = qw/ Exporter /;
@EXPORT = qw/
	make_hhm
	psiblast_to_hhm
	run_NN_filter
/;

use Rosetta::Util;
use Rosetta::Job;
use alignment;

# configuration options
$cleanup_tempfiles = 1;

sub make_hhm {
	my $fasta_fn = shift;
	my $hhm_dir  = shift;
	my $options  = shift;

	# brief hack for hhsearch scripts: all of these assume that the
	# > line in the .fasta file is not blank.
	my $id = id_from_fasta($fasta_fn);
	if ( $id eq '' || $id =~ /^\s*/ ) {
		use Rosetta::Sequence qw/ read_protein_fasta_fn write_fasta_file /;
		my $seq = read_protein_fasta_fn($fasta_fn);
		unlink $fasta_fn;
		write_fasta_file($fasta_fn,$seq,$fasta_fn);
	}

	if ( $options->{use_hhblits} ) {
		return make_hhm_hhblits($fasta_fn,$hhm_dir,$options);
	} else {
		return make_hhm_psiblast($fasta_fn,$hhm_dir,$options);
	}
}

sub id_from_fasta {
	my $fasta = shift;

	open FILE, "<$fasta" or die $!;
	my $header = '';
	while ( $header !~ /^>/ ) {
		$header = <FILE>;
	}
	close FILE or die $!;
	my $id = substr($line,1);

	return $id;
}

sub make_hhm_psiblast {
	my $fasta_fn = shift;
	my $hhm_dir  = shift;
	my $options  = shift;
	use File::Path;
	use File::Copy;
	use File::Basename;

	my ($id,undef) = basename($fasta_fn,'.fasta');
	my $local_hhm_dir = $hhm_dir;
	mkpath $local_hhm_dir;

	# run psi-blast
	if ( ! -f "$local_hhm_dir/$fasta_fn.tar.gz" ) {
		my $cmd = join ' ', (
			"$options->{script_base}/run-psiblast.pl",
			"--n_procs  $options->{n_procs}",
			"--n_rounds $options->{n_rounds}",
			"$local_hhm_dir/$fasta_fn --archive",
		);
		if ( exists $options->{nr_db} ) {
			$cmd = join ' ', ( $cmd, "--db", $options->{nr_db} );
		}
		exec_cmd($cmd);
	}

	# extract blast output from the .tar.gz file
	my $tar_file = "$id.fasta.tar.gz";
	my $psiblast_output = join '.', ( $id,"fasta",$options->{n_rounds},"psiblast" );
	my $cmd = "cd $local_hhm_dir; tar xzvf $tar_file $psiblast_output";
	exec_cmd($cmd);
	if ( ! -f "$local_hhm_dir/$psiblast_output" ) {
		die "Error: couldn't extract file $psiblast_output from file $tar_file!\n";
	}

	return psiblast_to_hhm( $id, "$local_hhm_dir/$psiblast_output", $options );
}

sub make_hhm_hhblits {
	my $fasta_fn = shift;
	my $hhm_dir  = shift;
	my $options  = shift;

	use Cwd;
	use File::Basename;

	my $numb_rounds = 3;
	my $mact = 0.5;

	my ($id,undef) = basename($fasta_fn,'.fasta');
	my $hhblits_path = '/work/robetta/src/rosetta_server/bin/hhblits/bin_64bit_static';
	my $local_hhm_dir = dirname  $hhm_dir;
	my $hhm_file = "$hhm_dir/$id.hhm";
	my $a3m_file = "$hhm_dir/$id.a3m";
	my $hhblits_out = "$hhm_dir/$id.hhblits";
	my @hhblits_jobs;
	my $job = Rosetta::Job->new(
		executable   => "$hhblits_path/hhblits",
		results_file => $hhm_file,
		args         => [
			"-i $fasta_fn",
			"-o $hhblits_out",
			"-ohhm $hhm_file",
			"-oa3m $a3m_file",
			"-n $numb_rounds",
			"-nopred",
			"-mact $mact",
			"-addss",
			"-cpu 4",
			"-realign"],
		dir          => '.',
		lockfile     => "$hhm_file.lock",
	);
	push @hhblits_jobs, $job;
	Rosetta::Job::run_jobs_in_parallel( \@hhblits_jobs, 4 );
	if ( ! -f "$id.hhm" ) {
		my $cmd = "$options->{hhsearch_path}/hhmake -i $a3mfile";
		exec_cmd($cmd);
	}
	return "$local_hhm_dir/$id.hhm";


}

sub exec_cmd {
	my $cmd = shift;
	#print "executing cmd: $cmd\n";
	system($cmd);
}

sub psiblast_to_hhm {
	my $id              = shift;
	my $psiblast_output = shift;
	my $options         = shift;

	use Cwd;
	use File::Basename;

	my $local_hhm_dir = dirname  $psiblast_output;
	$psiblast_output  = basename $psiblast_output;

	# Parse out a multiple alignment from the BLAST results.
	my $orig_dir = getcwd;
	chdir $local_hhm_dir;

	my $fasta_fn = "$id.fasta";
	my $a2mfile  = "$id.a2m";
	#print "running alignblast!\n";
	if ( ! -f "$local_hhm_dir/$a2mfile" ) {
		my $cmd = "cd $local_hhm_dir; $options->{hhsearch_path}/alignblast.pl -Q $fasta_fn $psiblast_output $a2mfile";
		exec_cmd($cmd);
	}

	# add psipred secondary structure information
	my $a3mfile = "$id.a3m";
	my $cmd = "$options->{hhsearch_path}/addpsipred.pl $a2mfile";
	exec_cmd($cmd);
	#print "adding psipred!\n";

	# make hhms
	#print "making hhms!\n";
	if ( ! -f "$id.hhm" ) {
		my $cmd = "$options->{hhsearch_path}/hhmake -i $a3mfile";
		exec_cmd($cmd);
	}

	# cleanup temporary files
	my @suffices = qw/
		chk sq psi pn sn in.a3m chk mtx mn blalog ss2 ss horiz a2m a3m blast aux
	/;
	if (0){
	#if ( $cleanup_tempfiles ) {
		foreach my $suffix (@suffices) {
			my $file = "$id.$suffix";
			unlink( $file );
		}
	}
	chdir $orig_dir;
	return "$local_hhm_dir/$id.hhm";
}

sub run_NN_filter {
    my $aln_dir = shift;
    my $NN_program = shift;
    my $targetid = shift;

    my @NN_jobs;
    my @hhr_files = glob("$aln_dir/*.hhr");	
    #-----Calc the NN ranking for each alignment file.
    foreach my $hhr_file(@hhr_files){
	my @path_split = split/\//, $hhr_file;
	my $hhr_file_no_path =$path_split[length($file_names)-1];
	my $hhr_file_no_ext = substr($hhr_file_no_path, 0,(length($hhr_file_no_path ) - 4));
	my $NN_dir = assemble_path( $aln_dir, "NN_dir_$hhr_file_no_ext");
	mkdir_safe($NN_dir);
	my $tab_file_no_path  = "$hhr_file_no_ext.start.tab";
	my $a3m_file_no_path = "$targetid.a3m"; 
	copy_safe("$aln_dir/$hhr_file_no_path","$NN_dir");
	copy_safe("$aln_dir/$tab_file_no_path","$NN_dir");
	system("cp $aln_dir/$a3m_file_no_path $NN_dir/$hhr_file_no_ext.a3m");
	system("cp $aln_dir/$hhr_file_no_path $NN_dir/$hhr_file_no_ext.start.hhr");
	my $results = "$hhr_file_no_ext.NN_results";
	my $job = Rosetta::Job->new(
				    executable => $NN_program,
				    results_file => "$results",
				    args       => [
						   "$NN_dir/$hhr_file_no_ext",
						   "> $results",
						   ],
				    dir => "$NN_dir",
				    lockfile     => "$results.lock",);
	push(@NN_jobs, $job);
    }
    Rosetta::Job::run_jobs_in_parallel( \@NN_jobs, 4);

    #----Move the result hhr files up 1 directory.
     foreach my $hhr_file(@hhr_files){
	my @path_split = split/\//, $hhr_file;
	my $hhr_file_no_path =$path_split[length($file_names)-1];
	my $hhr_file_no_ext = substr($hhr_file_no_path, 0,(length($hhr_file_no_path ) - 4));
	my $NN_dir = assemble_path( $aln_dir, "NN_dir_$hhr_file_no_ext");
	print "copying $NN_dir/$hhr_file_no_path\n";
	copy_safe("$NN_dir/$hhr_file_no_path","$aln_dir");
    }
}

return 1;
