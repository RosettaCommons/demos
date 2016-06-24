#!/usr/bin/perl -w	

# this script will take a *.jufo9d output file and create the span file necessary to run membrane mode in rosetta
#USAGE: perl jufo9d_span.pl *.jufo9d

$jufo9d_file = $ARGV[0];
chomp($jufo9d_file);

# counter for length of MC
$i = 0;

# counter for number of aa
$total_aa = 0;

# counter for number of MC
$number_of_MC = 0;
$total_AA_MC = 0;

# arrays for start and end positions
@start = ();
@length = ();


open (INPUT, $jufo9d_file) || die "cannot open sequence file";
while ($line = <INPUT>) {
	$total_aa++;
	@aa_line = split(' ', $line);
	$position = $aa_line[0];
	#$aa = $aa_line[1];
	#$ss = $aa_line[2];
	$environment = $aa_line[3];	
	
	# check to see that you have a MC (membrane core) AA
	if ($environment =~ /MC/) {
		# keep track of start of MC
		if ($i == 0) { 
			push(@start, $position);
		}	
		$i++;
	}
	
	else {
		next unless $i != 0;
		#print "length is $i\n";
		# keep track of length of MC
		push(@length, $i);
		# count the number of AA in MCs
		$total_AA_MC = $total_AA_MC + $i;
		# count the membrane core spans
		$number_of_MC++;
		# reset number in MC
		$i = 0;
	}
	
	# print out file in correct format
}

#print "-----------------------------------------\n";
print "TM region prediction for $jufo9d_file predicted using Jufo9D\n";
print "$number_of_MC $total_aa\n";
print "antiparallel\n";
print "n2c\n";
for ($a = 0; $a <= $number_of_MC-1; $a++) {
	# new numbers, renumbered so that beginning of TM span (pos 573) is changed to pos 1
	#$renum_start = $start[$a] - 572;
	#$end = $start[$a] + $length[$a] - 572;
	$end = $start[$a] + $length[$a];
	print " $start[$a]   $end   $start[$a]   $end\n";
}
