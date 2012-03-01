#!/usr/bin/perl

#

my $in_file = shift; #get path to input mutation file
my $offset = shift; #This is the offset for the number of residues that are missing in the PDB file, for example if the first 8 residues are missing in the pdb file, we have to subtract the mutation locations to match the renumbered pdb

my $output_path = "";
if( scalar( @ARGV ) > 0 ) {
  $output_path = shift;
}

open IN, "<", $in_file;

my $total = 0;

my $output = '';
my $mut_set ='';
my $num_mutations = 0;

while ($line = <IN>)		
{

if( $line =~ m/^(\D)(\d+)(\D)/g) { #get each input line and concatenate
	$mut_set = $1 . " " . ($2 - $offset) . " " . $3 . "\n".$mut_set;
	$total++;
	$num_mutations++;
}

elsif($line eq "\n") { #line break indicate the end of a mutation set, 
	$output = $output . "$num_mutations\n" . "$mut_set\n";
	$num_mutations=0;
	$mut_set='';
}


}

$output = $output."$num_mutations\n". "$mut_set\n"; #get last set of mutations

if ($output_path ne "") #print to outfile if specified
{

	open OUT, ">", $output_path;
	print OUT "total ", $total, "\n"; 
	print OUT $output;

}

else #print to standard output if no path
{
	print "total ", $total, "\n"; 
	print $output;
}

