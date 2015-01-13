#!/usr/bin/perl -w

# contactOrder.pl
# computes contact order of (protein) structures from PDB files
# Author: Eric Alm <ealm3141@users.sourceforge.net>

use strict;

my $helpStr = qq{

  Instruction:

    contactOrder.pl <options> <pdbfile>

  Options:

    -c: set contact cutoff in Angstroms, default = 6
    -r: returns relative CO
    -a: returns absolute CO

  Description:

    contactOrder.pl takes a PDB file and returns the contact order
    absolute = average sequence separation of contacting atoms
    relative = average sequence separation / protein length

  Warnings:

    PDB files that are not numbered sequentially or that contain
    long disordered regions (at the N- or C-termini) will give
    inaccurate results

};

#read command line options
#======================================================================
my ($absolute,$relative,$cutoff,$last_res);

my %opts = &getCommandLineOptions ();
$absolute = (exists $opts{'a'}) ? 1 : 0;
$relative = (exists $opts{'r'}) ? 1 : 0;
$relative = 1 if (!($absolute || $relative));
$cutoff = $opts{'c'} if exists $opts{'c'};
$cutoff ||= 6.0;
$cutoff *= $cutoff; #actually uses squared distance cutoff

unless (@ARGV == 1){die "$helpStr";}


#main function
{
    my %atm;
    my ($last_id,$chain_id,$max_res,$min_res,$non_seq_flag,$res,$last_res);
    my @atom_list;
    my $first_pass = 1;
		my $pdbf = $ARGV[0];
		chomp $pdbf;
		open(PDBFILE, $pdbf);

    #read in PDB
    #======================================================================
    while(<PDBFILE>){

	last if (/^TER/);
	if(/^ATOM  /){
	    %atm = &parseAtomLine($_);
	    next if $atm{atom_name} =~ /H/;
	    $last_id = $chain_id;
	    $chain_id = $atm{chain_id};
	    if(!$first_pass){
		last if ($chain_id ne $last_id);
	    }else{
		$max_res = $atm{res_num};
		$min_res = $atm{res_num};
	    }
	    push @atom_list, {%atm};
	    $last_res = $res;
	    $res = $atm{res_num};
	    if(!$first_pass){
		$non_seq_flag = 1 if (abs($res-$last_res)>1);
	    }
	    if($atm{res_num}>$max_res){
		$max_res = $atm{res_num};
	    }
	    if($atm{res_num}<$min_res){
		$min_res = $atm{res_num};
	    }
	    $first_pass = 0;
	}#if matches ATOM

    }#while PDBFILE


    #compute CO
    #======================================================================
    my $counts;
    my $order;
    for my $atom1 (@atom_list){
	for my $atom2 (@atom_list){
	    my $seq_dist = $atom1->{res_num} - $atom2->{res_num};
	    if($seq_dist > 0){
		if(&withinDist($atom1,$atom2,$cutoff)){
		    $counts++;
			$order += $seq_dist;
		}
	    }
	}
    }

    #output results
    #======================================================================
    print "Absolute Contact Order: ",$order/$counts,"\n" if $absolute;
    print "Relative Contact Order: ",$order/$counts/($max_res-$min_res+1),"\n" if $relative;
    print "Warning: nonsequential numbering in PDB!" if $non_seq_flag;

}#end program


#do atoms contact?
#======================================================================
sub withinDist{
    my ($atm1,$atm2,$cutoff) = @_;
    my $sqr_dist =
	($atm1->{x}-$atm2->{x})*($atm1->{x}-$atm2->{x}) +
	($atm1->{y}-$atm2->{y})*($atm1->{y}-$atm2->{y}) +
	($atm1->{z}-$atm2->{z})*($atm1->{z}-$atm2->{z});
    return $sqr_dist < $cutoff;
}


#read single atom from PDB
#======================================================================
sub parseAtomLine{
    my $line = shift;
    my %atom;

    $atom{atom_num} = substr($line,6,5);
    $atom{atom_name} = substr($line,12,4);
    $atom{res_type} = substr($line,17,3);
    $atom{res_num} = substr($line,22,4);
    $atom{x} = substr($line,30,8);
    $atom{y} = substr($line,38,8);
    $atom{z} = substr($line,46,8);
    $atom{chain_id} = substr($line,21,1);

    return %atom;
}


sub getCommandLineOptions {
    use Getopt::Long;

    # Get args
    #
    my %opts = ();
    &GetOptions (\%opts, "r", "c=i", "a");

    return %opts;
}



