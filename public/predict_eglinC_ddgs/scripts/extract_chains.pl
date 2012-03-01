#!/usr/bin/perl

use strict;
#extract specified chains 

if(scalar(@ARGV) != 2){
    print "usage: $0 <pdb> <chain>\n";
    print "extracts specified chain\n";
    exit(1);
}

my $pdb = shift;
my $chain = shift;

chomp $chain;
open(PDB,"<$pdb")
    or die " cannot open file: $pdb\n";

while(my $line = <PDB>){
    if($line =~ m/^ATOM/){
	if(substr($line,21,1) eq $chain){
	    $line =~ s/ $chain /   /;
	    print $line;
	}elsif($chain eq "_" && substr($line,21,1) eq " " || 
	                        substr($line,21,1) eq "_"){
	    print $line;
	}
    }
}
