#!/usr/bin/perl

if (!scalar@ARGV) {
        print "USAGE: $0 <alignment file> [format: 'hhsearch' or 'grishin' default: hhsearch]\n";
        exit(1);
}

use FindBin;
my $format = "hhsearch";
my $align = shift@ARGV;
if (scalar@ARGV) {
        $format = shift@ARGV;
}
chomp $align;
chomp $format;
($format =~ /^hhsearch|grishin$/) or die "ERROR! incorrect format specified\n";
open(F, $align) or die "ERROR! cannot open $align: $!\n";
while (<F>) {
        my $code;
        my $chain;
        if ($format eq "hhsearch" && /^>(\w{4})\_(\w)/) {
                $code = $1;
                $chain = $2;
        } elsif ($format eq "grishin" && /^##\s*\S+\s+(\w{4})(\w)\_/) {
                $code = $1;
                $chain = $2;
        }
        if (length($code) && length($chain) && !-s "$code$chain.pdb") {
                my $shell = "$FindBin::Bin/get_pdb.py $code $chain";
                print "$shell\n";
                (!system($shell)) or warn "WARNING! $shell command failed\n";
        }
}
close(F);
