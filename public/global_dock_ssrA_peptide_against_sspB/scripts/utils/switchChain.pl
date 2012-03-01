#!/usr/bin/perl
use FileHandle;
# Changes the name of a chain in the file

if(scalar @ARGV < 2){
	print "Usage: $0 <pdb-file> <chain1> <chain2> <part>\n\n";
	print "Changes the name of a chain in the file.\n";
	print "<part> can be 1(first), 2(second) or 3(both=default).\n";
	print "If chain2 is not defined, will switch space for chain1 in both chain parts.\n";
	die "\n";
}

#ATOM     91  CG2 VAL B  60      -2.997
$file = shift;
$chain1 = shift;
$chain2 = shift;
$takePart = shift; # 1 or 2 or 3(both)

if(! defined $chain2){
# switching from space to chain1
	$chain2 = $chain1;	
	$chain1 = " ";
	$takePart = 3; # 3 = both
}

if(! defined $takePart){
	$takePart = 3;
}

my $fh = new FileHandle;
if (!$fh->open("< $file")) {
        die "Could not open $file\n";
}

%newFile;
$i = 0;
$part = 1;

while(! $fh->eof){
	chomp( $line = $fh->getline() );
	if($line =~ m/^TER/){
		$part = 2;
	}
	if(($takePart == 3 || $part == $takePart) && $line =~ m/^ATOM.................(.).*/){
		if($1 eq $chain1){
			$line =~ s/^(ATOM.................)./$1$chain2/;
		}
	}
	$newFile{$i} = $line;
	$i++;
}

$fh->close;

# Create final file
my $fh = new FileHandle;
if (!$fh->open("> $file")) {
	die "Could not open $file\n";
}

for($j = 0; $j < $i; $j++){
	print $fh $newFile{$j}."\n";
}
$fh->close;
