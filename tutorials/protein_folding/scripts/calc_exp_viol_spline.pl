#!/usr/bin/perl

$filename = $ARGV[0];
$restraintfile = $ARGV[1];
$num_restr = $ARGV[2];
$num_args = 3;

if(scalar(@ARGV) < $num_args)
{
        print "USAGE:  <pdb_filename> <restraint_filename> <# restraints>\n";
        die;
}

#$filename = $ARGV[0];
#$restraintfile = $ARGV[1];
#$num_restr = $ARGV[2];
$atomA ;
$atomB ;
$aminoacidA ;
$aminoacidB ;
$current_viol = 0;
$sum_current_viol = 0;
$max_current_viol = 0;

$xA = 0.0;
$yA = 0.0;
$zA = 0.0;

$xB = 0.0;
$yB = 0.0;
$zB = 0.0;
@splitted_restraint_line[0] = 0;

open( READ_RESTRAINTS, $restraintfile) or die "unable to open".&restraintfile;

# skip over beginnig of file to get to the actual restraints
#while(@splitted_restraint_line[0] != $num_restr)
#{
#        @splitted_restraint_line = split /\s+/, <READ_RESTRAINTS>;
#}

# read in the restraints from the restraint file
while( !eof(READ_RESTRAINTS))
{
	@splitted_restraint_line = split /\s+/, <READ_RESTRAINTS>;
	$atomA = @splitted_restraint_line[1];
	$atomB = @splitted_restraint_line[3];
	$aminoacidA = @splitted_restraint_line[2];
	$aminoacidB = @splitted_restraint_line[4];
	$restraint_max = @splitted_restraint_line[7];
#	$restraint_min = @splitted_restraint_line[6];
	print "residueA: ".$aminoacidA."\t";
	print $atomA."\t";
	print "residueB: ".$aminoacidB."\t";
        print $atomB."\t";
	print "exp_val:\t$restraint_max\t"; 
	
#	print @splitted_restraint_line[0]." "."@splitted_restraint_line[1]"." ".@splitted_restraint_line[2]." "."@splitted_restraint_line[3]"." "."$atomA"." "."$atomB"." ".$aminoacidA." ".$aminoacidB."\n";
	
	open( READ, $filename) or die "unable to open".&filename;
	# read through the pdb file
	while( !eof( READ))
	{
		@splitted_line = split /\s+/, <READ>;
#		print "splitted line: ".@splitted_line[0]." ".@splitted_line[2]." ".@splitted_line[4]."\n";
		if( @splitted_line[0] eq "ATOM" && @splitted_line[2] eq $atomA && @splitted_line[5] == $aminoacidA)
		{
#			print "hello\n";
			$xA = @splitted_line[6];
			$yA = @splitted_line[7];
			$zA = @splitted_line[8];
#			print "xA: ".$xA."\n";
#			print "yA: ".$yA."\n";
#			print "zA: ".$zA."\n";
		}
	
	        if( @splitted_line[0] eq "ATOM" && @splitted_line[2] eq $atomB && @splitted_line[5] == $aminoacidB)
       		 {
                	$xB = @splitted_line[6];
                	$yB = @splitted_line[7];
	                $zB = @splitted_line[8];
#			print "xB: ".$xB."\n";
#			print "yB: ".$yB."\n";
#			print "zB: ".$zB."\n";
	        }

	}
	
	$distance = sqrt(($xA-$xB)*($xA-$xB) + ($yA-$yB)*($yA-$yB) + ($zA-$zB)*($zA-$zB));
	$rounded_distance = sprintf("%.1f",$distance);
#	print " distance ".$distance;	
	if( $distance < $restraint_max)
	{
		$current_viol =  $restraint_max - $distance;
		$sum_current_viol = $sum_current_viol + $current_viol;
	}
	
	if( $distance > $restraint_max)
	{
		$current_viol =  $distance - $restraint_max;
		$sum_current_viol = $sum_current_viol + $current_viol;
	}
	if( $current_viol > $max_current_viol)
	{
		$max_current_viol = $current_viol;
	}
	
#	print $aminoacidA;
#	print " ".$atomA;
#	print " ".$aminoacidB;
#	print " ".$atomB;
	print "Distance: ".$distance."\t";
#	print " Restraint_max: ".$restraint_max."\t";
#	print " Restraint_min: ".$restraint_min."\t";
	print " Current_viol: ".$current_viol."\n";
#	print " sum_current_viol: ".$sum_current_viol."\n";
	$current_viol = 0;	
}

print "Filename: ".$filename."\t";
print " Restraintfile: ".$restraintfile."\t";
print " sum_current_viol: ".$sum_current_viol."\t";
print " max_current_viol: ".$max_current_viol."\n";
