#!/usr/bin/perl

$filename = $ARGV[0];
$restraintfile = $ARGV[1];
$num_restraints = $ARGV[2];
$num_args = 3;

if(scalar(@ARGV) < $num_args)
{
	print "USAGE:  <pdb_filename> <restraint_filename> <# restraints>\n";
	die;
}


#$num_restr = $ARGV[2];
$atomA ;
$atomB ;
$aminoacidA ;
$aminoacidB ;
$current_viol = 0;
$pc = 0;
$pc_viol = 0;
$num_viol = 0;

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
	$restraint_min = @splitted_restraint_line[6];
#	print $atomA." ".@splitted_restraint_line[1]."\n";
#	print $atomB." ".@splitted_restraint_line[3]."\n";
#	print $aminoacidA." ".@splitted_restraint_line[2]."\n";
#	print @splitted_restraint_line[4]." ".$aminoacidB."\n"; 
	
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
		}
	
	        if( @splitted_line[0] eq "ATOM" && @splitted_line[2] eq $atomB && @splitted_line[5] == $aminoacidB)
       		 {
#			print "hello again\n";
                	$xB = @splitted_line[6];
                	$yB = @splitted_line[7];
	                $zB = @splitted_line[8];
	        }

	}
	
	$distance = sqrt(($xA-$xB)*($xA-$xB) + ($yA-$yB)*($yA-$yB) + ($zA-$zB)*($zA-$zB));
	$rounded_distance = sprintf("%.3f",$distance);
#	$rounded_current_viol = sprintf("%.3f", $current_viol);
	$rounded_pc = sprintf("%.3f",$pc);
	$rounded_pc_viol = sprintf("%.3f",$pc_viol);
#	print " distance ".$distance;	
	if( $distance < $restraint_min)
	{
		$current_viol =  $restraint_min - $distance;
		$pc = 	$pc + $current_viol;
		$num_viol++;
	}
	
	if( $distance > $restraint_max)
	{
		$current_viol =  $distance - $restraint_max;
		$pc =   $pc + $current_viol;
		$num_viol++;
	}
	if( $current_viol > $pc_viol)
	{
		$pc_viol = $current_viol;
	}
	$rounded_current_viol = sprintf("%.3f", $current_viol);
	
	print "resA: ".$aminoacidA;
	print " atomA: ".$atomA;
	print " resB: ".$aminoacidB;
	print " atomB: ".$atomB;
	print " distance: ".$rounded_distance;
	print " ub: ".$restraint_max;
	print " lb: ".$restraint_min;
	print " violation: ".$rounded_current_viol."\n";
#	print " ".$rounded_pc."\n";
	$current_viol = 0;	
}

print $filename;
#print " ".$restraintfile;
print " number_violations: ".$num_viol;
print " sum_violations: ".$rounded_pc;
print " max_violation: ".$rounded_pc_viol."\n";
