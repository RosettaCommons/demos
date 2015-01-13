#!/usr/bin/perl

if(scalar(@ARGV) != 2)
{
	print "<file with rmsds> <rmsd col. #>\n" and die;
}

my $rmsd_filename = $ARGV[0];
my $rmsd_column = $ARGV[1];

open(SCORE, "$rmsd_filename") or print "No $rmsd_filename" and die;
my @score_lines = <SCORE>;

foreach (my $x = 0; $x < scalar(@score_lines); ++$x)
{
	my @line = split(' ',$score_lines[$x]);
	if(@line[0] ne "tag")
	{
		push (@scores, $line[$rmsd_column]); 
	}
}

#print "@scores\n";
#print "$scores[0]\n";
$n0 = 0;
$n1 = 0;
$n2 = 0;
$n3 = 0;
$n4 = 0;
$n5 = 0;
$n6 = 0;
$n7 = 0;
$n8 = 0;
$n9 = 0;
$n10 = 0;
$n11 = 0;
$n12 = 0;
$n13 = 0;
$n14 = 0;
$n15 = 0;
$n16 = 0;
$n17 = 0;
$n18 = 0;
$n19 = 0;
$n20 = 0;

#Tally number of models of RMSD in bins
foreach $score(@scores)
{
	if (($score > 0.00) && ($score < 0.5))
	{
		++$n0
	} 
	if (($score >= 0.5) && ($score < 1.5))
        {
                ++$n1
        }
	if (($score >= 1.50) && ($score < 2.5))
        {
                ++$n2
        }
        if (($score >= 2.50) && ($score < 3.5))
        {
                ++$n3
        }
        if (($score >= 3.50) && ($score < 4.5))
        {
                ++$n4
        }
        if (($score >= 4.50) && ($score < 5.5))
        {
                ++$n5
        }
        if (($score >= 5.50) && ($score < 6.5))
        {
                ++$n6
        }
        if (($score >= 6.50) && ($score < 7.5))
        {
                ++$n7
        }
        if (($score >= 7.50) && ($score < 8.5))
        {
                ++$n8
        }
        if (($score >= 8.50) && ($score < 9.5))
        {
                ++$n9
        }
	if (($score > 9.50) && ($score < 10.5))
        {
                ++$n10
        }
        if (($score >= 10.50) && ($score < 11.5))
        {
                ++$n11
        }
        if (($score >= 11.50) && ($score < 12.5))
        {
                ++$n12
        }
        if (($score >= 12.50) && ($score < 13.5))
        {
                ++$n13
        }
        if (($score >= 13.50) && ($score < 14.5))
        {
                ++$n14
        }
        if (($score >= 14.50) && ($score < 15.5))
        {
                ++$n15
        }
        if (($score >= 15.50) && ($score < 16.5))
        {
                ++$n16
        }
        if (($score >= 16.50) && ($score < 17.5))
        {
                ++$n17
        }
        if (($score >= 17.50) && ($score < 18.5))
        {
                ++$n18
        }
	if( $score >= 18.5)
	{
		++$n19
		#die "RMSD is $score is way too large somethings wrong\n";
	}
	if( $score <= 8.00)
	{
		++$n20
	}
	if($score <= 4.00)
	{
		++$n21
	}
}
$total_frequency = $n0 + $n1 + $n2 + $n3 + $n4 + $n5 + $n6 + $n7 + $n8 + $n9 + $n10 + $n11 + $n12 + $n13 + $n14 + $n15 + $n16 + $n17 + $n18 + $n19;
$n0_tf = 100 * $n0 / $total_frequency;
$n1_tf = 100 * $n1 / $total_frequency;
$n2_tf = 100 * $n2 / $total_frequency;
$n3_tf = 100 * $n3 / $total_frequency;
$n4_tf = 100 * $n4 / $total_frequency;
$n5_tf = 100 * $n5 / $total_frequency;
$n6_tf = 100 * $n6 / $total_frequency;
$n7_tf = 100 * $n7 / $total_frequency;
$n8_tf = 100 * $n8 / $total_frequency;
$n9_tf = 100 * $n9 / $total_frequency;
$n10_tf = 100 * $n10 / $total_frequency;
$n11_tf = 100 * $n11 / $total_frequency;
$n12_tf = 100 * $n12 / $total_frequency;
$n13_tf = 100 * $n13 / $total_frequency;
$n14_tf = 100 * $n14 / $total_frequency;
$n15_tf = 100 * $n15 / $total_frequency;
$n16_tf = 100 * $n16 / $total_frequency;
$n17_tf = 100 * $n17 / $total_frequency;
$n18_tf = 100 * $n18 / $total_frequency;
$n19_tf = 100 * $n19 / $total_frequency;
$totalunder7_5A = $n0 + $n1 + $n2 + $n3 + $n4 + $n5 + $n6 + $n7;
$percentunder7_5A = 100 * $totalunder7_5A / $total_frequency;
$totalunder3_5A = $n0 + $n1 + $n2 + $n3;
$percentunder3_5A = 100 * $totalunder3_5A / $total_frequency;
$totalunder4A = $n21;
$percentunder4A = 100 * $totalunder4A / $total_frequency;
$totalunder8A = $n20;
$percentunder8A = 100 * $totalunder8A / $total_frequency;

#Report RMSD bin values
print "Range"."\t"."Center"."\t"."Freq"."\t"."Percentage"."\n";
print "(0.00,0.5)"."\t"."0"."\t"."$n0"."\t"."$n0_tf"."\n";
print "[0.5,1.5)"."\t"."1"."\t"."$n1"."\t"."$n1_tf"."\n";
print "[1.50,2.5)"."\t"."2"."\t"."$n2"."\t"."$n2_tf"."\n";
print "[2.50,3.5)"."\t"."3"."\t"."$n3"."\t"."$n3_tf"."\n";
print "[3.50,4.5)"."\t"."4"."\t"."$n4"."\t"."$n4_tf"."\n";
print "[4.50,5.5)"."\t"."5"."\t"."$n5"."\t"."$n5_tf"."\n";
print "[5.50,6.5)"."\t"."6"."\t"."$n6"."\t"."$n6_tf"."\n";
print "[6.50,7.5)"."\t"."7"."\t"."$n7"."\t"."$n7_tf"."\n";
print "[7.50,8.5)"."\t"."8"."\t"."$n8"."\t"."$n8_tf"."\n";
print "[8.50,9.5)"."\t"."9"."\t"."$n9"."\t"."$n9_tf"."\n";
print "[9.50,10.5)"."\t"."10"."\t"."$n10"."\t"."$n10_tf"."\n";
print "[10.50,11.5)"."\t"."11"."\t"."$n11"."\t"."$n11_tf"."\n";
print "[11.50,12.5)"."\t"."12"."\t"."$n12"."\t"."$n12_tf"."\n";
print "[12.50,13.5)"."\t"."13"."\t"."$n13"."\t"."$n13_tf"."\n";
print "[13.50,14.5)"."\t"."14"."\t"."$n14"."\t"."$n14_tf"."\n";
print "[14.50,15.5)"."\t"."15"."\t"."$n15"."\t"."$n15_tf"."\n";
print "[15.50,16.5)"."\t"."16"."\t"."$n16"."\t"."$n16_tf"."\n";
print "[16.50,17.5)"."\t"."17"."\t"."$n17"."\t"."$n17_tf"."\n";
print "[17.5,18.5)"."\t"."18"."\t"."$n18"."\t"."$n18_tf"."\n";
print "[18.5,19.5)"."\t"."19"."\t"."$n19"."\t"."$n19_tf"."\n";
print "Total"."\t"."$total_frequency"."\n";
print "\n";

#Compute average top 10% RMSD and best RMSD
my @ordered_scores = sort {$a <=> $b} (@scores);
my @ordered_scores_fin;
my $no_scores = scalar(@ordered_scores);
for (my $x = 0; $x < $no_scores; $x++)
{
	if ($ordered_scores[$x] > 0.00)
	{
		push (@ordered_scores_fin, $ordered_scores[$x]); 
	}
}
print "The Best model by RMSD: $ordered_scores_fin[0]\n";
my $no_scores_fin = scalar(@ordered_scores_fin);
my $ten_percent = int(0.1 * $no_scores_fin);
my $sum_for_avg = 0;
for (my $y = 0; $y < $ten_percent; $y++)
{
	$sum_for_avg += @ordered_scores_fin[$y];	
}
my $average_top_ten_rmsd = $sum_for_avg / $ten_percent;
print "Average Top Ten Percent RMSD: $average_top_ten_rmsd\n";
print "Percentage of Models under 7.5A: $percentunder7_5A\n";
print "Percentage of Models under 3.5A: $percentunder3_5A\n";
print "Percentage of Models under 8A:  $percentunder8A\n";
