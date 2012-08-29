#!/usr/bin/perl

use Getopt::Long;

&Getopt::Long::Configure('pass_through', 'no_autoabbrev');
&Getopt::Long::GetOptions('help|h' => \$help,
	'input|i=s' => \$fileName,
	'verbose|v' => \$verbose,
	'score|s:f' => \$scoreCut,
	'ddg|d:f' => \$ddgCut,
	'rmsd|r:f' => \$rmsdCut,
	'sasa:f' => \$sasaCut,
	'unsat|u:f' => \$unsatCut,
	'ddgSasa:f' => \$ddgSasaCut,
	'sasaUnsat:f' => \$sasaUnsatCut,
	'hbondE|b:f' => \$hbondECut,
	'packstat|p:f' => \$packstatCut,
	'scvalue:f' => \$scValueCut,
	);

if($help == 1){
	print "retrieve score file entries on any of the following metrics (assuming they exist in your score file): score file name (REQUIRED)
  -i
  --input

  -v
  --verbose

  score cutoff
  -s
  --score
  
  ddg cutoff
  -d
  --ddg
  
  RMSD cutoff sasa
  -r
  --rmsd

  SASA cutoff
  --sasa

  unsatisfied hydrogen bonds cutoff
  -u
  --unsat

  packstat cutoff
  -p
  --packstat

  DDG/Sasa cutoff
  --ddgSasa

  hbond energy ratio cutoff
  --hbondE

  shape complementarity cutoff
  --scvalue
  
  SASA/# Unsat cutoff
  --sasaUnsat\n";
	exit(0);
}

if(! $fileName){
	print "Error: must provide score file name using -i or --input\n";
	exit(1);
}

chomp $fileName; ## Chomp the last parameter because it has the newline character

open(my $inputFile, "$fileName") || die("Can't open $fileName:!\n");

my @colNames;
my $count = 0;
if(! $verbose){
	print "desc\tenergy\trmsd\tddg\tsasa\tp_sasa\th_sasa\tunsats\thbond_energy\tsc_HbondE\tburied_HbondE\tpackstat\tscvalue\t(dg/sasa)\t(sasa/unsat)\n";
}
while($line=<$inputFile>){
	chomp($line);
	if($line=~/^SCORE/){
		if($count==0){
			@colNames = split(/\s+/,$line);
		}
		else{
			@vals = split(/\s+/,$line);
			for($i=0; $i<=$#colNames; $i++){
				if($colNames[$i] eq "total_score"){
					$energy = $vals[$i];	
				}
				if($colNames[$i] eq "packstat"){
					$packstat = $vals[$i];	
				}
				if($colNames[$i] eq "dG_separated"){
					$ddg = $vals[$i];
				}
				if(grep $_ eq "rms", @colNames){
				  if($colNames[$i] eq "rms"){
				    $rmsd = $vals[$i];
				  }
				}
				else{
				  $rmsd ="N/A";
				}
				if($colNames[$i] eq "dSASA_int"){
					$sasa = $vals[$i];
				}
				if($colNames[$i] eq "dSASA_polar"){
					$pol_sasa = $vals[$i];
				}
				if($colNames[$i] eq "dSASA_hphobic"){
					$h_sasa = $vals[$i];
				}
				if($colNames[$i] eq "delta_unsatHbonds"){
					$unsats = $vals[$i];
				}
				if($colNames[$i] eq "description"){
					$desc = $vals[$i];
				}
				if($colNames[$i] eq "hbond_E_fraction"){
					$hbondEnergy = $vals[$i];
				}
				if($colNames[$i] eq "sc_value"){
					$scValue = $vals[$i];
				}
				if($colNames[$i] eq "buried_Hbond_E"){
					$burHbondValue = $vals[$i];
				}
				if($colNames[$i] eq "int_sc_Hbond_E"){
					$intSCHbondValue = $vals[$i];
				}
			}
			if($sasa != 0){
				$ddgSasa = ($ddg/$sasa);
			}
			if($unsats != 0){
				$sasaUnsat = ($sasa/$unsats);
			}
			else{
				$sasaUnsat = "N/A";
			}
			if( (!$scoreCut or $energy <= $scoreCut) and 
					(!$ddgCut or $ddg <= $ddgCut) and 
					(!$scValueCut or $scValue >= $scValueCut) and 
					(!$rmsdCut or $rmsd<=$rmsdCut) and 
					(!$unsatCut or $unsats<=$unsatCut) and 
					(!$sasaCut or $sasa>=$sasaCut) and 
					(!$packstatCut or $packstat>=$packstatCut) and 
					(!$hbondECut or $hbondEnergy<=$hbondECut) and 
					(!$ddgSasaCut or $ddgSasa <= $ddgSasaCut) and 
					($desc ne "") and 
					(!$sasaUnsatCut or $sasaUnsat eq "N/A" or $sasaUnsat >= $sasaUnsatCut)
					){
					if($verbose){
						print $desc . ":\n\tenergy:" . $energy .
							"\n\trmsd:" . $rmsd . 
							"\n\tddg:" . $ddg . 
							"\n\tsasa:" . sprintf("%.0f", $sasa) . 
							"\n\tsasa:" . sprintf("%.0f", $pol_sasa) . 
							"\n\tsasa:" . sprintf("%.0f", $h_sasa) . 
							"\n\tburied unsatisfied:" . $unsats . 
							"\n\tpackstat:" . $packstat . 
							"\n\tscvalue:" . $scValue . 
							"\n\thbond energy ratio:" . $hbondEnergy . 
							"\n\t(dg/sasa):" . sprintf("%.3f", $ddgSasa) . 
							"\n\t(sasa/unsat):" . sprintf("%.2f", $sasaUnsat) . "\n";
					}
					else{
						#$desc =~ s/.*_//g;
						print "$desc\t$energy\t$rmsd\t$ddg\t" . sprintf("%.0f", $sasa) . "\t" . sprintf("%.0f", $pol_sasa) . "\t" . sprintf("%.0f", $h_sasa) . "\t" . sprintf("%.0f", $unsats) ." \t $hbondEnergy \t $intSCHbondValue \t $burHbondValue \t $packstat \t $scValue \t" 
							. sprintf("%.4f", $ddgSasa) . " \t " . sprintf("%.1f", $sasaUnsat) . "\n";
					}
			}
		}
		$count++;
	}
}
