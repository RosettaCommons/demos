#!/usr/bin/perl
# extracts chains from pdb
#extract_chains.pl
# input: -p pdbin:  pdbfile
#        -c chains: string of chains to extract (orderdependent)
#        -o pdbout: output file
#        -n flag to exclude nonAA res: currently PCA, ACE

#      selects lines of following format:
#ATOM   2318  CA  VAL A 308 
#atom   number CA aa  chain position
# (if chain is A)
#new: AUTOMATICALLY converts MSE to MET
#     use -m Not to do so!!
# -t option: TER  added between chains 
#
# usage: extract_chains.pl -p <pdbin> -c <chains> [-o <pdbout> -t -h -m -n -r]
#
#note that if line contains two copies, the first is used (actually A/1 is used)
#ATOM    625  CB AARG A  72      20.109  13.758  14.197  0.50  8.07           C
#ATOM    625  CB BARG A  72      20.109  13.758  14.197  0.50  8.07           C
#
my ($opts,$args,$err) = &my_getopts('pcor:thnm');  # before : wants arg, after switch (true or false)
@ARGV=@{$args}; my %FLAGS = %{$opts};

$help = $opts->{opt_h};
die "Usage: $0 -p pdbfile -c chains [ -o output (default out) -t (flag for incl of terminus between chains) -r range x_y -n exclude PCA,ACE -m convert mse to met -h (flag for help)] " if ($help);

$pdbin  = $FLAGS{opt_p};
die "no pdb input is given (-p pdb)" if ($pdbin eq "");
$chains = $FLAGS{opt_c};
die "no chains given (-c chains)" if ($chains eq "");
$pdbout = $FLAGS{opt_o};
$pdbout = "out" if ($pdbout eq "");

$r = $FLAGS{opt_r};
if ($r){
  @range = split("-",$r);
print "extract range between residues $range[0] and $range[1]\n";
}else{
  @range = (0,10000);
}

$add_terminus_between_chains =  $opts->{opt_t};
print "TER added between chains\n" if ($add_terminus_between_chains);

$no_PCA =  $opts->{opt_n};

$convert_mse =  $opts->{opt_m};

open (PDBIN,"$pdbin") || die ("file $pdbin not found");
@pdbin = <PDBIN>; close(PDBIN); chomp(@pdbin);

@chainarray = split("",$chains);
for ($c=0;$c<@chainarray;$c++){
  $chainarray[$c] = " " if ($chainarray[$c] eq "0");
  $take_chains{$chainarray[$c]} = $c+1;
}

#print "chainarray @chainarray\n";

$adaptpos = 0;# set to T if position array too large: will lead to entry without space between chain and pos

@take_line = ();
for ($i=0;$i<@pdbin;$i++){
  $line = $pdbin[$i];
  $take_line[$i] = 0;
  if ($line =~ /^MODEL/){
    last if ($model>0);
    $model++;
  }
#convert MSE
  if (!$convert_mse){
    if ($line =~ /^HETATM\s+\d+\s+(\w+)\s+MSE\s.{1}\s*(\d+).{1}\s+.*$/){
      $atom= $1;
      $pos = $2;
      print "convert MSE at position $pos to MET\n" if ($atom eq "CA");
      $line =~  s/HETATM/ATOM  /g;
      $line =~  s/MSE/MET/g;
      $line =~  s/SE / SD/g;
      $pdbin[$i] = $line;
    }
  }
  
  #allow only ATOM correct format lines (or the first line in case of several 
  #ocupencies
  if (!(($line =~ /^ATOM\s+\d+\s+\w+\s+([A-Z]{3})\s(.{1})\s*(\d+)(.{1})\s+.*$/) ||
       ($line =~ /^ATOM\s+\d+\s+\w+\s*[A,1]([A-Z]{3})\s(.{1})\s*(\d+)(.{1})\s+.*$/))) {
    #print lines that are dumped
    if ($line =~ /^ATOM/) {
      print "This line was deleted:\n$line\n";
    }
    next;
  }
  if ($line =~ /^ATOM\s+\d+\s+\w+\s*A([A-Z]{3})\s(.{1})\s*(\d+)(.{1})\s+.*$/) {
    print "For position $pos used conformation A\n";
  }
  elsif ($line =~ /^ATOM\s+\d+\s+\w+\s*1([A-Z]{3})\s(.{1})\s*(\d+)(.{1})\s+.*$/) {
    print "For position $pos used conformation 1\n";
  }

#ATOM   4284  N   PHE D1001       8.859  35.692  13.572  1.00 95.55           N  
#ATOM   1478  CA  THR B 225      20.597  29.132  50.936  2     0 2.03          0.10
#ATOM   2317  N   VAL A 309     -13.102  87.192-160.762  1     0 1.63         -0.15
  $aa = $1;
  $ch = $2;
  $pos = $3;
#  $ad = $4;
  $tot++;

#  print "line $line\n";
#  print "ch $ch $take_chains{$ch}\n";
  $take_line[$i] = $take_chains{$ch};
  if ($pos>999 && $take_line[$i-1] != $take_line[$i] ){
    $adaptpos = 1;
  }
#
  $take_line[$i] = "" if ($no_PCA && ($aa eq "PCA"|| $aa eq "ACE"));
  $take_line[$i] = "" if ($pos<$range[0] || $pos>$range[1]);

}
close(PDB);

if ($tot==0) {
  print " no coordinates found!!\n";
  exit;
}

open (PDBOUT,">$pdbout") || die ("file $pdbout could not be opened");
for ($c=0;$c<@chainarray;$c++){
  for($i=0;$i<@pdbin;$i++){
    print PDBOUT "$pdbin[$i]\n" if ($take_line[$i] == $c+1);
  }
  print PDBOUT "TER                      \n" if ($add_terminus_between_chains);
}
close (PDBOUT);



###
### cems flag buster!!!!!!!!!!!!!! imported to my universe 8/99
   sub my_getopts {
 
# general purpose command line parser
# operates on global @ARGV
# call with string 'xyz:pdq' plus optional code to designate how to handle errors.
# where x y and z are command line switches that take arguments
# and p,d, &q are ones that do not take args.
# e.g. suppose contents of ARGV resulted from:
#  -x first_arg -pd -y second_arg not_an_arg also_not_arg -q -z third_arg  again_not_arg
# then return value is a hash containing:
# x => first_arg, y=>second_arg, z=>third_arg, p=>1,d=>1,q=>1
# plus an array containing (not_an_arg,also_not_arg,again_not_arg)
# $err flags if error occurs such as missing args
 
    my($argumentative,$ignore_err) = @_;
    local( $_); my($first,$rest);
    my $err = 0;
    my ($with_args,$without_args) = split /:/ , $argumentative;
 
    my $wal = length $with_args;
    my $i=0;my $pos;
    my (%arg_hash,@not_args);
    while ( $i<=$#ARGV ) {
 
    #  print "$i $ARGV[$i] \n";
    if ( $ARGV[$i] =~ /^-(.)(.*)/ ) # does it start with a dash?
         {
 
            ($first,$rest) = ($1,$2);
             $pos = index($argumentative,$first);
 
 
           if($pos < 0)
             {  # did not match any arg ?
                if ($ignore_err<2)
                 {  #
                        print STDERR "Switch -$first  no such switch \n";
                        $err = -1;  # a soft error
                        last unless ($ignore_err) ;
                 }
 
                $pos = $wal+1  #  pretend it is a non-argument type flag
             }
           if($pos <= $wal ) {  # needs an argument
                if (!$rest) #if $rest blank then arg is the next ARGV element
                  {
                        $i++; #will absorb this arg
                        if ($i>$#ARGV)
                           {  # got any more args?
                              print STDERR "Switch -$first requires an argument \n";   
                             $err=1;
                             last unless ($ignore_err >1) ;
                             $ARGV[$i] = 1;  # make up the arg for it}
                           }
                          $rest = $ARGV[$i];
                   }
              }
                else # does not need an arg
              {
                   if ($rest) # if more tags then recylce the argument remainder
                     {
                        $ARGV[$i--] = "-$rest";
                     }
                       $rest = 1;
              }
              $arg_hash{"opt_$first"} = $rest;   # remember the arg
         } # did not start with dash
           else # is not an arg, so save it
         {
           push @not_args, $ARGV[$i];
         }
           $i++;  # next arg
       } #next while
     return (\%arg_hash, \@not_args,$err);
} #end sub

