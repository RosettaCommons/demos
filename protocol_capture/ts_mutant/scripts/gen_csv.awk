#!/usr/bin/awk -f

function usage(msg) {
  if(msg) print msg > "/dev/stderr";
  print "usage: gen_csv --species <species> --pdb <file.pdb> [--aafile <aa.txt>] [<reslist.txt>]" > "/dev/stderr";
  print "generates a csv file of mutations to file.pdb from residue spec in reslist.txt" > "/dev/stderr";
  print "residue spec has residue number in first column, 3-char residue name in second column" > "/dev/stderr";
  print "will read residue spec from stdin if reslist.txt is not specified" > "/dev/stderr";
  print "  --pdb: required name of base pdb file" > "/dev/stderr";
  print "  --species: required name of species (e.g. Scer)" > "/dev/stderr";
  bail = 1;
  exit(1);
}
BEGIN {
  # args processing
  i = 1;
  while(i < ARGC && ARGV[i] ~ /^-./) {
    if(ARGV[i] ~ /^-+(usage|help)/) {
      usage();
    } else if(ARGV[i] == "--species") {
      species = ARGV[i+1];
      delete ARGV[i];
      delete ARGV[i+1];
      i += 2;
    } else if(ARGV[i] == "--pdb") {
      pdb = ARGV[i+1];
      delete ARGV[i];
      delete ARGV[i+1];
      i += 2;
    } else if(ARGV[i] == "--aafile") {
      aafile = ARGV[i+1];
      delete ARGV[i];
      delete ARGV[i+1];
      i += 2;
    } else {
      usage("error: unknown argument \"" ARGV[i] "\"");
    }
  }
  # args checking
  if(species == "") usage("must specify species file using --species");
  if(pdb == "") usage("must specify pdb file using --pdb");
  # read amino acid name conversion file
  if(aafile == "")
    aafile = ENVIRON["MODELLER_ROOT"] "/ts/predict1_scripts/aa.txt"
  if(system("test -f \"" aafile "\"") != 0)
    usage("internal error: missing required file \"" aafile "\"");
  while( (getline < aafile) > 0) {
    aastol[$1] = $2;
    aaltos[$2] = $1;
  }
  close(aafile);
  # write csv file header
  base = gensub(/\.pdb$/, "", 1, pdb);
  outfile = base ".csv";
  printf("generating %s\n", outfile);
  printf("species=%s\n", species) > outfile;
  printf("Name,Label,Mutation\n") > outfile;
  printf(base "-WT,wt,\n") > outfile;
}
$1 ~ /^[[:digit:]]*$/ {
  resi = $1;
  resn = $2;
  nat = aaltos[resn];
  for(m in aastol)
    if(m != nat) {
      mutstr = nat resi m;
      printf("+%s-%s,?,%s\n", base, mutstr, mutstr) > outfile;
    }
}
END {
  close(outfile);
  if(!bail)
    ;
}
