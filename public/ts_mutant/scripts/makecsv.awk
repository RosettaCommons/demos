# script to produce ml input file from new-format .csv of name/label/mutation
# called from makecsv.sh
# ASSUMPTION: only accepts one .csv as input
# usage: awk -f makecsv.awk input.csv sasa-file.txt native.sc aa.txt wt-nat-label wt-per-group ssfile
#        assumes all .sc files are in working directory
# output to stdout
#
# if wt-per-group is non-zero, wt lines are expected throughout the file at the start of each group
# otherwise an initial wt series is stored, and spit out at every group change
#
# Args index:
# 1: input.csv
# 2: sasa-file
# 3: native.sc
# 4: aa.txt
# 5: wt/nat-label
# 6: wt-per-group
# 7: merge
# 8: radius for local run series
# 9: ssfile
function error(msg) {
  print "error: " msg > "/dev/stderr";
  bail = 1;
  exit(1);
}
BEGIN {
  FS="[[:space:]]+|,";
  do_debug = 0;
  # read s-to-l aa conversion file
  aafile = ARGV[4];
  delete ARGV[4];
  checkfile(aafile);
  while( (getline < aafile) > 0)
    aa[$1] = $2;
  close(aafile);
  # read in and index sasa data
  safn = ARGV[2];
  delete ARGV[2];
  checkfile(safn);
  while( (getline < safn) > 0) {
    if($0 ~ /^RES/) {
      # starting with columns 4-5, columns are in area,total pairs
      for(i = 3; i <= NF; i+=2)
	sahead = sahead (sahead=="" ? "" : ",") gensub("(.*)ABS$", "ACC\\1", 1, $i);
    } else {
      a = $1 + 0;
      sasa[a] = $3;
      for(i = 4; i <= NF; i++)
	sasa[a] = sasa[a] "," $i;
    }
  }
  close(safn);
  # get merge flag
  merge = ARGV[7] + 0;
  delete ARGV[7];
  # get radius
  rad = ARGV[8];
  delete ARGV[8];
  # read secondary structure file if provided
  ssf = ARGV[9];
  if(ssf != "") {
    checkfile(ssf);
    sshead = "ss_S,ss_H,ss_L";
    while(getline < ssf)
      ss[$1+0] = $3;
    close(ssf);
  } else {
    sshead = "";
  }
  delete ARGV[9];
  # get label for NAT and wt lines
  wt_label = ARGV[5];
  delete ARGV[5];
  # is there one wt line for the whole file, or one per mutation group (e.g. mutations at a position)
  wt_per_group = ARGV[6] + 0;
  delete ARGV[6];
  # get source dir
  cmd = "pwd";
  cmd | getline dir;
  close(cmd);
  mod = ENVIRON["MODELLER_ROOT"];
  if(mod == dir)
    cwd = ".";
  else
    cwd = gensub(mod "/", "", 1, dir);
  delete ARGV[3];
  # init vars
  lastname = "";
  have_wt = 0;
}
/^\#/ { next; }
/^[Ss]pecies/ {
  species = gensub("[Ss]pecies(: |=)(.*)$", "\\2", 1);
  next;
}
/^corresp/ {
  corresp = 1;
  next;
}
/^Name/ {
  if(species == "")
    error("missing species name");
  if($1 == "Name+")
    global_single_mut = 1;
  else
    global_single_mut = 0;
}
!/^Name/ {
  # possible per-line single_mut determination
  if(substr($1, 1, 1) == "+") {
    name_in = substr($1, 2);
    single_mut = 1;
  } else {
    name_in = $1;
    single_mut = global_single_mut;
  }
  # strip last char (mut single-char) from name, reduce to position so NAT/wt triggering works
  if(single_mut && name_in ~ /[[:alpha:]]$/)
    name_out = substr(name_in, 1, length(name_in)-1);
  else
    name_out = name_in;
  label_in = $2;
  label_out = (label_in == "wt") ? wt_label : label_in;
  # remove ending "WT" from multiple-mutant wt runs
  if(!single_mut && label_in == "wt")
    name_out = gensub(/WT$/, "", 1, name_out);
  # define mutation fields
  for(i = 3; i <= NF; i++) {
    if($i != "")
      mut[i-2] = $i;
  }
  debug(sprintf("mut: %s", length(mut)));
  fn = name_in (rad != "0" ? ".d" rad : "") (merge ? "merge.sc" : "rescore.sc");
#  print " -- " fn > "/dev/stderr";
  checkfile(fn, "scorefile");
  # check for starting wt line
  if(!wt_per_group && ((wtc == 0 && label_in != "wt") || (wtc != 0 && label_in == "wt")))
    error("*** wt line not first line in file or extra wt line found, aborting");
  if(label_in == "wt" && !wt_per_group) {  # process wt lines the first time they appear
    if(!have_wt) {
      wtc = 0;
      while(getline ln < fn)
	if(ln ~ /description$/) {
	  if(!header) {
	    header = fixline(ln, "RES", "MUT", "", "PT", "SPECIES", sshead, sahead, "PT_val");
	    print header;
	    print "##dir:" cwd;
	    if(corresp)
	      print "##corresp";
	  }
	} else {
	  wtc++;
	  wtf[wtc] = ln;
	  have_wt = 1;
	}
    }
  } else {  # otherwise normal output
    natres = mutres = mutpos = "";
    mutcnt = 0;
    acc_act = acc_pot = prsa_act = prsa_pot = 0;
    ss_S = ss_H = ss_L = 0;
    for(i = 1; i <= length(mut); i++) {
      # skip blank entries left by saving as spreadsheet
      if(mut[i] == "")
	continue;
      mutcnt++;
      # find long-form values of nat and mutated residues
      nr = substr(mut[i], 1, 1);
      mr = substr(mut[i], length(mut[i]), 1);
      pos = gensub(/[[:alpha:]]([[:digit:]]+)[[:alpha:]]?/, "\\1", 1, mut[i]);
      if(mr ~ /[[:digit:]]/)      # wt line of local run file may lack "mutated res" character
	mr = nr;
      if(aa[nr] == "" || aa[mr] == "")
	error("*** missing residue \"" nr "\" or \"" mr "\"\nline: " $0);
      natres = natres (natres=="" ? "" : ":") aa[nr];
      mutres = mutres (mutres=="" ? "" : ":") aa[mr];
      mutpos = mutpos (mutpos=="" ? "" : ":") pos;
      # improved ACC calc: (sum(mut_actual)/sum(mut_expected))
      # single-residue value will be same as regular ACC
      n = gensub("[[:alpha:]]", "", "g", mut[i]);
      s = sasa[n];
      if(s == "")
	error("*** could not find residue \"" n "\" in sasa table, aborting");
      ss_S = or(ss_S, ssf != "" && ss[n] == "S");
      ss_H = or(ss_H, ssf != "" && ss[n] == "H");
      ss_L = or(ss_L, ssf != "" && ss[n] == "L");
      l = split(s, sav);
      for(j = 1; j <= l; j++) 
	acc[j] += sav[j];
    }
    sa = "";
    for(i = 1; i <= length(acc); i += 2)
      sa = sprintf("%s%.2f", (sa == "" ? "" : sa ","), (100*acc[i]/acc[i+1]));
    delete acc;
    if(ssf != "")
      sso = ss_S "," ss_H "," ss_L;
    else
      sso = "";
    # for now: simple NAT/wt output on change in run name
    if(name_out != lastname && !wt_per_group)
      output_wt(name_out, sso, sa);
    while(getline ln < fn)
      if(ln ~ /description$/) {
	if(!header) {
	  header = fixline(ln, "RES", "MUT", "", "PT", "SPECIES", sshead, sahead, "PT_val");
	  print header;
	  print "##dir:" cwd;
	  if(corresp)
	    print "##corresp";
	}
      } else {
	print fixline(ln, name_out, mutres, mutpos, label_in, species, sso, sa, label_out);
      }
  }
  close(fn);
  delete mut;
  lastname = name_out;
}
END {
  if(!bail)
    ;
}
function output_wt(name_out, sso, sa) {
  for(i = 1; i <= length(wtf); i++)
    print fixline(wtf[i], name_out, natres, mutpos, "wt", species, sso, sa, wt_label);
}  
function fixline(line, res, mut, pos, pt, species, sso, acclist, pt_val) {
  sub("^SCORE:[[:space:]]*", "", line);
  gsub("[[:space:]][[:space:]]*", ",", line);
  return sprintf("%s,%s,%s,%s,%s%s,%s,%s", res, (pos == "" ? mut : mut "|" pos), pt, species, (sso == "" ? "" : sso ","), acclist, line, pt_val);
}
function checkfile(fn, src) {
  if(system("test -f \"" fn "\"") != 0)
    error("*** missing required file " fn ", aborting" (src == "" ? "" : " (source: " src ")"));
}
function debug(msg) {
  if(do_debug)
    printf("=== %d === %s\n", NR, msg);
}
