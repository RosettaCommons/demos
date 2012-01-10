# run this:
# ../../../build/demo/release/macos/10.4/32/x86/gcc/centroid_scores -database ~/svn/branches/minirosetta_database/ | grep ^SCORE > mini_scores.sc

assert = function(cond) 
{
  if( ! cond )
    print("WHINE WHINE WHINE")
}
check.scores = function() 
{
  rs = read.table('./rosetta++_scores.sc',header=T)
  ms = read.table('./mini_scores.sc',header=T)
  assert( all(rs$pdb==ms$pdb) )
  par(mfrow=c(3,4))
  par(mar=c(2,1,1,0))
  par(tck=0)
  par(mgp=c(1,0,0))
  scores = c("env",
             "pair",
             "vdw",
             "hs",
             "ss",
             "sheet",
             "cbeta",
             "rsigma",
             "hb_srbb",
             "hb_lrbb",
             "rg",
             "rama")
  for( s in scores ) {
    plot( rs[[s]]
        , ms[[s]]
        , main=paste(s,"corr:", format(cor(rs[[s]],ms[[s]]),digits=4) )
        , xlab="rosetta++ score",ylab="minirosetta score")
  }
}

pdf( file="./comp_scores.pdf", width=10, height=7.5 )
check.scores()
q(F)
