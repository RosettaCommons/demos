pdf("../output_files/score_vs_rmsd.pdf")
score_file <- read.table('../output_files/score.sc',skip=1,header=T)
plot(score_file$Irms, score_file$I_sc,xlab="Interface RMSD",ylab="Rosetta interface score",pch=19)
dev.off()
