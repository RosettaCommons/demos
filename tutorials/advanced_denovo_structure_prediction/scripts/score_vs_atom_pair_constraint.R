#!/blue/meilerlab/apps/Linux2/x86_64/R/2.15.0/bin/Rscript

library(ggplot2,quietly=TRUE)

args <- commandArgs(TRUE)
if(length(args)!=4)
{
        stop("USAGE:  score_vs_atom_pair_constraint.R <silent file or score_vs_atom_pair_constraint_table> <output.pdf> <max_y_value> <plot_title>")
}

infile<-args[1]
outfile<-args[2]
max_y_value<-as.numeric(args[3])
plot_title<-args[4]

data<-read.table(infile, header=TRUE)
pdf(outfile)
plot(data$atom_pair_constraint,data$score,pch=4,cex.axis=1.5,xlab="",ylab="",main="",ylim=range(min(data$score),max_y_value))
title(main=plot_title,xlab="atom_pair_constraint score (REU)",ylab="Score (REU)",cex.lab=1.5,cex.main=2)
dev.off()
