#!/blue/meilerlab/apps/Linux2/x86_64/R/2.15.0/bin/Rscript --vanilla
library(ggplot2,quietly=TRUE)

args <- commandArgs(TRUE)

data <-read.table(args[1],header=TRUE,skip=1)

data["protocol"] <- NA
data$protocol <- sapply(strsplit( as.character(data$description),'_'),"[",5)
png( paste ( args[2], "vs_rmsd.png" , sep="_" ) )

ggplot(data) + 
	aes_string(x="rmsd", y=args[2]) +
	geom_point(aes(colour=data$protocol)) + 
	xlim(0,max(data["rmsd"])) +
	xlab("CA-RMSD") +
	theme(title = element_text(paste ( args[2], "vs. rmsd"))) +
	labs(colour="protocol") +
	theme_bw()

dev.off()
