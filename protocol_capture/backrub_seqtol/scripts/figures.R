# To run this script, you need to define rosetta_analysis_dir either before
# you source it, or uncomment the line below and update the path
#rosetta_analysis_dir <- "/kortemmelab/home/colin/svn/trunk/mini/analysis"

### New Plotting/Calculation Functions ###

plot.charmat <- function(charmat, col = NULL, bg = NULL, cex=1, xlim=par("usr")[1:2], ylim=par("usr")[3:4]) {
	
	xlines <- seq(xlim[1], xlim[2], length.out=ncol(charmat)+1)
	xleft <- rep(xlines[-length(xlines)], each=nrow(charmat))
	xright <- rep(xlines[-1], each=nrow(charmat))
	
	ylines <- seq(ylim[1], ylim[2], length.out=nrow(charmat)+1)
	ybottom <- rep(ylines[-length(ylines)], nrow(charmat))
	ytop <- rep(ylines[-1], nrow(charmat))
	
	xcenter <- (xleft+xright)/2
	ycenter <- (ybottom+ytop)/2
	
	if (!is.null(bg)) {
		rect(xleft, ybottom, xright, ytop, col=bg, border=NA)
	}
	
	text(xcenter, ycenter, charmat, col = col, cex = cex)
}

exp_freq_mat <- function(exp_seq, codon_multiplicity = NULL) {

	freq_mat <- matrix(nrow=length(aa1), ncol=ncol(exp_seq), dimnames=list(unname(aa1), NULL))
	
	for (i in seq_len(ncol(freq_mat))) {
		for (j in seq_len(nrow(freq_mat))) {
			freq_mat[j,i] <- sum(exp_seq[,i] == rownames(freq_mat)[j])/nrow(exp_seq)
		}
	}
	
	if (!is.null(codon_multiplicity)) {
		correctfunc <- function(x, codon_multiplicity) {
			x <- x/codon_multiplicity
			x/sum(x)
		}
		freq_mat <- apply(freq_mat, 2, correctfunc, codon_multiplicity[rownames(freq_mat)])
	}
	
	freq_mat
}

plot_seqrank <- function(freq_mat, exp_freq_mat = NULL, wt_seq = NULL, star_mat = NULL, rank_line = 0, wt_col = "red", other_col = "black") {

	char_mat <- matrix(nrow=nrow(freq_mat), ncol=ncol(freq_mat))
	col_mat <- matrix(other_col, nrow=nrow(freq_mat), ncol=ncol(freq_mat))
	bg_freq_mat <- matrix(nrow=nrow(freq_mat), ncol=ncol(freq_mat))
	
	for (i in seq_len(ncol(freq_mat))) {

		char_mat[,i] <- rownames(freq_mat)[order(freq_mat[,i])]
		if (!is.null(star_mat)) {
			star_mat[,i] <- rev(star_mat[char_mat[,i],i])
		}
		for (j in seq_len(nrow(freq_mat))) {
		
			if (is.null(exp_freq_mat)) {
				bg_freq_mat[j,i] <- freq_mat[char_mat[j,i],i]
			} else {
				bg_freq_mat[j,i] <- exp_freq_mat[char_mat[j,i],i]
			}
			
			if (!is.null(wt_seq)) {
				if (char_mat[j,i] == wt_seq[i]) col_mat[j,i] <- wt_col
			}
		}
	}
	
	col_levels <- seq(0,1,by=.1)
	col_levels <- seq(0,ceiling(max(bg_freq_mat)/.1)*.1,by=.1)
	#cols <- gray(seq(1,0,length.out=length(col_levels)-1))
	cols <- rev(c(topo.colors(length(col_levels)-2), "white"))
	
	bg_mat <- matrix(cols[pmin(floor((bg_freq_mat)*length(cols)/max(col_levels))+1, length(cols))], nrow=nrow(bg_freq_mat))
	
	op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    mar1 <- c(0.6, 2.7, 2.7, 0.4)
    mar2 <- c(0.6, 0.4, 2.7, 3.4)
    
    devwidth <- par("din")[1]*2.54
    charheight <- par("cin")[2]*2.54
    width1 <- (mar1[2]+mar1[4])*charheight
    width2 <- (mar2[2]+mar2[4])*charheight
    boxwidth <- (devwidth - sum(width1+width2))/(ncol(freq_mat)+1)
    layout(matrix(1:2, nrow=1,ncol=2), widths=c(width1+boxwidth*ncol(freq_mat),width2+boxwidth))
	
	par(mar=mar1, mgp=c(1.5, .25, 0), cex=1)
	
	plot(0, 0, type="n", xlim=c(0.5,0.5+ncol(freq_mat)), ylim=c(20.5,0.5), xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="", ylab="")
	plot.charmat(char_mat, col_mat, bg_mat)
	mtext("Predicted Rank", 2, par("mgp")[1])
	axis(2, 1:20, tick=FALSE, las=2)
	mtext("Residue", 3, par("mgp")[1])
	residx <- seq(1, ncol(freq_mat), by=2)
	axis(3, residx, colnames(freq_mat)[residx], tick=FALSE)
	residx <- seq(2, ncol(freq_mat), by=2)
	axis(3, residx, colnames(freq_mat)[residx], tick=FALSE)
	box(lwd=.5)
    
    if (!is.null(star_mat)) {
    	points(t(t(which(t(star_mat), arr.ind=TRUE))+c(.3,0)), pch="*")
    }
    
    if (rank_line) {
    	abline(h=rank_line+0.5, lty="dashed")
    }
	
	maradj <- (1-length(cols)/nrow(col_mat))*0.5*par("pin")[2]/par("cin")[2]
	mar2[1] <- mar2[1]+maradj
	mar2[3] <- mar2[3]+maradj
	
	par(mar=mar2, mgp=c(2.2, .25, 0), cex=1)
	
	plot.new()
    plot.window(xlim = c(0, 1), ylim = range(col_levels), xaxs = "i", yaxs = "i")
    rect(0, col_levels[-length(col_levels)], 1, col_levels[-1L], col = cols, lwd=.5)
    axis(4, col_levels[seq(1,length(col_levels),1)], paste(round(col_levels[seq(1,length(col_levels),1)]*100), "%", sep=""), tick=FALSE, las=2)
    bg_title <- "Predicted Frequency"
    if (!is.null(exp_freq_mat)) bg_title <- "Experimental Frequency"
    mtext(bg_title, 4, par("mgp")[1])
    box(lwd=.5)
    
    invisible(bg_freq_mat)
}

plot_gen_contrib <- function(entitieslist, generationslist, 
                             fitness_coef = c(1/2.5, 1/2.5, 1/2.5, 1),
                             temp_or_thresh = 0.228, 
                             type = c("boltzmann", "cutoff"),
                             main = "") {

	type <- match.arg(type)

	names(fitness_coef) <- paste("state1_fitness_comp", seq_along(fitness_coef), sep="")

	gen_contrib <- matrix(nrow=length(entitieslist), ncol=length(generationslist[[1]]))

	for (i in seq_along(entitieslist)) {

		fitness <- entities_fitness(entitieslist[[i]], fitness_coef)
		min_fitness <- min(fitness)
		
		if (type == "cutoff") {
			weight <- fitness <= min_fitness+temp_or_thresh
		} else {
			if (temp_or_thresh != 0) {
				weight <- exp(-(fitness-min_fitness)/temp_or_thresh)
			} else {
				weight <- fitness == min_fitness
			}
		}
		
		weight <- weight/sum(weight)
		first_gen <- integer(length(fitness))
		for (j in rev(seq_along(generationslist[[i]]))) {
			first_gen[generationslist[[i]][[j]]] <- j
		}
		for (j in seq_along(generationslist[[i]])) {
			gen_contrib[i,j] <- sum(weight[first_gen == j])
		}
	}
	
	colnames(gen_contrib) <- seq_along(generationslist[[1]])
	
	boxplot(as.data.frame(gen_contrib), xlab="Generation", ylab="Sequence Contribution", main=main, ylim=c(0,1), yaxt="n")
	axis(2, seq(0, 1, by=.25), labels=FALSE)
	axis(2, seq(0, 1, by=.5), seq(0, 1, by=.5), tick=FALSE, line=FALSE)
}

write.seqdata <- function(seqs, filename, genename, domainnum = 1, domainrange = "79-83") {

	write(paste(paste("Gene Name", genename, sep = "\t"),
                "Accession	Refseq:NP_004078 UniProt:Q12959",
                "Organism	",
                "NCBITaxonomyID	9606",
                paste("Domain Number", domainnum, sep = "\t"),
                "Domain Type	",
                "Interpro ID	",
                "Technique	",
                "Domain sequence	",
                paste("Domain Range", domainrange, sep = "\t"),
                "Comment	
PeptideName	Peptide	CloneFrequency	QuantData	ExternalIdentifier",
                sep = "\n"),
          filename)
	seqs <- seqs[seqs != "NANANANANA"]
	write(paste(seq_along(seqs), seqs, rep(1, length=length(seqs)), sep = "\t", collapse = "\n"), filename, append = TRUE)
}

write.plastarray <- function(plastarray, dirpath, numseq = 100) {

	if (!file.exists(dirpath)) dir.create(dirpath, recursive = TRUE)
	
	protnames <- dimnames(plastarray)$prot
	
	for (protname in protnames) {
	
		domainname <- sub("^(.+)-([0-9]+)$", "\\1", protname)
		domainnum <- sub("^(.+)-([0-9]+)$", "\\2", protname)
		if (domainname == domainnum) domainnum <- 1
		#print(c(domainname, domainnum))
		pwm <- plastarray[,,protname]
		seqmat <- matrix(character(), nrow=numseq, ncol=ncol(pwm))
		for (i in seq_len(ncol(pwm))) {
			colfun <- stepfun(c(0,cumsum(pwm[,i])), c(1,seq_along(pwm[,i]),length(pwm[,i])))
			funx <- seq(0, 1, length.out=numseq+1)
			funx <- funx[-1] - mean(diff(funx))/2
			seqmat[,i] <- names(pwm[,i])[colfun(funx)]
		}
		seqs <- apply(seqmat, 1, paste, collapse = "")
		write.seqdata(seqs, file.path(dirpath, paste(protname, ".txt", sep = "")), paste("Predicted", domainname), domainnum)
		
	}
	
	filenames <- paste(protnames, ".txt", sep = "")
	
	write(paste(c("#ProjectFile", filenames), sep = "\n"), file.path(dirpath, "project.txt"))
}

numbits <- function(x) sum((x*log(x, 2))[x!=0])-log(1/length(x), 2)

roc <- function(truth, pred=NULL, thresh, plot_roc=FALSE) {

	if (is.complex(truth) && is.null(pred)) {
		pred <- Im(truth)
		truth <- Re(truth)
	}

	truth <- (truth >= thresh)[order(pred, decreasing=TRUE)]
	predrank <- sort(rank(-pred))
	
	tp <- cumsum(truth)/sum(truth)
	fp <- cumsum(!truth)/sum(!truth)
	
	tp <- tp[rev(!duplicated(rev(predrank)))]
	fp <- fp[rev(!duplicated(rev(predrank)))]
	
	tp <- c(0, tp)
	fp <- c(0, fp)
	
	if (plot_roc) plot(fp, tp, xlim=c(0,1), ylim=c(0,1), type="b")
	
	if (sum(truth)) sum(diff(fp)*(tp[-1]-.5*diff(tp))) else NA
}

auc_array <- function(truefreq, predfreq, thresh, margin) {

	comb <- truefreq
	comb[seq_along(comb)] <- complex(real=truefreq, imaginary=predfreq)
	
	apply(comb, margin, roc, thresh=thresh)
}

toprank <- function(exparray, predarray) {

	expbest <- apply(exparray, 2:3, function(x) paste(which(x==max(x)), collapse=" "))
	predrank <- apply(-predarray, 2:3, rank, ties.method="max")
	
	topmat <- exparray[1,,]
	
	for (i in seq_len(dim(exparray)[2])) {
		for (j in seq_len(dim(exparray)[3])) {
			x <- exparray[,i,j]
			bestpos <- which(x==max(x))
			topmat[i,j] <- mean(predrank[bestpos,i,j])
		}
	}
	
	topmat
}

ranktop <- function(truefreq, predfreq = NULL) {

	if (is.complex(truefreq) && is.null(predfreq)) {
		predfreq <- Im(truefreq)
		truefreq <- Re(truefreq)
	}

	predrank <- rank(-predfreq, ties.method="max")
	bestpos <- which(truefreq==max(truefreq))
	mean(predrank[bestpos])
}

ranktop_array <- function(truefreq, predfreq, margin) {

	comb <- truefreq
	comb[seq_along(comb)] <- complex(real=truefreq, imaginary=predfreq)
	
	apply(comb, margin, ranktop)
}

pred_stats <- function(truefreq, predfreq) {

	stats_mat <- matrix(nrow=4, ncol=ncol(truefreq))
	rownames(stats_mat) <- c("aad", "auc", "frobenius", "ranktop")
	colnames(stats_mat) <- colnames(truefreq)
	
	stats_mat["auc",] <- auc_array(truefreq, predfreq, 0.1, 2)
	stats_mat["ranktop",] <- ranktop_array(truefreq, predfreq, 2)
	
	# Make sure experimental frequencies are normalized
	truefreq <- apply(truefreq, 2, function(x) x/sum(x))
	
	stats_mat["aad",] <- apply(abs(truefreq-predfreq), 2, mean)
	stats_mat["frobenius",] <- sqrt(apply((truefreq-predfreq)^2, 2, sum))
	
	stats_mat
}

stats_table <- function(group_expfreq_list, group_predfreq_list) {

	dframe <- data.frame()
	
	for (i in seq_along(group_expfreq_list)) {
	
		stats_list <- vector("list", length(group_expfreq_list[[i]]))
		
		for (j in seq_along(group_expfreq_list[[i]])) {
			stats_list[[j]] <- pred_stats(group_expfreq_list[[i]][[j]], group_predfreq_list[[i]][[j]])
		}
		
		dnames <- c(dimnames(group_expfreq_list[[i]][[1]]), list(sim=names(group_expfreq_list[[i]])))
		expfreqarray <- array(do.call(c, group_expfreq_list[[i]]), sapply(dnames, length), dnames)
		predfreqarray <- array(do.call(c, group_predfreq_list[[i]]), sapply(dnames, length), dnames)
		
		repfrac <- apply(apply(predfreqarray, 2:3, function(x) x > sort(x)[20-5]) & expfreqarray >= .1, 2:3, sum)/apply(expfreqarray >= .1, 2:3, sum)
		# Exclude those without any with frequencies greater than 10%
		repfrac <- repfrac[!is.nan(repfrac)]
		
		# Make sure experimental frequencies are normalized
		expfreqarray <- apply(expfreqarray, 2:3, function(x) x/sum(x))
		
		dfrow <- data.frame(
			proteins=length(group_expfreq_list[[i]]),
			res_pos=sum(sapply(group_expfreq_list[[i]], ncol)),
			bits_phage=mean(apply(expfreqarray, 2:3, numbits)),
			bits_pred=mean(apply(predfreqarray, 2:3, numbits)),
			frac_top5=mean(repfrac)*100,
			aad=mean(sapply(stats_list, "[", "aad", seq_len(ncol(stats_list[[1]]))))*100,
			auc=mean(sapply(stats_list, "[", "auc", seq_len(ncol(stats_list[[1]]))), na.rm=TRUE),
			frobenius=mean(sapply(stats_list, "[", "frobenius", seq_len(ncol(stats_list[[1]])))),
			ranktop=mean(sapply(stats_list, "[", "ranktop", seq_len(ncol(stats_list[[1]]))))
		)
		
		dframe <- rbind(dframe, dfrow)
	}
	
	rownames(dframe) <- names(group_expfreq_list)
	
	dframe
}

# This is a stats_table test function which regenerates JMB 2010 Table 2
jmb2010_table2 <- function() {

	groupfiles <- rev(list.files("../data/jmb2010_groups", full.names=TRUE))
	groups <- vector("list", length(groupfiles))
	names(groups) <- sub(".+/(.+).txt", "\\1", groupfiles)
	for (i in seq_along(groupfiles)) groups[[i]] <- read.table(groupfiles[i], as.is=TRUE)
	
	group_expfreq_list <- vector("list", length(groupfiles))
	group_predfreq_list <- vector("list", length(groupfiles))
	
	names(group_expfreq_list) <- names(groups)
	names(group_predfreq_list) <- names(groups)
	
	for (i in seq_along(groups)) {
	
		group_expfreq_list[[i]] <- vector("list", nrow(groups[[i]]))
		group_predfreq_list[[i]] <- vector("list", nrow(groups[[i]]))
		
		names(group_expfreq_list[[i]]) <- groups[[i]][,1]
		names(group_predfreq_list[[i]]) <- groups[[i]][,1]
		
		for (j in seq_len(nrow(groups[[i]]))) {
		
			group_expfreq_list[[i]][[j]] <- as.matrix(read.table(paste("../data/exp_pwm/", groups[[i]][j,1], ".txt", sep=""), check.names=FALSE))
			group_predfreq_list[[i]][[j]] <- as.matrix(read.table(paste("../data/jmb2010_pwm/", groups[[i]][j,1], ".txt", sep=""), check.names=FALSE))
		}
	}
	
	stats_table(group_expfreq_list, group_predfreq_list)
}

source(file.path(rosetta_analysis_dir, "apps", "sequence_tolerance.R"))

### Updated Functions from sequence_tolerance.R ###

# Read a list of *.ga.entities checkpoint files out of a directory. By default,
# the parsed data is saved in the R format

read_ga_entities_list <- function(dirpath, filepattern = NULL, recompute = FALSE, savedata = FALSE, readgen = FALSE) {

	filename <- file.path(dirpath, paste("entities", filepattern, ".Rda", sep = ""))
	
	if (file.exists(filename) && !recompute) {
		load(filename)
	} else {
		simpattern <- if (is.null(filepattern)) "" else filepattern
		
		entitiesfiles <- list.files(dirpath, pattern = paste(simpattern, ".*\\.ga\\.entities", sep = ""),
		                            full.names = TRUE, recursive = TRUE)
		entitiesfiles <- entitiesfiles[file.info(entitiesfiles)$size > 0]
		
		entitieslist <- vector("list", length(entitiesfiles))
		generationslist <- vector("list", length(entitiesfiles))
		
		for (i in seq(along = entitiesfiles)) {
			print(entitiesfiles[i])
			entitieslist[[i]] <- read_ga_entities(entitiesfiles[i], unname(aa1))
			if (readgen) generationslist[[i]] <- read_ga_generations(sub("entities", "generations", entitiesfiles[i]), entitieslist[[i]])
		}
		
		if (readgen) attr(entitieslist, "generations") <- generationslist
		
		if (savedata == TRUE) save(entitieslist, file = filename)
	}
	
	entitieslist
}

entities_list_pwms <- function(entities, fitness_coef = c(1/2.5, 1/2.5, 1/2.5, 1),
                               temp_or_thresh = 0.228, 
                               type = c("boltzmann", "cutoff")) {

	type <- match.arg(type)
	if (is.null(names(fitness_coef))) {
		names(fitness_coef) <- paste("state1_fitness_comp", seq_along(fitness_coef), sep="")
	}

	pwms <- entities_pwms(entities, temp_or_thresh, fitness_coef, type)
	
	posnames <- colnames(entities[[1]])[seq_along(pwms)]
	posnames <- gsub("AA", "", posnames)
	
	pwmsdimnames <- list(aa=rownames(pwms[[1]]), rep=NULL, pos=posnames)
	
	pwms <- array(do.call("c", pwms), dim = c(dim(pwms[[1]]), length(pwms)))
	dimnames(pwms) <- pwmsdimnames
	
	pwms
}

collapse_pwms <- function(pwms, percentile = .5) {

	pwm <- apply(pwms, c(1,3), quantile, percentile)
	
	minnotzero <- function(x) {
		x <- x[x!=0]
		if (length(x)) return(min(x))
		NA
	}
	plastmin <- apply(pwms, c(1,3), minnotzero)
	correcteddist <- apply(plastmin, 2, function(x) as.numeric(!is.na(x) & x==min(x, na.rm = TRUE)))
	
	for (i in which(colSums(pwm) == 0)) {
		#print(paste("Correcting", i))
		pwm[,i] <- correcteddist[,i]
	}
	
	pwm <- apply(pwm, 2, function(x) x/sum(x))
	
	pwm
}

# This function is the main data processing procedure. It takes a directory 
# path which contains *.ga.entities files. It reads all those files and
# produces a set of boxplots in several different file formats. It also 
# generates a position weight matrix and FASTA file for producing a sequence 
# logo.

process_specificity <- function(dirpath = ".", fitness_coef = c(1/2.5, 1/2.5, 1/2.5, 1),
                                temp_or_thresh = 0.228, 
                                type = c("boltzmann", "cutoff"),
                                percentile = .5) {

	type <- match.arg(type)
	names(fitness_coef) <- paste("state1_fitness_comp", seq_along(fitness_coef), sep="")
	
	entities <- read_ga_entities_list(dirpath)
	pwms <- entities_list_pwms(entities, fitness_coef, temp_or_thresh, type)
	pwm <- collapse_pwms(pwms, percentile)
	posnames <- colnames(pwm)
	
	backruboutput <- file.path(dirpath, "backrub_0_stdout.txt")
	if (file.exists(backruboutput)) {
		backrubcmd <- readLines(backruboutput, 1)
		startpdbfile <- gsub("^.+ -s ([^ ]+) .+$", "\\1", backrubcmd)
		pdbseq <- pdb_sequence(file.path(dirpath, startpdbfile))
		colnames(pwm) <- names(pdbseq)[as.integer(posnames)]
	}
	
	write.table(pwm, "specificity_pwm.txt", quote=FALSE, sep="\t", col.names=NA)
	
	seqmat <- pwm_to_seqmat(pwm)
	
	cat(paste(">", seq_len(nrow(seqmat)), "\n", apply(seqmat, 1, paste, collapse=""), sep=""), file="specificity_sequences.fasta", sep="\n")
	
	plotwidth <- 7
	plotheight <- 3
	pointsize <- 12
	
	pdf("specificity_boxplot.pdf", width=plotwidth, height=plotheight, pointsize=pointsize)
	pdfdev <- dev.cur()
	png("specificity_boxplot.png", width=plotwidth*72, height=plotheight*72*length(posnames), pointsize=3/2*pointsize)
	pngdev <- dev.cur()
	par(mfrow=c(length(posnames), 1))
	
	for (i in seq_along(posnames)) {
		
		for (type in c("pdf", "png", "pngsep")) {
			
			if (type == "pdf") dev.set(pdfdev)
			if (type == "png") dev.set(pngdev)
			if (type == "pngsep") png(paste("specificity_boxplot_", posnames[i],".png", sep=""), width=plotwidth*72, height=plotheight*72, pointsize=pointsize)
			
			par(mar = c(2.8, 2.8, 1.5, 0.1), mgp = c(1.7, 0.6, 0))
			main <- paste("Residue", posnames[i], "Specificity Boxplot")
			plot(0, 0, type="n", xlim=c(1,20), ylim=c(0,1), main=main, xlab="Amino Acid", ylab="Predicted Frequency", axes=FALSE)
			abline(h=seq(0, 1, by=.2), col="gray")
			boxplot(as.data.frame(t(pwms[,,i])), col="white", add=TRUE)
			points(1:20, pwm[,i], pch=4, col="blue")
			
			if (type == "pngsep") dev.off()
		}
	}
	dev.off(pdfdev)
	dev.off(pngdev)
}

plot.scorearray <- function(scorearray, main = "", levels = seq(0,10,length.out=21), xlab = "Percentile Cutoff", ylab = "Score Cutoff", zlab="Mean %\nDifference", params = NULL, fitparams = NULL, lowpoint = TRUE) {

	xpoints <- as.numeric(gsub("[^0-9.]", "", dimnames(scorearray)[[1]]))
	ypoints <- as.numeric(gsub("[^0-9.]", "", dimnames(scorearray)[[2]]))
	lowind <- which(scorearray == min(scorearray), arr.ind=TRUE)
	#if (is.null(zlim)) zlim <- range(scorearray, finite=TRUE)
	filled.contour(xpoints, ypoints, scorearray, xaxt="n", yaxt="n", levels=levels,
	               color.palette = topo.colors, zlim = range(levels),
	               plot.axes = { abline(v = xpoints, h = ypoints, col = "black", lwd=.2)
	                             if (!is.null(params)) {
	                             	segments(params[1], params[3], params[2], params[4], lwd=2, col="red")
	                             }
	                             contour(xpoints, ypoints, scorearray, zlim = range(levels), levels = levels, labcex = 1, add=TRUE)
	                             if (lowpoint) {
	                             	points(xpoints[lowind[1]], ypoints[lowind[2]], cex=3, pch=20)
	                             	text(xpoints[lowind[1]], ypoints[lowind[2]], round(scorearray[lowind[1],lowind[2]], 2), adj=-c(.4,.4))
	                             }
	                             axis(1)#, xpoints, dimnames(scorearray)[[1]]) 
	                             axis(2)#, ypoints, dimnames(scorearray)[[2]])
	                             if (!is.null(fitparams)) {
	                             	points(fitparams[1], fitparams[2], cex=3, pch=20, col="red")
	                             	text(fitparams[1], fitparams[2], round(interp.array(scorearray, list(xpoints, ypoints), matrix(fitparams, nrow=1)), 2), adj=-c(.4,.4))
	                             
	                             }
	                           },
	               plot.title = title(main, xlab = xlab, ylab = ylab), 
	               key.title = title(main=zlab, xpd=NA),
	               key.axes = axis(4, levels) )
	xlen <- length(dimnames(scorearray)[[1]])
	ylen <- length(dimnames(scorearray)[[2]])
	#axis(1, seq(0,1,length.out=xlen), dimnames(scorearray)[[1]])
	#axis(2, seq(0,1,length.out=ylen), dimnames(scorearray)[[2]])
	#print(lowind)
	#print(c(xpoints[lowind[1]], ypoints[lowind[2]]))
	
	lowind
}

dftomat <- function(df, xcol, ycol, zcol) {

	df <- df[!duplicated(df[,c(xcol,ycol)]),]
	xvals <- sort(unique(as.character(df[,xcol])))
	yvals <- sort(unique(as.character(df[,ycol])))
	mat <- matrix(NA, nrow=length(xvals), ncol=length(yvals))
	rownames(mat) <- xvals
	colnames(mat) <- yvals
	mat[cbind(match(as.character(df[,xcol]),xvals),match(as.character(df[,ycol]),yvals))] <- df[,zcol]
	mat
}

pdb_sequence <- function(pdbpath) {

	if (length(grep("\\.gz$", pdbpath))) {
		con <- gzfile(pdbpath)
		pdblines <- readLines(con)
		close(con)
	} else {
		pdblines <- readLines(pdbpath)
	}
	calines <- grep("^ATOM  .....  CA  ... ......", pdblines, value=TRUE)
	chains <- substr(calines, 22, 22)
	aatypes <- substr(calines, 18, 20)
	resnums <- gsub(" ", "", substr(calines, 23, 27))
	
	seqs <- vector("list", length(unique(chains)))
	names(seqs) <- unique(chains)
	
	for (i in seq_along(seqs)) {
		idx <- which(chains == names(seqs)[i])
		seqs[[i]] <- structure(aatypes[idx], names=resnums[idx])
	}
	
	seqs
}

start_sequence <- function(dirpath) {

	pdbfile <- list.files(dirpath, "_low.pdb", full.names=TRUE)[1]
	pdbseq <- pdb_sequence(pdbfile)
	
	entfile <- list.files(dirpath, ".ga.entities", full.names=TRUE)[1]
	entities <- read_ga_entities(entfile)
	seqnums <- as.integer(sub("AA", "", grep("^AA", colnames(entities), value=TRUE)))
	
	structure(aa1[do.call(c, pdbseq)[seqnums]], names=seqnums)
}

naive_pwm <- function(startseq) {

	naive_groups <- list(
		c("D", "E", "N", "Q"), 
		c("R", "K", "H"), 
		c("L", "I", "V", "M"),
		c("F", "Y", "W"),
		c("P", "A", "G"),
		c("S", "T"),
		c("C")
	)

	pwm <- matrix(0, nrow=length(aa1), ncol=length(startseq))
	rownames(pwm) <- unname(aa1)
	colnames(pwm) <- names(startseq)
	
	for (i in seq_along(startseq)) {
		for (j in seq_along(naive_groups)) {
			if (startseq[i] %in% naive_groups[[j]]) {
				pwm[naive_groups[[j]],i] <- 1/length(naive_groups[[j]])
			}
		}
	}
	
	pwm
}

if (length(list.files("1N7T", "ga.entities")) != 2000) {

### Figure 2 ###

cwd <- getwd()

setwd("2QMT")
if (!file.exists("specificity_pwm.txt")) process_specificity(fitness_coef=c(1/2.5,1/2.5))
setwd(cwd)

pred_freq_2QMT <- as.matrix(read.table("2QMT/specificity_pwm.txt"))
colnames(pred_freq_2QMT) <- gsub("[^0-9]", "", colnames(pred_freq_2QMT))

seqpos <- c(5,7,16,18,30,33)

wtseq <- strsplit("PNRFRGKDLAGSPGGTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTEGGGSSGGGADFDYEKMANAN", "")[[1]]

round3seq <- readLines("../data/gb1_round_3_seq.txt")

round3mat <- strsplit(round3seq, "")
round3mat <- do.call(rbind, round3mat[sapply(round3mat, length) == 90])

round3good <- rep(TRUE, nrow(round3mat))

for (i in seq_along(round3good)) {
	if (any(round3mat[i,seqpos+14] == "X")) round3good[i] <- FALSE
	if (any(round3mat[i,-(seqpos+14)] != wtseq[-(seqpos+14)])) round3good[i] <- FALSE
}

round3mat <- round3mat[round3good,seqpos+14]

exp_freq_2QMT <- exp_freq_mat(round3mat)
colnames(exp_freq_2QMT) <- seqpos

pdf("Figure2_2QMTSeqRank.pdf", FALSE, "Arial", width=3.27, height=4, pointsize=10)
#postscript("Figure2_2QMTSeqRank.eps", FALSE, "Arial", width=3.27, height=4, pointsize=10, horizontal=FALSE)
ranked_exp_freq <- plot_seqrank(pred_freq_2QMT, exp_freq_2QMT, wtseq[seqpos+14], rank_line=5)
dev.off()


### Figure 1 & 2 Sequence Logos ###

if (! "ent_gen_2QMT" %in% ls()) ent_gen_2QMT <- read_ga_entities_list("2QMT", readgen = TRUE)

if ("ent_gen_2QMT" %in% ls()) {

if (! "kT_2QMT" %in% ls()) {

	optfunc_2QMT <- function(kT) {
		pred_freq <- collapse_pwms(entities_list_pwms(ent_gen_2QMT, c(0.4, 0.4), kT))
		aad <- mean(abs(pred_freq-exp_freq_2QMT))
		cat("kT=", kT, " AAD=", aad*100, "%\n", sep="")
		aad
	}
	kT_2QMT <- optimize(optfunc_2QMT, c(0.2,0.8))$minimum
	kT_2QMT_label <- format(kT_2QMT, digits=2)
}

pred_freq_2QMT_opt <- collapse_pwms(entities_list_pwms(ent_gen_2QMT, c(0.4, 0.4), kT_2QMT))

freqarray_2QMT <- array(numeric(1), dim=c(dim(pred_freq_2QMT), length(ent_gen_2QMT)+2), dimnames=c(dimnames(pred_freq_2QMT), list(prot=c(paste("2QMT", c(seq_along(ent_gen_2QMT), 0), sep="-"), "2QMT_Opt"))))

freqarray_2QMT[rownames(pred_freq_2QMT),,"2QMT-0"] <- as.matrix(pred_freq_2QMT)
freqarray_2QMT[rownames(pred_freq_2QMT),,"2QMT_Opt"] <- pred_freq_2QMT_opt

fitness_coef <- c(0.4, 0.4)
names(fitness_coef) <- paste("state1_fitness_comp", seq_along(fitness_coef), sep="")

for (i in seq_along(ent_gen_2QMT)) {
	freqarray_2QMT[,,i] <- entities_pwm(ent_gen_2QMT[[i]], 0.228, fitness_coef)
}

write.plastarray(freqarray_2QMT, "2QMT_LOLA")
write.seqdata(apply(round3mat, 1, paste, collapse=""), "2QMT_LOLA/GB1.txt", "GB1", domainnum = 1, domainrange = "79-83")
cat("GB1.txt", file="2QMT_LOLA/project.txt", append=TRUE)

pdf("Figure2_2QMTSeqRank_Opt.pdf", FALSE, "Arial", width=3.27, height=4, pointsize=10)
#postscript("Figure2_2QMTSeqRank_Opt.eps", FALSE, "Arial", width=3.27, height=4, pointsize=10, horizontal=FALSE)
ranked_exp_freq <- plot_seqrank(pred_freq_2QMT_opt, exp_freq_2QMT, wtseq[seqpos+14], rank_line=5)
dev.off()

} # end if


### Figure 3 and S4 ###

hgh_freq <- vector("list", 6)
for (i in 1:6) {
	dname <- paste("1A22_", i, sep="")
	setwd(dname)
	if (!file.exists("specificity_pwm.txt")) process_specificity(fitness_coef=c(1/2.5, 1/2.5, 1/2.5, 1))
	setwd(cwd)
	hgh_freq[[i]] <- read.table(file.path(dname, "specificity_pwm.txt"))
	colnames(hgh_freq[[i]]) <- gsub("[^0-9]", "", colnames(hgh_freq[[i]]))
}
hgh_freq <- as.matrix(do.call(cbind, hgh_freq))
hgh_freq <- hgh_freq[,order(as.numeric(colnames(hgh_freq)))]

hgh_bind_freq <- t(as.matrix(read.table("../data/hgh_bind_freq.txt", TRUE, row=1)))/100
hgh_wt_seq <- gsub("[0-9]", "", colnames(hgh_bind_freq))
colnames(hgh_bind_freq) <- gsub("[A-Z]", "", colnames(hgh_bind_freq))
colnames(hgh_freq) <- colnames(hgh_bind_freq) 

postscript("FigureS4_1A22SeqRank.eps", FALSE, "Arial", width=13.83, height=4, pointsize=10, horizontal=FALSE, paper="special")
hgh_bind_freq_ordered <- plot_seqrank(hgh_freq, hgh_bind_freq, hgh_wt_seq, rank_line=5)
dev.off()

hgh_table_pos <- match(c("67", "171", "176", "61", "178", "183", "167", "172", "25", "179", "21", "45", "64", "22", "48", "42"), colnames(hgh_freq))

hgh_comp_lib <- list(
	"67" =c("A", "V", "S", "G", "D"),
	"171"=c("H", "S"),
	"176"=c("Y", "H", "M", "L"),
	"61" =c("A", "G", "S"),
	"178"=c("M", "K", "Q", "L"),
	"183"=c("F", "Y", "H", "M"),
	"167"=c(),
	"172"=c("M", "I", "R", "Q", "A", "V"),
	"25" =c("Y", "W", "R"),
	"179"=c("H", "V", "K", "M"),
	"21" =c("F", "Y"),
	"45" =c("Y", "M", "F", "H"),
	"64" =c("F", "K", "H"),
	"22" =c("H", "N", "S"),
	"48" =c("F", "Y", "H"),
	"42" =c("F", "H", "K", "R")
)

hgh_star_mat <- matrix(nrow=nrow(hgh_freq), ncol=length(hgh_comp_lib))
rownames(hgh_star_mat) <- rownames(hgh_freq)
for (i in seq_len(ncol(hgh_star_mat))) {
	hgh_star_mat[,i] <- rownames(hgh_freq) %in% hgh_comp_lib[[i]]
}

postscript("Figure3_1A22SeqRankTable.eps", FALSE, "Arial", width=6.83, height=4, pointsize=10, horizontal=FALSE, paper="special")
hgh_bind_freq_ordered <- plot_seqrank(hgh_freq[,hgh_table_pos], hgh_bind_freq[,hgh_table_pos], hgh_wt_seq[hgh_table_pos], hgh_star_mat, rank_line=5)
dev.off()


### Figure 4 ###

pdz_runs <- c("2I0L_A_C_V2006", "2IWP_B_A_V1927", "2FNE_A_C_V2048", "1N7T", "1N7T_V83K")

pdz_freq <- vector("list", length(pdz_runs))
names(pdz_freq) <- pdz_runs
for (i in seq_along(pdz_runs)) {
	dname <- pdz_runs[i]
	setwd(dname)
	boltzmann_factor <- 0.228
	if (dname == "1N7T_V83K") boltzmann_factor <- boltzmann_factor + 0.021
	if (!file.exists("specificity_pwm.txt")) process_specificity(fitness_coef=c(1/2.5, 1/2.5, 1/2.5, 1), temp_or_thresh=boltzmann_factor)
	setwd(cwd)
	pdz_freq[[i]] <- as.matrix(read.table(file.path(dname, "specificity_pwm.txt")))
	colnames(pdz_freq[[i]]) <- -4:0
}

dnames <- c(dimnames(pdz_freq[[1]]), list(prot=pdz_runs))
freqarray_pdz <- array(do.call(c, pdz_freq), sapply(dnames, length), dnames)

write.plastarray(freqarray_pdz, "PDZ_LOLA")


### Figure 5 ###

if (!"ent_gen_2I0L" %in% ls()) ent_gen_2I0L <- read_ga_entities_list("2I0L_A_C_V2006", readgen = TRUE)

if ("ent_gen_2I0L" %in% ls()) {

pdf("Figure5_2I0LGenContrib.pdf", FALSE, "Arial", width=2.2, height=2.2, pointsize=10)
#postscript("Figure5_2I0LGenContrib.eps", FALSE, "Arial", width=3.27, height=3.27, pointsize=10, horizontal=FALSE)
par(mar=c(2.7,2.7,1.5,0.2), mgp=c(1.7, 0.6, 0))
plot_gen_contrib(ent_gen_2I0L, attr(ent_gen_2I0L, "generations"), main="PDZ/Peptide Interface")
mtext("A", 3, .4, adj=-0.26, cex=12/10, font=2)
dev.off()

} # end if

if ("ent_gen_2QMT" %in% ls()) {

pdf("Figure5_2QMTGenContrib.pdf", FALSE, "Arial", width=2.2, height=2.2, pointsize=10)
#postscript("Figure5_2QMTGenContrib.eps", FALSE, "Arial", width=3.27, height=3.27, pointsize=10, horizontal=FALSE)
par(mar=c(2.7,2.7,1.5,0.2), mgp=c(1.7, 0.6, 0))
plot_gen_contrib(ent_gen_2QMT, attr(ent_gen_2QMT, "generations"), c(1/2.5, 1/2.5), main="GB1 Fold (kT = 0.23)")
mtext("B", 3, .4, adj=-0.26, cex=12/10, font=2)

dev.off()

pdf("Figure5_2QMTGenContrib_Opt.pdf", FALSE, "Arial", width=2.2, height=2.2, pointsize=10)
#postscript("Figure5_2QMTGenContrib_Opt.eps", FALSE, "Arial", width=3.27, height=3.27, pointsize=10, horizontal=FALSE)
par(mar=c(2.7,2.7,1.5,0.2), mgp=c(1.7, 0.6, 0))
plot_gen_contrib(ent_gen_2QMT, attr(ent_gen_2QMT, "generations"), c(1/2.5, 1/2.5), 0.5515, main=paste("GB1 Fold (kT = ", kT_2QMT_label, ")", sep=""))
mtext("C", 3, .4, adj=-0.26, cex=12/10, font=2)
dev.off()

} # end if

} # end if

if (length(list.files("1N7T", "ga.entities")) == 2000) {


### Figure S1 ###

all_runs <- c("1A22_1", "1A22_2", "1A22_3", "1A22_4", "1A22_5", "1A22_6", "1N7T", "1N7T_V83K", "2FNE_A_C_V2048", "2I0L_A_C_V2006", "2IWP_B_A_V1927", "2QMT")

sizes <- c(200, 100, 50, 20)

freqhist <- numeric(20)

freqmaster <- rankmaster <- vector("list", length(all_runs))

multfreqstats <- rep(list(rep(list(numeric()), 20)), 4)
multrankstats <- rep(list(rep(list(numeric()), 20)), 4)

for (i in seq_along(all_runs)) {
	
	run <- all_runs[i]

	print(run)
	fitness_coef <- c(1/2.5, 1/2.5, 1/2.5, 1)
	names(fitness_coef) <- paste("state1_fitness_comp", seq_along(fitness_coef), sep="")
	if (run == "2QMT") fitness_coef <- fitness_coef[1:2]
	temp_or_thresh <- 0.228
	if (run == "1N7T_V83K") temp_or_thresh <- temp_or_thresh + 0.021
	type <- "boltzmann"
	percentile <- 0.5
	
	entities <- read_ga_entities_list(run, savedata=TRUE)
	pwms <- entities_list_pwms(entities, fitness_coef, temp_or_thresh, type)
	pwm <- collapse_pwms(pwms, percentile)
	pwmranks <- apply(-pwm, 2, rank)
	
	freqmaster[[i]] <- pwm
	rankmaster[[i]] <- apply(-pwm, 2, rank)
	
	emptymat <- pwm
	emptymat[seq_along(emptymat)] <- 0
	
	freqhist <- freqhist + hist(pwm, plot=FALSE, breaks=seq(0,1,by=.05))$counts
	
	for (j in 1:4) {
	
		size <- sizes[j]
	
		frmse <- emptymat
		rrmse <- emptymat
		
		print(size)
		for (start_idx in seq(1,2000,by=size)) {
			stop_idx <- start_idx+size-1
			pwm <- collapse_pwms(pwms[,start_idx:stop_idx,], percentile)
			rnk <- apply(-pwm, 2, rank)
			ifrmse <- (freqmaster[[i]]-pwm)^2
			irrmse <- (rankmaster[[i]]-rnk)^2
			frmse <- frmse + ifrmse
			rrmse <- rrmse + irrmse
			
			freqidx <- pmax(ceiling(pwm*20), 1)
			rankidx <- ceiling(rnk)
			
			for (idx in 1:20) {
		
				multfreqstats[[j]][[idx]] <- c(multfreqstats[[j]][[idx]], ifrmse[freqidx == idx])
				multrankstats[[j]][[idx]] <- c(multrankstats[[j]][[idx]], irrmse[rankidx == idx])
			}
		}
	}
}

freqstatsmat <- sqrt(sapply(multfreqstats, sapply, mean))
rankstatsmat <- sqrt(sapply(multrankstats, sapply, mean))

cols <- c("red", "orange", "cyan", "purple")

pdf("FigureS1_FreqRMSE.pdf", FALSE, "Arial", width=3.27, height=3, pointsize=10)
#postscript("FigureS1_FreqRMSE.eps", FALSE, "Arial", width=3.27, height=3.27, pointsize=10, horizontal=FALSE)

par(mar=c(2.7,2.7,0.2,0.2), mgp=c(1.7, 0.6, 0))

plot(0, 0, type="n", xlim=c(0,1), ylim=range(0,freqstatsmat), xlab="Predicted Frequency", ylab="Predicted Frequency RMSE")
for (i in seq_along(sizes)) {
	points(seq(0.025,0.975,by=.05), freqstatsmat[,i], type="l", col=cols[i])
}
exfreq <- seq(0.025,0.975,by=.05)[9]
exerror <- freqstatsmat[9,2]
segments(-1, exerror, exfreq, exerror, lty="dashed")
segments(exfreq, -1, exfreq, exerror, lty="dashed")
cat("Example Frequency: ", exfreq, " Error: ", exerror, "\n", sep="")

dev.off()

pdf("FigureS1_RankRMSE.pdf", FALSE, "Arial", width=3.27, height=3, pointsize=10)
#postscript("FigureS1_RankRMSE.eps", FALSE, "Arial", width=3.27, height=3.27, pointsize=10, horizontal=FALSE)

par(mar=c(2.7,2.7,0.2,0.2), mgp=c(1.7, 0.6, 0))

plot(0, 0, type="n", xlim=c(0,20), ylim=range(0,rankstatsmat+.02), xlab="Predicted Rank", ylab="Predicted Rank RMSE")
for (i in seq_along(sizes)) {
	points(seq(1,20), rankstatsmat[,i], type="l", col=cols[i])
}
exrank <- 3
exerror <- rankstatsmat[3,4]
segments(-1, exerror, exrank, exerror, lty="dashed")
segments(exrank, -1, exrank, exerror, lty="dashed")
cat("Example Rank: ", exrank, " Error: ", exerror, "\n", sep="")

dev.off()

num_backbones <- 2^seq(0,7)
num_backbones <- seq(1,40)
offsets <- seq(0, 1999, by=40)


### Figure S2 ###

bb_pwm_list <- vector("list", length(all_runs))
names(bb_pwm_list) <- all_runs

for (i in seq_along(all_runs)) {

	run <- all_runs[i]

	print(run)
	fitness_coef <- c(1/2.5, 1/2.5, 1/2.5, 1)
	names(fitness_coef) <- paste("state1_fitness_comp", seq_along(fitness_coef), sep="")
	if (run == "2QMT") fitness_coef <- fitness_coef[1:2]
	temp_or_thresh <- 0.228
	if (run == "1N7T_V83K") temp_or_thresh <- temp_or_thresh + 0.021
	type <- "boltzmann"
	percentile <- 0.5
	
	entities <- read_ga_entities_list(run, savedata=TRUE)
	pwms <- entities_list_pwms(entities, fitness_coef, temp_or_thresh, type)
	
	backbone_data <- vector("list", length(num_backbones))
	names(backbone_data) <- num_backbones
	
	for (j in seq_along(num_backbones)) {
	
		nbbs <- num_backbones[j]
		
		print(nbbs)
		
		rep_data <- vector("list", length(offsets))
		
		for (k in seq_along(offsets)) {
				
			rep_data[[k]] <- collapse_pwms(pwms[,seq_len(nbbs)+offsets[k],,drop=FALSE], percentile)
		}
		
		backbone_data[[j]] <- rep_data
	}
	
	bb_pwm_list[[i]] <- backbone_data
}

seqpos <- c(5,7,16,18,30,33)
wtseq <- strsplit("PNRFRGKDLAGSPGGTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTEGGGSSGGGADFDYEKMANAN", "")[[1]]
round3seq <- readLines("../data/gb1_round_3_seq.txt")

round3mat <- strsplit(round3seq, "")
round3mat <- do.call(rbind, round3mat[sapply(round3mat, length) == 90])

round3good <- rep(TRUE, nrow(round3mat))

for (i in seq_along(round3good)) {
	if (any(round3mat[i,seqpos+14] == "X")) round3good[i] <- FALSE
	if (any(round3mat[i,-(seqpos+14)] != wtseq[-(seqpos+14)])) round3good[i] <- FALSE
}

round3mat <- round3mat[round3good,seqpos+14]

exp_freq_2QMT <- exp_freq_mat(round3mat)
colnames(exp_freq_2QMT) <- seqpos

hgh_bind_freq <- t(as.matrix(read.table("../data/hgh_bind_freq.txt", TRUE, row=1)))/100
hgh_wt_seq <- gsub("[0-9]", "", colnames(hgh_bind_freq))
colnames(hgh_bind_freq) <- gsub("[A-Z]", "", colnames(hgh_bind_freq))

pdz_runs <- c("2I0L_A_C_V2006", "2IWP_B_A_V1927", "2FNE_A_C_V2048", "1N7T", "1N7T_V83K")
pdz_expfreq <- vector("list", length(pdz_runs))
names(pdz_expfreq) <- pdz_runs
for (i in seq_along(pdz_expfreq)) {
	pdz_expfreq[[i]] <- as.matrix(read.table(paste("../data/exp_pwm/", names(pdz_expfreq)[[i]], ".txt", sep=""), check.names=FALSE))
}

group_expfreq_list <- list(
	"GB1 (kT=0.23)"=list(GB1=exp_freq_2QMT),
	"hGH/hGHR (kT=0.23)"=list(hGH=hgh_bind_freq),
	"PDZ/Peptide"=pdz_expfreq
)

auc_mat <- matrix(ncol=length(num_backbones), nrow=length(offsets))
colnames(auc_mat) <- num_backbones

for (j in seq_along(num_backbones)) {
	
	nbbs <- num_backbones[j]
	
	print(nbbs)
	
	for (k in seq_along(offsets)) {
		
		hgh_freq <- vector("list", 6)
		for (i in 1:6) {
			dname <- paste("1A22_", i, sep="")
			hgh_freq[[i]] <- bb_pwm_list[[dname]][[j]][[k]]
			colnames(hgh_freq[[i]]) <- gsub("[^0-9]", "", colnames(hgh_freq[[i]]))
		}
		hgh_freq <- as.matrix(do.call(cbind, hgh_freq))
		hgh_freq <- hgh_freq[,order(as.numeric(colnames(hgh_freq)))]
		
		pdz_freq <- vector("list", length(pdz_runs))
		names(pdz_freq) <- pdz_runs
		for (i in seq_along(pdz_runs)) {
			pdz_freq[[i]] <- bb_pwm_list[[pdz_runs[i]]][[j]][[k]]
		}
		
		group_predfreq_list <- list(
			"GB1 (kT=0.23)"=list(GB1=bb_pwm_list[["2QMT"]][[j]][[k]]),
			"hGH/hGHR (kT=0.23)"=list(hGH=hgh_freq),
			"PDZ/Peptide"=pdz_freq
		)
		
		stats_tab <- stats_table(group_expfreq_list, group_predfreq_list)
		
		auc_mat[k,j] <- mean(stats_tab[,"auc"])
	}
}

postscript("FigureS2_AUCBox.eps", FALSE, "Arial", width=6.83, height=4, pointsize=10, horizontal=FALSE, paper="special")
par(mar=c(2.7,2.7,0.2,0.2), mgp=c(1.7, 0.6, 0))
boxplot(as.data.frame(auc_mat), xlab="Number of Backbones", ylab="Area Under ROC Curve", xaxt="n")
axis(1, at=seq(1, 40), labels=FALSE)
axis(1, at=seq(2, 40, by=2), labels=TRUE, tick=FALSE)
dev.off()

} else {


### Figure S3 ###

if (! "ent_1A22" %in% ls()) {
	ent_1A22 <- vector("list", 6)
	for (i in 1:6) ent_1A22[[i]] <- read_ga_entities_list(paste("1A22_", i, sep=""))
}

hgh_bind_freq <- t(as.matrix(read.table("../data/hgh_bind_freq.txt", TRUE, row=1)))/100
hgh_bind_freq_norm <- apply(hgh_bind_freq, 2, function(x) x/sum(x))
target_bits <- mean(apply(hgh_bind_freq_norm, 2, numbits))

param_mat <- as.matrix(expand.grid(
	intraweight=seq(0.2, 1.0, by=.1),
	percentile=seq(0.25, 0.95, by=.05)
))

if (! "result_mat" %in% ls()) result_mat <- matrix(nrow=0, ncol=9, dimnames=list(NULL, c("intraweight", "percentile", "kT", "bits_pred", "frac_top5",  "aad", "auc", "frobenius", "ranktop")))

for (i in seq_len(nrow(param_mat))) {

	intraweight <- unname(param_mat[i,1])
	percentile <- param_mat[i,2]

	idx <- which(abs(result_mat[,1]-intraweight) < 0.0001 & abs(result_mat[,2]-percentile) < 0.0001)
	
	best_objective <- Inf
	best_freq <- NULL
	
	if (length(idx)) next
	
	optfunc_1A22_bits <- function(x) {
		kT <- x
		pred_freq <- vector("list", 6)
		for (i in 1:6) pred_freq[[i]] <- collapse_pwms(entities_list_pwms(ent_1A22[[i]], c(intraweight, intraweight, intraweight, 1), kT), percentile)
		pred_freq <- as.matrix(do.call(cbind, pred_freq))
		pred_freq <- pred_freq[,order(as.numeric(colnames(pred_freq)))]
		
		nbits <- mean(apply(pred_freq, 2, numbits))
		objective <- (target_bits-nbits)^2
		if (objective < best_objective) {
			best_objective <<- objective
			best_freq <<- pred_freq
		}
		cat("kT=", kT, " numbits=", nbits, " objective=", objective, "\n", sep="")
		objective
	}
	optdata <- optimize(optfunc_1A22_bits, c(0.2,3))
	aad <- mean(abs(best_freq-hgh_bind_freq))
	nbits <- mean(apply(best_freq, 2, numbits))
	kT <- optdata$minimum
	print(c(optdata$objective, best_objective))
	result_mat <- rbind(result_mat, c(intraweight, percentile, kT, as.numeric(stats_table(list(list(hGH=hgh_bind_freq)), list(list(hGH=best_freq)))[,4:9])))
	
	cat("intraweight=", intraweight, " kT=", kT, " percentile=", percentile, " AAD=", aad*100, "%\n", sep="")
}

plot_result_mat <- result_mat[result_mat[,1]>=min(param_mat[,1]) & result_mat[,2]>=min(param_mat[,2]) & result_mat[,2]<=max(param_mat[,2]),]

kT_mat <- dftomat(plot_result_mat, "intraweight", "percentile", "kT")
aad_mat <- dftomat(plot_result_mat, "intraweight", "percentile", "aad")

pdf("FigureS3_kTArray.pdf", FALSE, "Arial", width=3.27, height=3, pointsize=10)
par(mar=c(2.7,2.7,1.5,0.2), mgp=c(1.7, 0.6, 0))
plot.scorearray(kT_mat, xlab="Intramolecular Weight", ylab="Percentile", zlab="", main="kT", levels=pretty(range(kT_mat), 8), lowpoint=FALSE)
dev.off()

pdf("FigureS3_AADArray.pdf", FALSE, "Arial", width=3.27, height=3, pointsize=10)
par(mar=c(2.7,2.7,1.5,0.2), mgp=c(1.7, 0.6, 0))
plot.scorearray(aad_mat, xlab="Intramolecular Weight", ylab="Percentile", zlab="", main="AAD (%)", levels=pretty(range(aad_mat), 8), lowpoint=FALSE)
dev.off()

intraweight <- 0.4
percentile <- 0.5
kT_1A22 <- kT_mat[as.character(intraweight),as.character(percentile)]
kT_1A22_label <- format(kT_1A22, digits=2)
hgh_freq_opt <- vector("list", 6)
for (i in 1:6) hgh_freq_opt[[i]] <- collapse_pwms(entities_list_pwms(ent_1A22[[i]], c(intraweight, intraweight, intraweight, 1), kT_1A22), percentile)
hgh_freq_opt <- as.matrix(do.call(cbind, hgh_freq_opt))
hgh_freq_opt <- hgh_freq_opt[,order(as.numeric(colnames(hgh_freq_opt)))]

print(c(aad_mat[as.character(intraweight),as.character(percentile)], min(aad_mat), (aad_mat[as.character(intraweight),as.character(percentile)] - min(aad_mat))/10))

### Table 1 ###

pdz_expfreq <- vector("list", length(pdz_freq))
names(pdz_expfreq) <- pdz_runs
pdz_freq_jmb2010 <- vector("list", length(pdz_freq))
names(pdz_freq_jmb2010) <- pdz_runs
for (i in seq_along(pdz_expfreq)) {
	pdz_expfreq[[i]] <- as.matrix(read.table(paste("../data/exp_pwm/", names(pdz_expfreq)[[i]], ".txt", sep=""), check.names=FALSE))
	pdz_freq_jmb2010[[i]] <- as.matrix(read.table(paste("../data/jmb2010_pwm/", names(pdz_expfreq)[[i]], ".txt", sep=""), check.names=FALSE))
}

group_expfreq_list <- list(
	"GB1 (kT=0.23)"=list(GB1=exp_freq_2QMT),
	"GB1 (kT=GB1kT)"=list(GB1=exp_freq_2QMT),
	"hGH/hGHR Table"=list(hGH=hgh_bind_freq[,hgh_table_pos]),
	"hGH/hGHR All"=list(hGH=hgh_bind_freq),
	"PDZ/Peptide"=pdz_expfreq,
	"PDZ/Peptide JMB2010"=pdz_expfreq
)

group_predfreq_list <- list(
	"GB1 (kT=0.23)"=list(GB1=pred_freq_2QMT),
	"GB1 (kT=GB1kT)"=list(GB1=pred_freq_2QMT_opt),
	"hGH/hGHR Table"=list(hGH=hgh_freq[,hgh_table_pos]),
	"hGH/hGHR All"=list(hGH=hgh_freq),
	"PDZ/Peptide"=pdz_freq,
	"PDZ/Peptide JMB2010"=pdz_freq_jmb2010
)

names(group_expfreq_list) <- names(group_predfreq_list) <- sub("GB1kT", kT_2QMT_label, names(group_expfreq_list))
names(group_expfreq_list) <- names(group_predfreq_list) <- sub("hGHkT", kT_1A22_label, names(group_expfreq_list))

stats_tab <- stats_table(group_expfreq_list, group_predfreq_list)

write.table(stats_tab, "Table1.txt", quote=FALSE, sep="\t", col.names=NA)


### Table S2 ###

gb1_freq_naive <- naive_pwm(wtseq[seqpos+14])
hgh_freq_naive <- naive_pwm(hgh_wt_seq)
pdz_freq_naive <- lapply(lapply(names(pdz_freq), start_sequence), naive_pwm)
names(pdz_freq_naive) <- names(pdz_freq)

group_expnaive_list <- list(
	"GB1"=list(GB1=exp_freq_2QMT),
	"hGH/hGHR Table"=list(hGH=hgh_bind_freq[,hgh_table_pos]),
	"hGH/hGHR"=list(hGH=hgh_bind_freq),
	"PDZ/Peptide"=pdz_expfreq
)

group_prednaive_list <- list(
	"GB1"=list(GB1=gb1_freq_naive),
	"hGH/hGHR Table"=list(hGH=hgh_freq_naive[,hgh_table_pos]),
	"hGH/hGHR"=list(hGH=hgh_freq_naive),
	"PDZ/Peptide"=pdz_freq_naive
)

naive_stats_tab <- stats_table(group_expnaive_list, group_prednaive_list)

write.table(naive_stats_tab, "TableS2.txt", quote=FALSE, sep="\t", col.names=NA)

} # end if
