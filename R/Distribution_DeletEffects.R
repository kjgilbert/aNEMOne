#'
#' Takes an output file from Nemo of deleterious loci genotypes, and plots the distribution of both the homozygous and heterozygous effects of these loci.
#'
#' @title Examine effect size distribution of deleterious loci from Nemo
#'
#'
#'  @param file The file containing deleterious loci output from Nemo.
#'
#'  @param num.loci The number of deleterious loci that were simulated in that file.
#'
#'  @param xlim.ho The limits of the x-axis for the homozygous effect distribution, default is from 0 to 1.
#'
#'  @param xlim.he The limits of the x-axis for the heterozygous effect distribution, default is from 0 to 1.
#'
#'  @return
#'
#'  Creates a plot of two histograms showing the distributions.
#'
#' @author Kimberly J Gilbert
#'
#' @references \href{http://nemo2.sourceforge.net/index.html}{Nemo} is created and maintained by Fred Guillaume. The manual and source files are available online.
#'
#' @export dist.delet.effects


dist.delet.effects <- function(file, num.loci, xlim.ho=c(0,1), xlim.he=c(0,1)){
#	delet.traits <- read.table(file, header=TRUE, sep=" ", stringsAsFactors=FALSE)
  delet.traits <- matrix(scan(file, skip=1, nlines=2, what="numeric()"), ncol=num.loci+5, byrow=TRUE)
  # strip the -1's from nemo's extra information
  # because there are 1000 loci, get rid of spot 1, pop ID and last 4 spots - age, sex, ped, origin
  delet.traits <- delet.traits[, -c(1,(num.loci+2):(num.loci+5))]
  
  # delet loci effect sizes:
	ho <- as.numeric(delet.traits[1,])
	he <- as.numeric(delet.traits[2,])
	
	par(mfrow=c(1,2))
	hist(as.matrix(ho), col="steelblue1", breaks=50, xlab="Homozygous Effect Size", main="", xlim=xlim.ho)
	hist(as.matrix(he), col="steelblue3", breaks=50, xlab="Heterozygous Effect Size", main="", xlim=xlim.he)

}



#'
#' Takes an output file from Nemo of deleterious loci genotypes, and plots the mean number of deleterious mutations per individual in a patch over the landscape.
#'
#' @title Look at numbers of deleterious mutations across the landscape
#'
#'  @param del.file The file containing deleterious loci output from Nemo.
#'
#'  @param num.loci The number of deleterious loci that were simulated in that file.
#'
#'  @param patches.x The number of patches on the landscape along the x-axis (parallel to expansion).
#'
#'  @param patches.y The number of patches on the landscape along the y-axis (perpendicular to expansion).
#'
#'  @param count.type Whether to count homozygous, heterozygous, or total number of mutations.
#'
#'  @return
#'
#'  Creates a plot over the landscape (heat map style) for the mean number of deleterious mutations per individual in a patch.
#'
#' @author Kimberly J Gilbert
#'
#' @references \href{http://nemo2.sourceforge.net/index.html}{Nemo} is created and maintained by Fred Guillaume. The manual and source files are available online.
#'
#' @export delet.muts.over.landscape


delet.muts.over.landscape <- function(del.file, patches.x, patches.y, num.loci, count.type="total"){

	# custom colors
	purple=rgb(1,0,1)
	red=rgb(1,0,0)
	yellow=rgb(1,1,0)
	green=rgb(0,1,0)
	teal=rgb(0,1,1)
	blue=rgb(0,0,1)
	white=rgb(1,1,1)
	whiteToWhite <- colorRampPalette(white)
	whiteToYellow <- colorRampPalette(c(white, yellow))
	yellowToRed <- colorRampPalette(c(yellow, red))
	redToPurple <- colorRampPalette(c(red, purple))
	greenToWhite <- colorRampPalette(c(green, white))

	delet.muts <- matrix(scan(del.file, skip=3, what="character()"), ncol=num.loci+6, byrow=TRUE)
	# strip the -1's from nemo's extra information
	# because there are 1000 loci, last 5 spots - age, sex, ped, origin, and some other number
	delet.muts <- delet.muts[, -c((num.loci+2):(num.loci+6))]

	pop.list <- as.numeric(delet.muts[,1])
	delet.muts <- data.frame(delet.muts[,-1])
	
	num.zero <- apply(delet.muts, MARGIN=1, FUN=function(x) length(which(x == "00")))
	num.hets <- apply(delet.muts, MARGIN=1, FUN=function(x) length(which(x == "01" | x=="10")))
	num.homs <- apply(delet.muts, MARGIN=1, FUN=function(x) length(which(x == "11")))
	total.muts <- num.hets + (2*num.homs)

	mut.counts <- data.frame(cbind(pop.list, num.zero, num.hets, num.homs, total.muts))
	avg.mut.counts <- aggregate(mut.counts, by=list(mut.counts$pop.list), FUN=mean)

	total.num.patches <- patches.x*patches.y
		# there shouldn't be any ghost patches in the list because they're culled to pop size zero, but if errors arise down the line, that could be at fault
	if(dim(avg.mut.counts)[1] > total.num.patches) avg.mut.counts <- avg.mut.counts[- (total.num.patches+1),]
	if(dim(avg.mut.counts)[1] > total.num.patches) avg.mut.counts <- avg.mut.counts[- (total.num.patches+1),]
		# if some patches are empty:
	if(dim(avg.mut.counts)[1] < total.num.patches){
		empty.patches <- setdiff(1:total.num.patches, avg.mut.counts$pop.list)
		empty.rows <- data.frame(matrix(NA, ncol=6, nrow=length(empty.patches)))
		empty.rows[,1] <- empty.patches
		empty.rows[,2] <- empty.patches
		names(empty.rows) <- names(avg.mut.counts)
		avg.mut.counts <- rbind(avg.mut.counts, empty.rows)
		avg.mut.counts <- avg.mut.counts[order(avg.mut.counts$pop.list), ]	
	}
	  
	if(count.type == "total" | count.type == "all"){
		# make total fitness into a matrix matched to the landscape
		total.mut.mat <- matrix(avg.mut.counts$total.muts, nrow=patches.y, ncol=patches.x, byrow=FALSE)
		heatmap.2(x=total.mut.mat, dendrogram='none', labRow=NA, labCol=NA, margins=c(2,1), trace='none', na.color="blue", keysize=1, key.ylab=NA, key.title=NA, key.par=list(yaxt="n"), xlab="axis of expansion ->", col=c(whiteToWhite(40), whiteToYellow(60), yellowToRed(60), redToPurple(15)), main="Mean total number delet muts per ind (within a patch)")
	}
	if(count.type == "homozygous" | count.type == "all"){
		# make quanti fitness into a matrix matched to the landscape
		hom.mut.mat <- matrix(avg.mut.counts$num.homs, nrow=patches.y, ncol=patches.x, byrow=FALSE)
		heatmap.2(x=hom.mut.mat, dendrogram='none', labRow=NA, labCol=NA, margins=c(2,1), trace='none', na.color="blue", keysize=1, key.ylab=NA, key.title=NA, key.par=list(yaxt="n"), xlab="axis of expansion ->", col=c(whiteToWhite(40), whiteToYellow(60), yellowToRed(60), redToPurple(15)), main="Mean number homozygous delet muts per ind (within a patch)")
	}
	if(count.type == "heterozygous" | count.type == "all"){
		# make delet fitness into a matrix matched to the landscape
		het.mut.mat <- matrix(avg.mut.counts$num.hets, nrow=patches.y, ncol=patches.x, byrow=FALSE)
		heatmap.2(x=het.mut.mat, dendrogram='none', labRow=NA, labCol=NA, margins=c(2,1), trace='none', na.color="blue", keysize=1, key.ylab=NA, key.title=NA, key.par=list(yaxt="n"), xlab="axis of expansion ->", col=c(whiteToWhite(40), whiteToYellow(60), yellowToRed(60), redToPurple(15)), main="Mean number heterozygous delet muts per ind (within a patch)")
	}

}

