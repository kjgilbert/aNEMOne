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
#' @export


dist.delet.effects <- function(file, num.loci, xlim.ho=c(0,1), xlim.he=c(0,1)){
	delet.traits <- read.table(file, header=TRUE, sep=" ", stringsAsFactors=FALSE)
	# delet loci effect sizes:
	# because there are 1000 loci, get rid of spot 1, pop ID and last 4 spots - age, sex, ped, origin
	ho <- delet.traits[1,2:(num.loci+1)]
	he <- delet.traits[2,2:(num.loci+1)]
	
	par(mfrow=c(1,2))
	hist(as.matrix(ho), col="steelblue1", breaks=50, xlab="Homozygous Effect Size", main="", xlim=xlim.ho)
	hist(as.matrix(he), col="steelblue3", breaks=50, xlab="Heterozygous Effect Size", main="", xlim=xlim.he)

}