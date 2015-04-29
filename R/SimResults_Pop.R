#'
#' Takes an output file from Nemo of deleterious loci genotypes, and plots the distribution of both the homozygous and heterozygous effects of these loci.
#'
#' @title Examine some population parameters from Nemo stat output
#'
#'
#'  @param input The stat file output from Nemo.
#'
#'  @param out.file The name of the file to output the .pdf of figures to, currently not enabled in the function.
#'
#'  @param which.rep If there are stats from multiple replicates, state which one to examine here, otherwise leave as NULL.
#'
#'  @param add If adding to an existing plor, change to TRUE, otherwise leave at default of FALSE I think this will probably not work as is.
#'
#'  @return
#'
#'  Creates four plots of simulation results, showing mean fitness, total number of adults in the metapopulation, mean phenotypic value, and variances across generations.
#'
#' @author Kimberly J Gilbert
#'
#' @references \href{http://nemo2.sourceforge.net/index.html}{Nemo} is created and maintained by Fred Guillaume. The manual and source files are available online.
#'
#' @export


sim.results.pop <- function(input, out.file="SimResults.pdf", which.rep=NULL, add=FALSE){

	# separate replicates
	if(is.null(which.rep)){
		dat <- input
	}else{
		dat <- subset(input, input$replicate == which.rep)
	}


#	pdf(out.file, width=10, height=10)

	par(mfrow=c(2,2))


	# MEAN FITNESS

	plot(dat$generation, dat$fitness.mean, ylim=c(0,1), type="o", xlab="Generation", ylab="Mean Fitness", col="blue", add=add)


	# TOTAL NUMBER OF ADULTS IN THE METAPOP

	plot(dat$generation, dat$adlt.nbr, type="o", col="blue", xlab="Generation", ylab="Total number adults in metapop", add=add)


	# GENETIC STATs

	plot(dat$generation, dat$adlt.q1, ylim=range(dat$adlt.q1), xlab="Generation", ylab="Mean phenotypic value", type="o", col="blue", add=add)		# mean phenotypic value of the trait in the whole metapop
	points(dat$generation, dat$off.q1, col="lightblue", type="o", pch=2)
	legend("topright", c("Adult", "Offspring"), col=c("blue", "lightblue"), pch=c(1,2), pt.lwd=2)


	# Va, Vb, Vp

	plot(dat$generation, dat$adlt.q1.Va, type="o", ylim=c(0,8), col="black", pch=0, xlab="Generation", ylab="", add=add)	# avg w/n patch additive gen. var.
	points(dat$generation, dat$adlt.q1.Vb, type="o", col="blue", pch=1) # among patch gen. var.
	points(dat$generation, dat$adlt.q1.Vp, type="o", col="red", pch=5) # avg w/n patch pheno. var., only if Ve != 0
	#points(dat$generation, dat$off.q1.Va, type="o", col="black", pch=2, cex=0.5)
	#points(dat$generation, dat$off.q1.Vb, type="o", col="red", pch=2, cex=0.5)
	#points(dat$generation, dat$off.q1.Vp, type="o", col="green", pch=2, cex=0.5)
	legend("topright", c("Va", "Vb", "Vp"), col=c("black", "blue", "red"), pch=c(0,1,5))


#	dev.off()

}
