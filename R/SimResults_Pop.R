#'
#' Takes an output file from Nemo of deleterious loci genotypes, and plots the distribution of both the homozygous and heterozygous effects of these loci.
#'
#' @title Examine some population parameters from Nemo stat output
#'
#'
#'  @param input.stats.file The stat file output from Nemo.
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


sim.results.pop <- function(input.stats.file, out.file="SimResults.pdf", which.rep=NULL, par.size=c(2,2), stats.to.plot=c("fitness.mean", "adlt.nbr", "adlt.q1", "adlt.q1.Va")){

  input <- read.table(input.stats.file, header=TRUE)
  
	# separate replicates
	if(is.null(which.rep)){
		dat <- input
	}else{
		dat <- subset(input, input$replicate == which.rep)
	}

	par(mfrow=par.size)
	
	for(i in stats.to.plot){
	  plot(dat$generation, unlist(dat[i]), type="o", xlab="Generation", ylab=i, col="blue")
	}
}
