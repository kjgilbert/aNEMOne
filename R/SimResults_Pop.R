#'
#' Takes an output file from Nemo of deleterious loci genotypes, and plots the distribution of both the homozygous and heterozygous effects of these loci.
#'
#' @title Examine some population parameters from Nemo stat output
#'
#'  @param input.stats.file The stat file output from Nemo.
#'
#'  @param which.rep If there are stats from multiple replicates, state which one to examine here, otherwise leave as NULL.
#'
#'  @param par.size Create the layout/number of plotting windows for making the desired output plots.
#'
#'  @param stats.to.plot List the stats to be plotted.
#'
#'  @return
#'
#'  Creates four plots of simulation results, showing mean fitness, total number of adults in the metapopulation, mean phenotypic value, and variances across generations.
#'
#' @author Kimberly J Gilbert
#'
#' @references \href{http://nemo2.sourceforge.net/index.html}{Nemo} is created and maintained by Fred Guillaume. The manual and source files are available online.
#'
#' @export sim.results.pop


sim.results.pop <- function(input.stats.file, which.rep=NULL, par.size=c(2,2), stats.to.plot=c("fitness.mean", "adlt.nbr", "adlt.q1", "adlt.q1.Va")){

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



#' Takes an "ind_seln" output file (with quanti and delet fitness values per individual) and plots fitness over the landscape in several ways.
#'
#' @title Look at landscape fitness from .fit files
#'
#'  @param input.fit.file The .fit file output from Nemo.
#'
#'  @param patches.x The number of patches on the landscape along the x-axis (parallel to expansion)
#'
#'  @param patches.y The number of patches on the landscape along the y-axis (perpendicular to expansion)
#'
#'  @param fitness.type Which component of fitness to plot: for the quantitative trait, for deleterious mutations, or total fitness.
#'
#'  @return
#'
#'  Creates a plot over the landscape (heat map style) for fitness of the specified component.
#'
#' @author Kimberly J Gilbert
#'
#' @references \href{http://nemo2.sourceforge.net/index.html}{Nemo} is created and maintained by Fred Guillaume. The manual and source files are available online.
#'
#' @export fit.results.landscape


fit.results.landscape <- function(input.fit.file, patches.x, patches.y, fitness.type="total"){

	input <- read.table(input.fit.file, header=TRUE)
	# column 1 is pop
	# column 2 is fitness for trait 1 = delet
	# column 3 is fitness for trait 2 = quanti
	input$total.fitness <- input$trait1 * input$trait2
  
	# average within patches
	per.patch.fitness <- aggregate(input, by=list(input$pop), FUN=mean)
	total.num.patches <- patches.x*patches.y
		# there shouldn't be any ghost patches in the list because they're culled to pop size zero, but if errors arise down the line, that could be at fault
	if(dim(per.patch.fitness)[1] > total.num.patches) per.patch.fitness <- per.patch.fitness[- (total.num.patches+1),]
	if(dim(per.patch.fitness)[1] > total.num.patches) per.patch.fitness <- per.patch.fitness[- (total.num.patches+1),]
		# if some patches are empty:
	if(dim(per.patch.fitness)[1] < total.num.patches){
		empty.patches <- setdiff(1:total.num.patches, per.patch.fitness$pop)
		empties <- matrix(0, ncol=5, nrow=length(empty.patches))
		empty.rows <- data.frame(matrix(cbind(empty.patches, empty.patches, empties), ncol=7, byrow=TRUE))
		names(empty.rows) <- names(per.patch.fitness)
		per.patch.fitness <- rbind(per.patch.fitness, empty.rows)
		per.patch.fitness <- per.patch.fitness[order(per.patch.fitness$pop), ]	
	}
	  
	if(fitness.type == "total"){
		# make total fitness into a matrix matched to the landscape
		total.fit.mat <- matrix(per.patch.fitness$total.fitness, nrow=patches.y, ncol=patches.x, byrow=FALSE)
		heatmap(x=total.fit.mat, col=heat.colors(256), scale="column", margins=c(2,2), xlab="x axis", ylab="y axis", Rowv=NA, Colv=NA, labRow=NA, labCol=NA, main="Mean total fitness per patch")
	}
	if(fitness.type == "quanti"){
		# make quanti fitness into a matrix matched to the landscape
		quanti.fit.mat <- matrix(per.patch.fitness$trait2, nrow=patches.y, ncol=patches.x, byrow=FALSE)
		heatmap(x=quanti.fit.mat, col=heat.colors(256), scale="column", margins=c(2,2), xlab="x axis", ylab="y axis", Rowv=NA, Colv=NA, labRow=NA, labCol=NA, main="Mean quanti fitness per patch")
	}
	if(fitness.type == "delet"){
		# make delet fitness into a matrix matched to the landscape
		delet.fit.mat <- matrix(per.patch.fitness$trait1, nrow=patches.y, ncol=patches.x, byrow=FALSE)
		heatmap(x=delet.fit.mat, col=heat.colors(256), scale="column", margins=c(2,2), xlab="x axis", ylab="y axis", Rowv=NA, Colv=NA, labRow=NA, labCol=NA, main="Mean delet fitness per patch")
	}
}