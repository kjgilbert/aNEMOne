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

	# custom colors
	purple=rgb(1,0,1)
	red=rgb(1,0,0)
	yellow=rgb(1,1,0)
	green=rgb(0,1,0)
	teal=rgb(0,1,1)
	blue=rgb(0,0,1)
	white=rgb(1,1,1)
	whiteToWhite <- colorRampPalette(c(white,white))
	redToYellow <- colorRampPalette(c(red,yellow))
	whiteToYellow <- colorRampPalette(c(white, yellow))
	yellowToTeal <- colorRampPalette(c(yellow, teal))
	greenToTeal <- colorRampPalette(c(green, teal))
	tealToBlue <- colorRampPalette(c(teal, blue))
	yellowToGreen <- colorRampPalette(c(yellow, green))


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
		empty.rows <- data.frame(matrix(0, ncol=7, nrow=length(empty.patches)))
		empty.rows[,1] <- empty.patches
		empty.rows[,2] <- empty.patches
		names(empty.rows) <- names(per.patch.fitness)
		per.patch.fitness <- rbind(per.patch.fitness, empty.rows)
		per.patch.fitness <- per.patch.fitness[order(per.patch.fitness$pop), ]	
	}
	  
	if(fitness.type == "total" | fitness.type == "all"){
		# make total fitness into a matrix matched to the landscape
		total.fit.mat <- matrix(per.patch.fitness$total.fitness, nrow=patches.y, ncol=patches.x, byrow=FALSE)
		# make the scale always the same:
		total.fit.mat[1,1] <- 1; total.fit.mat[2,1] <- 0
		image.plot(x=1:patches.x, y=1:patches.y, t(total.fit.mat), col=c(whiteToWhite(1), redToYellow(60), yellowToGreen(15), greenToTeal(15), tealToBlue(15)), ylab="", xlab="Axis of expansion", main="Mean total fitness per patch")
	}
	if(fitness.type == "quanti" | fitness.type == "all"){
		# make quanti fitness into a matrix matched to the landscape
		quanti.fit.mat <- matrix(per.patch.fitness$trait2, nrow=patches.y, ncol=patches.x, byrow=FALSE)
		# make the scale always the same:
		quanti.fit.mat[1,1] <- 1; quanti.fit.mat[2,1] <- 0
		image.plot(x=1:patches.x, y=1:patches.y, t(quanti.fit.mat), col=c(whiteToWhite(1), redToYellow(60), yellowToGreen(15), greenToTeal(15), tealToBlue(15)), ylab="", xlab="Axis of expansion", main="Mean quanti fitness per patch")
	}
	if(fitness.type == "delet" | fitness.type == "all"){
		# make delet fitness into a matrix matched to the landscape
		delet.fit.mat <- matrix(per.patch.fitness$trait1, nrow=patches.y, ncol=patches.x, byrow=FALSE)
		# make the scale always the same:
		delet.fit.mat[1,1] <- 1; delet.fit.mat[2,1] <- 0
		image.plot(x=1:patches.x, y=1:patches.y, t(delet.fit.mat), col=c(whiteToWhite(1), redToYellow(60), yellowToGreen(15), greenToTeal(15), tealToBlue(15)), ylab="", xlab="Axis of expansion", main="Mean delet fitness per patch")
	}
}








#' Makes 1-D plots over the landscape for population size, fitness, numbers of deleterious mutations, and environmental optimum match to the quantitative trait. Simulations must have been in 1 dimension for current implementation.
#'
#' @title Look at 1-D simulation results over landscape
#'
#'  @param input.fit.file The .fit file output from Nemo.
#'
#'  @param input.quanti.file The .quanti file output from Nemo which contains genotypes for the quantitative trait.
#'
#'  @param input.del.file The .del file output from Nemo which contains deleterious mutation genotypes.
#'
#'  @param patches.x The number of patches on the landscape along the x-axis (parallel to expansion)
#'
#'  @param patches.y The number of patches on the landscape along the y-axis (perpendicular to expansion)
#'
#'  @param slope.opt The slope of the environmental optimum defined in the input file.
#'
#'  @param opt.rate.change The rate of change of the environment per generation defined in the input file.
#'
#'  @param generation The current generation being plotted.
#'
#'  @param del.loci The number of deleterious loci simulated.
#'
#'  @param xlimits A smaller range on the x-axis to plot, if desired.
#'
#'  @return
#'
#'  Creates 1-D plots over the landscape for population size, fitness, numbers of deleterious mutations, and environmental optimum match to the quantitative trait.
#'
#' @author Kimberly J Gilbert
#'
#' @references \href{http://nemo2.sourceforge.net/index.html}{Nemo} is created and maintained by Fred Guillaume. The manual and source files are available online.
#'
#' @export plot1D.results


plot1D.results <- function(input.fit.file, input.quanti.file=NULL, input.del.file=NULL, patches.x, patches.y, slope.opt=NULL, opt.rate.change=NULL, generation, del.loci, xlimits=NULL){
  
  total.num.patches <- patches.x * patches.y
  if(is.null(xlimits)){
  	xlimits <- c(1, patches.x)
  }
  
  
  # make some pop size and fitness plots
    # column 1 is pop
    # column 2 is fitness for trait 1 = delet
    # column 3 is fitness for trait 2 = quanti
  fit.input <- read.table(input.fit.file, header=TRUE)
  fit.input$total.fitness <- fit.input$trait1 * fit.input$trait2
  
  # average within patches
  per.patch.fitness <- aggregate(fit.input, by=list(fit.input$pop), FUN=mean)

  pop.size <- aggregate(fit.input, by=list(fit.input$pop), FUN=length)
  # there shouldn't be any ghost patches in the list because they're culled to pop size zero, but if errors arise down the line, that could be at fault
  if(dim(per.patch.fitness)[1] > total.num.patches) per.patch.fitness <- per.patch.fitness[- (total.num.patches+1),]
  if(dim(per.patch.fitness)[1] > total.num.patches) per.patch.fitness <- per.patch.fitness[- (total.num.patches+1),]
  # if some patches are empty:
  if(dim(per.patch.fitness)[1] < total.num.patches){
    empty.patches <- setdiff(1:total.num.patches, per.patch.fitness$pop)
    empty.rows <- data.frame(matrix(NA, ncol=7, nrow=length(empty.patches)))
    empty.rows[,1] <- empty.patches
    empty.rows[,2] <- empty.patches
    names(empty.rows) <- names(per.patch.fitness)
    per.patch.fitness <- rbind(per.patch.fitness, empty.rows)
    per.patch.fitness <- per.patch.fitness[order(per.patch.fitness$pop), ]
    
    empty.rows[,2] <- rep(NA, length(empty.rows[,1]))
    per.patch.pop.size <- rbind(pop.size, empty.rows)
    per.patch.pop.size <- per.patch.pop.size[order(per.patch.pop.size$Group.1), ]
  }
  
  plot(1:total.num.patches, per.patch.pop.size$pop, type="o", col="darkorange", pch=".", lwd=2, ylim=c(0,100), xlim=xlimits, xlab="Landscape x position", ylab="Population size", main=paste(c("Generation ", generation), collapse=""))
  
  plot(per.patch.fitness$Group.1, per.patch.fitness$total.fitness, type="o", col="black", pch=".", lwd=2, ylim=c(0,1), xlim=xlimits, xlab="Landscape x position", ylab="Mean fitness", main=paste(c("Generation ", generation), collapse=""))
  points(per.patch.fitness$Group.1, per.patch.fitness$trait2, type="o", col="blue", pch=".", lwd=2, ylim=c(0,1), xlim=xlimits)
  points(per.patch.fitness$Group.1, per.patch.fitness$trait1, type="o", col="red", pch=".", lwd=2, ylim=c(0,1), xlim=xlimits)
  
  
  
  
  
  
  # make some deleterious mutation plots
  if(!is.null(input.del.file)){
    delet.input <- matrix(scan(input.del.file, skip=3, what="character()"), ncol=del.loci+6, byrow=TRUE)
    delet.muts <- delet.input[, -c((del.loci+2):(del.loci+6))]
    
    pop.list <- as.numeric(delet.muts[,1])
    delet.muts <- data.frame(delet.muts[,-1])
    
    num.zero <- apply(delet.muts, MARGIN=1, FUN=function(x) length(which(x == "00")))
    num.hets <- apply(delet.muts, MARGIN=1, FUN=function(x) length(which(x == "01" | x=="10")))
    num.homs <- apply(delet.muts, MARGIN=1, FUN=function(x) length(which(x == "11")))
    total.muts <- num.hets + (2*num.homs)
    
    mut.counts <- data.frame(cbind(pop.list, num.zero, num.hets, num.homs, total.muts))
    avg.mut.counts <- aggregate(mut.counts, by=list(mut.counts$pop.list), FUN=mean)
    
    if(dim(avg.mut.counts)[1] < total.num.patches){
      empty.patches <- setdiff(1:total.num.patches, avg.mut.counts$pop.list)
      empty.rows <- data.frame(matrix(NA, ncol=6, nrow=length(empty.patches)))
      empty.rows[,1] <- empty.patches
      empty.rows[,2] <- empty.patches
      names(empty.rows) <- names(avg.mut.counts)
      avg.mut.counts <- rbind(avg.mut.counts, empty.rows)
      avg.mut.counts <- avg.mut.counts[order(avg.mut.counts$pop.list), ]
    }
    
    plot(1:patches.x, avg.mut.counts$total.muts, xlim=xlimits, type="l", lwd=2, xlab="Landscape x position", ylab="Mean number deleterious mutations", main=paste(c("Generation ", generation), collapse=""))
    points(avg.mut.counts$pop.list, avg.mut.counts$num.hets, xlim=xlimits, type="l", lwd=2, col="orange")
    points(avg.mut.counts$pop.list, avg.mut.counts$num.homs, xlim=xlimits, type="l", lwd=2, col="red")
  }
 
  
   
  # make some environment and quanti geno/pheno plots
  if(!is.null(input.quanti.file) & !is.null(slope.opt)){
    quanti.input <- read.table(input.quanti.file, header=TRUE)
    # col G1 is geno, col P1 is pheno
    quanti.input <- quanti.input[,c("pop", "G1", "P1")]

    avg.quanti <- aggregate(quanti.input, by=list(quanti.input$pop), FUN=mean)
    
    if(dim(avg.quanti)[1] < total.num.patches){
      empty.patches <- setdiff(1:total.num.patches, avg.quanti$pop)
      empty.rows <- data.frame(matrix(NA, ncol=4, nrow=length(empty.patches)))
      empty.rows[,1] <- empty.patches
      empty.rows[,2] <- empty.patches
      names(empty.rows) <- names(avg.quanti)
      avg.quanti <- rbind(avg.quanti, empty.rows)
      avg.quanti <- avg.quanti[order(avg.quanti$pop), ]
    }
    
    # remake the landscape:
    first.col <- -(slope.opt*patches.x)/2
    last.col <- (slope.opt*patches.x)/2
    num.steps <- patches.x   # all of width 1
    step.values <- seq(first.col, last.col, length.out=num.steps)
    total.num.patches <- patches.x*patches.y
    patches.per.step <- total.num.patches/num.steps
    
    env <- NULL
    for(i in 1:num.steps){
      one.step <- rep(step.values[i], patches.per.step)
      env <- c(env, one.step)
    }
    
    # plot the landscape and genotypes and phenotypes on it
    env.change.time <- generation*opt.rate.change
    env <- env + env.change.time
    
    plot(1:patches.x, env, xlim=xlimits, type="l", lwd=1.5, xlab="Landscape x position", ylab="Quanti trait & env. optimum", main=paste(c("Generation ", generation), collapse=""))
    points(avg.quanti$pop, avg.quanti$G1, xlim=xlimits, type="l", lwd=2, col="blue")
    points(avg.quanti$pop, avg.quanti$P1, xlim=xlimits, type="l", lwd=2, col="green3")
  }
  
}