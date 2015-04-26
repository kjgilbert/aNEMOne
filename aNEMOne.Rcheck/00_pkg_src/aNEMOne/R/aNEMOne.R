# Functions to create Nemo input and analyse Nemo output 
#	by Kimberly J. Gilbert (kgilbert@zoology.ubc.ca)

# made to work with code edits by Kim for breeding window and connectivity matrix with dspersal kernel, see: https://github.com/kjgilbert/NemoDispersalKernel


#'
#' Compares the QST of a single phenotypic trait to the mean FST of series of marker loci. 
#' It calculates the distribution of QST - FST under a model assuming neutrality of both 
#' the phenotypic trait and the genetic markers from which FST is estimated.
#' Returns the simulated estimates of Qst - Fst under neutrality following the procedure
#' described in Gilbert and Whitlock (2014) and Whitlock & Guillaume (2009). Also returns 
#' the simulated estimates of Fst and Qst used to compute the null distribution.
#'
#' @title Create Nemo input file
#'
#'
#'  @param run.mode
#'
#'  @param random.seed Set the random seed for the run, default is 12345.
#'
#'  @param log.file Name of the logfile to be output.
#'
#'  @param root.dir
#'
#'  @param filename 
#'
#'  @param reps
#'
#'  @param gens
#'
#'  @param num.patches
#'
#'  @param patch.capacity
#'
#'  @param patch.capacity.fem
#'
#'  @param patch.capacity.mal
#'
#'  @param cap.temp
#'
#'  @param LCE.order
#'
#'  @param mating.system
#'
#'  @param mean.fec
#'
#'  @param seln.trait
#'
#'  @param seln.model
#'
#'  @return
#'
#'  Write the file specified into the working directory.
#'
#' @author Kimberly J Gilbert
#'
#'
#'
#' @examples
#'
#' # population
#' cap <- patch.cap(6,1000, 0,2)
#' num <- 1002
#' 
#' # simulation components
#' life.cycles <- c("breed", "disperse", "selection", "aging")
#' 
#' 
#' # mating
#' 
#' breed.kernel <- make.kernel.and.matrix(cell.size=50, land.x=1000, land.y=2000, dist.mean=0, dist.sd=25, breed.window=TRUE, two.kernels=FALSE, second.dist.mean=0, second.dist.sd=NULL)
#' 
#' 
#' # dispersal
#' 
#' disp.kernel <- make.kernel.and.matrix <- function(cell.size=50, land.x=1000, land.y=2000, dist.mean=0, dist.sd=25, breed.window=FALSE, two.kernels=FALSE, second.dist.mean=0, second.dist.sd=NULL)
#' 
#' 
#' make.input(root.dir="test", filename="test2", 
#' 	reps=10, gens=100, num.patches=num, patch.capacity=cap, 
#' 	LCE.order= life.cycles, mating.system=1, mean.fec=5
#' 	)
#'  
#' 
#' @export



make.input <- function(
	run.mode="overwrite",
	random.seed="12345",
	log.file="logfile.log",
	root.dir=NULL,
	filename=NULL,
	reps=NULL,
	gens=NULL,
	num.patches=NULL,
	patch.capacity=NULL,
	patch.capacity.fem=NULL,
	patch.capacity.mal=NULL,
	cap.temp=FALSE,
	LCE.order=NULL,
	mating.system=NULL,
	mean.fec=NULL,
	seln.trait=NULL,
	seln.model=NULL
	){
		row1 <- paste(c("run_mode", run.mode), collapse=" ")
		row2 <- paste(c("random_seed", random.seed), collapse=" ")
		row3 <- paste(c("logfile", log.file), collapse=" ")
		row4 <- paste(c("root_dir", root.dir), collapse=" ")
		row5 <- paste(c("filename", filename), collapse=" ")
		row6 <- paste(c("replicates", reps), collapse=" ")
		row7 <- paste(c("generations", gens), collapse=" ")
		row8 <- paste(c("patch_number", num.patches), collapse=" ")
		if(is.null(patch.capacity.fem)){
			row9 <- paste(c("patch_capacity", patch.capacity), collapse=" ")
		}else{
			row9 <- paste(c(paste(c("patch_nbfem", patch.capacity.fem), collapse=" "), paste(c("patch_nbmal", patch.capacity.mal), collapse=" ")), collapse="\n")
		}
		if(cap.temp==TRUE) row9 <- paste("(", row9, ")")
		
		row10 <- NULL
		for(i in 1:length(LCE.order)){
			temp <- paste(paste(c(LCE.order[i], i), collapse=" "), sep="\n")
			row10 <- paste(c(paste(row10), paste(temp), sep="\n"))
		}
		
		row11 <- paste(c("mating_system", mating.system), collapse=" ")
		row12 <- paste(c("mean_fecundity", mean.fec), collapse=" ")
		row13 <- paste(c("selection_trait", seln.trait), collapse=" ")
		row14 <- paste(c("selection_model", seln.model), collapse=" ")
		row15 <- paste(c("selection_fitness_model"), collapse=" ")
		row16 <- paste(c("selection_variance"), collapse=" ")
		row17 <- paste(c("selection_trait_dimension"), collapse=" ")
		row18 <- paste(c("selection_local_optima"), collapse=" ")
		row19 <- paste(c("quanti_init_trait_values"), collapse=" ")
		row20 <- paste(c("ntrl_init_patch_freq"), collapse=" ")
		row21 <- paste(c(), collapse=" ")
		row22 <- paste(c(), collapse=" ")
		row23 <- paste(c(), collapse=" ")
		row24 <- paste(c(), collapse=" ")
		row25 <- paste(c(), collapse=" ")
		row26 <- paste(c(), collapse=" ")
		row27 <- paste(c(), collapse=" ")
		row28 <- paste(c(), collapse=" ")
		row29 <- paste(c(), collapse=" ")
		row30 <- paste(c(), collapse=" ")
		row31 <- paste(c(), collapse=" ")
		row32 <- paste(c(), collapse=" ")
		row33 <- paste(c(), collapse=" ")
		row34 <- paste(c(), collapse=" ")
		row35 <- paste(c(), collapse=" ")
		row36 <- paste(c(), collapse=" ")
		row37 <- paste(c(), collapse=" ")
		row38 <- paste(c(), collapse=" ")
		row39 <- paste(c(), collapse=" ")
		row40 <- paste(c(), collapse=" ")
		
		
	
	init.file <- paste(c(row1, row2, row3, row4, row5, "\n", row6, row7, row8, row9, row10, row11, row12), collapse="\n")
	
		
	writeLines(init.file, paste(c(getwd(), "/", filename, ".ini"), collapse=""))
		
	return(paste(c("File written to ", getwd(), "/", filename, ".ini"), collapse=""))
}




#'
#' @title Analyze a Nemo stat output file
#'
#' @description Analyze output from Nemo specified by the stat options from a Nemo .init file 
#'
#'
#'  @param filename The name and path to the file to be analyzed.
#'
#'  @param params.set A lsit of parameters that were set to be included in the stat output file
#'
#'
#'
#'  @return
#'
#' Returns either a concise list of a subset of results or a full list with all possible results. Both output
#'
#' @author Kimberly J Gilbert
#'
#'
#'
#'
#'
#' @examples
#' 
#' 
#' @export read.output

read.output <- function(filename, params.set){
	if("test" %in% params.set) print("do things with that output statistic")
	print("foo")
}
