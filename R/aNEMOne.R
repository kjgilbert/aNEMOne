# Functions to create Nemo input and analyse Nemo output 
#	by Kimberly J. Gilbert (kgilbert@zoology.ubc.ca)

# made to work with code edits by Kim for breeding window and connectivity matrix with dspersal kernel, see: https://github.com/kjgilbert/NemoDispersalKernel


#'
#' Makes an input file of type .ini for \href{http://nemo2.sourceforge.net/index.html}{Nemo}, an individual-based, forward time simulation program created and maintained by Fred Guillaume.
#'
#' @title Create Nemo input file
#'
#'
#'  @param run.mode The mode to run the simulation in; default is "overwrite", for other options see the Nemo manual. 
#'
#'  @param random.seed Set the random seed for the run, default is 12345.
#'
#'  @param log.file Name of the logfile to be output.
#'
#'  @param root.dir Root directory to put run outputs.
#'
#'  @param filename Base name of the simulation file.
#'
#'  @param reps Number of replicates to perform.
#'
#'  @param gens Number of generations to run the simulation.
#'
#'  @param num.patches Number of patches on the landscape.
#'
#'  @param patch.capacity The carrying capacity of each patch, can be specified by an integer value or an array of length num.patches.
#'
#'  @param patch.capacity.fem If sex-specific carrying capacity, this is the female carrying capacity per patch.
#'
#'  @param patch.capacity.mal If sex-specific carrying capacity, this is the male carrying capacity per patch.
#'
#'  @param cap.temp If carrying capacity changes at a given generation, this is the generation at which that occurs.
#'
#'  @param LCE.order The order of life cycle events in the simulation. Should be given as multiple character strings, in order.
#'
#'  @param mating.system See the nemo manual for details, random mating is specified by a '1'.
#'
#'  @param mating.proportion See the nemo manual for details, sets the proportion of non-random mating.
#'
#'  @param mean.fec Mean fecundity per mother, may also be set to be patch-specific with an array. By default follows a Poisson distribution, see the Nemo manual for details.
#'
#'  @param self.if.alone Boolean, if true, an individual will self if it finds no mate.
#'
#'  @param always.breed.window Boolean, if present, the breeding window is always used. Otherwise, the breeding window is only called if no mate can be found in the focal patch.
#'
#'  @param breeding.connectivity.matrix The connectivity matrix of patches matched to the breeding kernel.
#'
#'  @param breeding.kernel The array of probabilities of searching for a mate within a given patch.
#'
#'  @param dispersal.connectivity.matrix The connectivity matrix of patches matched to the dispersal kernel.
#'
#'  @param dispersal.kernel The array of dispersal probabilities for forward migration.
#'
#'  @param seln.trait The trait(s) specified to be under selection. Either "quant" or "delet" for quantitative or deleterious traits.
#'
#'  @param seln.model The model of selection to use. See Nemo manual section 4.7 for all details on selection.
#'
#'  @param seln.fitness.model Options are absolute, relative_global, and relative_local. Default is absolute. See Nemo manual section 4.7 for all details on selection.
#'
#'  @param seln.var Variance in selection.
#'
#'  @param seln.trait.dim See Nemo manual section 4.7 for all details on selection. Default is 1.
#'
#'  @param seln.local.optima Optimal trait value per patch or for the entire landscape
#'
#'  @param quanti.init The initial QTL value to be set; see Nemo manual section 5.3 for all details on quantitative traits.
#'
#'  @param num.quanti.traits How many QTL traits to model, default is 1.
#'
#'  @param num.quanti.loci How many loci underly the quantitative trait.
#'
#'  @param quanti.mut.rate Mutation rate per QTL.
#'
#'  @param quanti.mut.var Variance in mutation rate per QTL.
#'
#'  @param quanti.recomb.rate Recombination rate among QTL, default is freely recombining, 0.5.
#'
#'  @param quanti.init.model See Nemo manual section 5.3 for all details on quantitative traits. Default is 1.
#'
#'  @param quanti.env.var Environmental variance, default is 1.
#'
#'  @param num.ntrl.loci The number of neutral loci to simulate. See Nemo manual section 5.2 for all details on neutral markers.
#'
#'  @param num.ntrl.alleles The number of alleles per neutral locus.
#'
#'  @param ntrl.mut.rate The mutation rate per neutral locus.
#'
#'  @param ntrl.recomb.rate The recombination rate among neutral loci.
#'
#'  @param ntrl.mut.model The mutation model for neutral loci: 1 = single step, 2 = K allele model. See Nemo manual section 5.2 for all details on neutral markers.
#'
#'  @param ntrl.init.model How to initiate the neutral allele frequencies: 0 = no initial variance, 1 = maximum initial variance. See Nemo manual section 5.2 for all details on neutral markers.
#'
#'  @param save.ntrl How often to save the neutral marke genotype files. Files are automatically output to subdirectory "ntrl_geno". If not present, no neutral genotype output is saved.
#'
#'  @param save.quanti How often to save the quantitative trait genotype files. Files are automatically output to subdirectory "quanti_geno". If not present, no quantitative genotype output is saved.
#'
#'  @param save.stats How often to save the values of the parameters defined by "stats". Files are automatically output to subdirectory "stats". *If this paramter is not present, no stats are output, even if "stats" is defined.
#'
#'  @param stats Population and simlation parameters to return, see Nemo manual section 7 "Output Statistics".
#'
#'  @return
#'
#'  Write the file specified into the working directory.
#'
#' @author Kimberly J Gilbert
#'
#' @references \href{http://nemo2.sourceforge.net/index.html}{Nemo} is created and maintained by Fred Guillaume. The manual and source files are available online.
#'
#' @examples
#'
#' # population
#' horiz.patches <- 100
#' vert.patches <- 20
#' num.patches <- (horiz.patches*vert.patches) + 2	# plus two for the ghost and absorbing patches
#' cell.size <- 50
#' sigma <- 25
#' 
#' land.x <- cell.size*horiz.patches
#' land.y <- cell.size* vert.patches
#' 
#' cap <- patch.cap(6,1000, 0,2)
#'
#' # simulation components
#' life.cycles <- c("breed", "disperse", "selection", "aging")
#'
#' # mating
#' breed.kernel <- make.kernel.and.matrix(cell.size= cell.size, horizontal.land= land.x,
#'  vertical.land= land.y, dist.mean=0, dist.sd= sigma, breed.window=TRUE, two.kernels=FALSE,
#'  second.dist.mean=0, second.dist.sd=NULL)
#' breed.file <- "Breeding_ConnMatrix.txt"
#' breeding.connectivity.matrix <- readChar(breed.file, file.info(breed.file)$size)
#'
#' # dispersal
#' disp.kernel <- make.kernel.and.matrix(cell.size= cell.size, horizontal.land= land.x,
#'  vertical.land= land.y, dist.mean=0, dist.sd= sigma, breed.window=FALSE, two.kernels=FALSE,
#'  second.dist.mean=0, second.dist.sd=NULL)
#' disp.file <- "Dispersal_ConnMatrix.txt"
#' dispersal.connectivity.matrix <- readChar(disp.file, file.info(disp.file)$size)
#'
#' # landscape
#' landscape.init.optima <- make.landscape(horizontal.patches= horiz.patches,
#'  vertical.patches= vert.patches, range=5)
#' landscape.file <- "Landscape.txt"
#' landscape <- readChar(landscape.file, file.info(landscape.file)$size)
#' 
#'
#' make.input(
#'	root.dir="test", filename="test",
#'	reps=100, gens=1000, num.patches= num.patches, patch.capacity=cap,
#'	LCE.order= life.cycles,	
#'	mating.system=1, mean.fec=7,
#'	breeding.connectivity.matrix= breeding.connectivity.matrix,
#'	breeding.kernel= breed.kernel,
#'	dispersal.connectivity.matrix= dispersal.connectivity.matrix,
#'	dispersal.kernel= disp.kernel,
#'	seln.trait="quanti", seln.model="gaussian", seln.var=7.5, seln.trait.dim=1,	seln.local.optima= landscape, 
#'	quanti.init= landscape.init.optima, num.quanti.loci=100, quanti.mut.rate=0.001, quanti.mut.var=0.01,
#'	num.ntrl.loci=100, num.ntrl.alleles=2, ntrl.mut.rate=0.001, ntrl.mut.model=1, ntrl.init.model=1, 
#'	save.ntrl=50, save.quanti=10, save.stats=10,
#'	stats=c("demography", "fecundity", "migrants")
#' )
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
	mating.proportion=NULL,
	mean.fec=NULL,
	self.if.alone=FALSE,
	always.breed.window=FALSE,
	breeding.connectivity.matrix=NULL,
	breeding.kernel=NULL,
	dispersal.connectivity.matrix=NULL,
	dispersal.kernel=NULL,
	seln.trait=NULL,
	seln.model=NULL,
	seln.fitness.model="absolute",
	seln.var=NULL,
	seln.trait.dim=1,
	seln.local.optima=NULL,
	quanti.init=NULL,
	num.quanti.traits=1,
	num.quanti.loci=NULL,
	quanti.mut.rate=NULL,
	quanti.mut.var=NULL,
	quanti.recomb.rate=0.5,
	quanti.init.model=1,
	quanti.env.var=1,
	num.ntrl.loci=NULL,
	num.ntrl.alleles=NULL,
	ntrl.mut.rate=NULL,
	ntrl.recomb.rate=0.5,
	ntrl.mut.model=NULL,	# 1 = single step, 2 = K allele model
	ntrl.init.model=NULL,	# 0 = no initial variance, 1 = max. initial variance
	save.ntrl=NULL,
	save.quanti=NULL,
	save.stats=NULL,
	stats=NULL
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
		if(!is.null(quanti.init)) row10 <- paste("quanti_init 0")
		for(i in 1:length(LCE.order)){
			temp <- paste(paste(c(LCE.order[i], i), collapse=" "), sep="\n")
			row10 <- paste(c(paste(row10), paste(temp), sep="\n"))
		}
		
		row11 <- paste(c("mating_system", mating.system), collapse=" ")

		if(is.null(mating.proportion)){
			row11 <- paste(c("mating_system", mating.system), collapse=" ")
		}else{
			row11 <- paste(c(paste(c("mating_system", mating.system), collapse=" "), paste(c("mating_proportion", mating.proportion), collapse=" ")), collapse="\n")
		}
		row12 <- paste(c("mean_fecundity", mean.fec), collapse=" ")
		if(self.if.alone==TRUE){
			row12 <- paste(c(row12, paste("self_if_alone")), collapse="\n")
		}
		if(always.breed.window ==TRUE){
			row12 <- paste(c(row12, paste("always_breed_window")), collapse="\n")
		}
		row13 <- paste(c("breeding_connectivity_matrix {", breeding.connectivity.matrix, "}"), collapse=" ")
		row14 <- paste(c("breeding_kernel {", breeding.kernel, "}"), collapse=" ")
		row15 <- paste(c("dispersal_connectivity_matrix {", dispersal.connectivity.matrix, " }"), collapse=" ")
		row16 <- paste(c("dispersal_kernel {", dispersal.kernel, "}"), collapse=" ")
	
		row17 <- paste("\n## SELECTION TRAITS")
		row18 <- paste(c("selection_trait", seln.trait), collapse=" ")
		row19 <- paste(c("selection_model", seln.model), collapse=" ")
		row20 <- paste(c("selection_fitness_model", seln.fitness.model), collapse=" ")
		row21 <- paste(c("selection_variance", seln.var), collapse=" ")
		row22 <- paste(c("selection_trait_dimension", seln.trait.dim), collapse=" ")
		row23 <- paste(c("selection_local_optima", seln.local.optima), collapse=" ")
	
		row24 <- paste("\n## QUANTI TRAITS")
		row25 <- paste(c("quanti_init_value {{", quanti.init, "}}"), collapse=" ")
		row26 <- paste(c("quanti_traits", num.quanti.traits), collapse=" ")
		row27 <- paste(c("quanti_loci", num.quanti.loci), collapse=" ")
		row28 <- paste(c("quanti_mutation_rate", quanti.mut.rate), collapse=" ")
		row29 <- paste(c("quanti_mutation_variance", quanti.mut.var), collapse=" ")
		row30 <- paste(c("quanti_recombination_rate", quanti.recomb.rate), collapse=" ")
		row31 <- paste(c("quanti_init_model", quanti.init.model), collapse=" ")
		row32 <- paste(c("quanti_environmental_variance", quanti.env.var), collapse=" ")
	
		row33 <- paste("\n## NEUTRAL TRAITS")
		row34 <- paste(c("ntrl_loci", num.ntrl.loci), collapse=" ")
		row35 <- paste(c("ntrl_all", num.ntrl.alleles), collapse=" ")
		row36 <- paste(c("ntrl_mutation_rate", ntrl.mut.rate), collapse=" ")
		row37 <- paste(c("ntrl_recombination_rate", ntrl.recomb.rate), collapse=" ")
		row38 <- paste(c("ntrl_mutation_model", ntrl.mut.model), collapse=" ")
		row39 <- paste(c("ntrl_init_model", ntrl.init.model), collapse=" ")

		row40 <- paste("\n## OUTPUT")
		row41 <- NULL
		if(!is.null(save.stats)){
			row41 <- paste(c(
			row41,
			paste(c("stat_dir stats", 
				paste(c("stat_log_time", save.stats), collapse=" "), 
				paste(c("stat", stats), collapse=" "),
				paste("stat_output_compact"),
				paste("stat_output_CSV"),
				paste("stat_output_precision 4")), collapse="\n"),
				collapse="\n"
			))
		}
		if(!is.null(save.ntrl)){
			row41 <- paste(c(
			row41,
			paste(c("ntrl_save_genotype FSTAT", "ntrl_output_dir ntrl_geno", 
				paste(c("ntrl_output_logtime", save.ntrl), collapse=" ")), collapse="\n"), 
				collapse="\n"
			))
		}
		if(!is.null(save.quanti)){
			row41 <- paste(c(
			row41,
			paste(c("quanti_output genotypes", "quanti_dir quanti_geno", 
				paste(c("quanti_logtime", save.ntrl), collapse=" ")), collapse="\n"), 
				collapse="\n"
			))
		}			
	
	init.file <- paste(c(row1, row2, row3, row4, row5, "\n", row6, row7, row8, row9, row10, row11, row12, row13, row14, row15, row16, row17, row18, row19, row20, row21, row22, row23, row24, row25, row26, row27, row28, row29, row30, row31, row32, row33, row34, row35, row36, row37, row38, row39, row40, row41), collapse="\n")
	
		
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
#' 
#' 
#' @export read.output

read.output <- function(filename, params.set){
	if("test" %in% params.set) print("do things with that output statistic")
	print("foo")
}
