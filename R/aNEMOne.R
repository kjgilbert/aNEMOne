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
#' @param fst.dat  A data frame with the first column indicating population of origin and
#'          the following columns representing genotypes at loci; see the 
#'
#' @param qst.dat  the input table of the breeding data
#'  \itemize{
#'      \item    For \code{breeding.design = "half.sib.sire"}:  qst.dat should have four
#'           columns in this order: population, sire, dam, and the trait value of
#'           the individual. Each population, sire, and dam should have unique
#'           names or numbers.
#'      \item    For \code{breeding.design = "half.sib.dam"}:  qst.dat should have three
#'           columns in this order: population, dam, and the trait value of
#'           the individual. Each population and dam should have unique
#'           names or numbers.
#'  }
#'
#' @param numpops  number of populations in the sample
#'
#' @param nsim     number of simulation replicates to perform to create the null
#'          distributions and bootstraps
#' 
#'  @param AFLP whether or not to use AFLP data
#'
#'  @param breeding.design	the breeding design used when collecting the trait data
#'  There are two options for breeding design:
#'  \enumerate{
#'  \item 	"half.sib.sire" is a half sib design with dam nested within sire nested within
#'  		population which works for either balanced or unbalanced sampling designs
#'  \item 	"half.sib.dam" is a half sib design with dam nested within population
#'  		which works for either balanced or unbalanced sampling designs
#'  }
#'
#'  @param dam.offspring.relatedness relatedness between offspring in the dam model, default is 1/4, i.e. half-sib
#'
#'  @param output whether to output full or concise results, see details below
#'
#'  @return
#'
#' Returns either a concise list of a subset of results or a full list with all possible results. Both output
#'  options write the vector of Qst-Fst values to a text file.
#'
#'  Concise list returns (default)
#'  \itemize{
#'  \item 	the calculated difference between Qst and Fst with 95\% critical values,
#'  \item 	one- and two- tailed p-values for this difference,
#'  \item 	the Fst value estimated from the data with 95\% confidence intervals,
#'  \item 	the Qst value estimated from the data with 95\% confidence intervals, and
#'  \item 	the additive genetic variance for the trait with 95\% confidence intervals
#'  }
#'  Full list returns 
#'  \itemize{
#'  		\item  the calculated difference between Qst and Fst with 95\% critical values,
#'  		\item  one- and two- tailed p-values for this difference,
#'  		\item  a list of all Qst-Fst values for plotting the distribution of Qst-Fst,
#'  		\item  the Fst value estimated from the data with 95\% confidence intervals,
#'  		\item  the resampled Fst as calculated from bootstrapping across simulations, with standard deviation and 95\% confidence intervals,
#'  		\item  a list of all resampled Fst values for plotting the distribution of Fst,
#'  		\item  the Qst value estimated from the data with 95\% confidence intervals,
#'  		\item  the resampled neutral Qst as calculated from bootstrapping across simulations, with standard deviation and 95\% critical values,
#'  		\item  a list of all resampled Qst values for plotting the distribution of the neutral Qst,
#'  		\item  the ANOVA table for the calculated means squared, n coefficients, and degrees of freedom,
#'  		\item  the additive genetic variance for the trait with 95\% confidence intervals, and
#'  		\item  the coefficient of additive genetic variance for the trait with 95\% confidence intervals
#'  }
#'  Be mindful of the fact that with a small number of loci, bootstrapped confidence intervals of Fst can be less accurate.
#'
#' @author Kimberly J Gilbert
#'
#' @import hierfstat
#'
#' @references Gilbert KJ 
#'
#'
#'
#' @examples
#'
#' ## using balanced half-sib sire trait data and biallelic marker data 
#' data(hssire) # trait data
#' data(biallelic) # marker data
#' QstFstComp(biallelic, hssire, numpops=15, nsim=100, breeding.design="half.sib.sire", output="full")
#' 
#' data(hsdam)
#' data(aflp)
#' QstFstComp(aflp, hsdam, numpops=15, nsim=100, AFLP=TRUE, breeding.design="half.sib.dam", output="concise")
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
	LCE.order=NULL,
	mating.system=NULL,
	mean.fec=NULL,
	
	){
		row1 <- paste(c("run_mode", run.mode), collapse=" ")
		row2 <- paste(c("random_seed", random.seed), collapse=" ")
		row3 <- paste(c("logfile", log.file), collapse=" ")
		row4 <- paste(c("root_dir", root.dir), collapse=" ")
		row5 <- paste(c("filename", filename), collapse=" ")
		row6 <- paste(c("replicates", reps), collapse=" ")
		row7 <- paste(c("generations", gens), collapse=" ")
		row8 <- paste(c("patch_number", num.patches), collapse=" ")
		row9 <- NULL
		for(i in 1:length(LCE.order)){
			temp <- paste(paste(c(LCE.order[i], i), collapse=" "), sep="\n")
			row9 <- rbind(row9, temp)
		}
		
		row9 <- cat(c(), collapse=" ")
		row10 <- paste(c(), collapse=" ")
		row11 <- paste(c(), collapse=" ")
		row12 <- paste(c(), collapse=" ")
		row13 <- paste(c(), collapse=" ")
		row14 <- paste(c(), collapse=" ")
		row15 <- paste(c(), collapse=" ")
		row16 <- paste(c(), collapse=" ")
		row17 <- paste(c(), collapse=" ")
		row18 <- paste(c(), collapse=" ")
		row19 <- paste(c(), collapse=" ")
		row20 <- paste(c(), collapse=" ")
		
	
}




#'
#' @title Analyze a Nemo stat output file
#'
#' @description Compares the QST of a single phenotypic trait to the mean FST of series of marker loci. 
#'
#' @param fst.dat  A data frame with the first column indicating population of origin and
#'          the following columns representing genotypes at loci; see the 
#'
#' @param qst.dat  the input table of the breeding data
#'  \itemize{
#'      \item    For \code{breeding.design = "half.sib.sire"}:  qst.dat should have four
#'           columns in this order: population, sire, dam, and the trait value of
#'           the individual. Each population, sire, and dam should have unique
#'           names or numbers.
#'      \item    For \code{breeding.design = "half.sib.dam"}:  qst.dat should have three
#'           columns in this order: population, dam, and the trait value of
#'           the individual. Each population and dam should have unique
#'           names or numbers.
#'  }
#'
#' @param numpops  number of populations in the sample
#'
#' @param nsim     number of simulation replicates to perform to create the null
#'          distributions and bootstraps
#' 
#'  @param AFLP whether or not to use AFLP data
#'
#'  @param breeding.design	the breeding design used when collecting the trait data
#'  There are two options for breeding design:
#'  \enumerate{
#'  \item 	"half.sib.sire" is a half sib design with dam nested within sire nested within
#'  		population which works for either balanced or unbalanced sampling designs
#'  \item 	"half.sib.dam" is a half sib design with dam nested within population
#'  		which works for either balanced or unbalanced sampling designs
#'  }
#'
#'  @param dam.offspring.relatedness relatedness between offspring in the dam model, default is 1/4, i.e. half-sib
#'
#'  @param output whether to output full or concise results, see details below
#'
#'  @return
#'
#' Returns either a concise list of a subset of results or a full list with all possible results. Both output
#'  options write the vector of Qst-Fst values to a text file.
#'
#'  Concise list returns (default)
#'  \itemize{
#'  \item 	the calculated difference between Qst and Fst with 95\% critical values,
#'  \item 	one- and two- tailed p-values for this difference,
#'  \item 	the Fst value estimated from the data with 95\% confidence intervals,
#'  \item 	the Qst value estimated from the data with 95\% confidence intervals, and
#'  \item 	the additive genetic variance for the trait with 95\% confidence intervals
#'  }
#'  Full list returns 
#'  \itemize{
#'  		\item  the calculated difference between Qst and Fst with 95\% critical values,
#'  		\item  one- and two- tailed p-values for this difference,
#'  		\item  a list of all Qst-Fst values for plotting the distribution of Qst-Fst,
#'  		\item  the Fst value estimated from the data with 95\% confidence intervals,
#'  		\item  the resampled Fst as calculated from bootstrapping across simulations, with standard deviation and 95\% confidence intervals,
#'  		\item  a list of all resampled Fst values for plotting the distribution of Fst,
#'  		\item  the Qst value estimated from the data with 95\% confidence intervals,
#'  		\item  the resampled neutral Qst as calculated from bootstrapping across simulations, with standard deviation and 95\% critical values,
#'  		\item  a list of all resampled Qst values for plotting the distribution of the neutral Qst,
#'  		\item  the ANOVA table for the calculated means squared, n coefficients, and degrees of freedom,
#'  		\item  the additive genetic variance for the trait with 95\% confidence intervals, and
#'  		\item  the coefficient of additive genetic variance for the trait with 95\% confidence intervals
#'  }
#'  Be mindful of the fact that with a small number of loci, bootstrapped confidence intervals of Fst can be less accurate.
#'
#' @author Kimberly J Gilbert
#'
#' @import hierfstat
#'
#' @references Gilbert KJ 
#'
#'
#'
#' @examples
#'
#' ## using balanced half-sib sire trait data and biallelic marker data 
#' data(hssire) # trait data
#' data(biallelic) # marker data
#' QstFstComp(biallelic, hssire, numpops=15, nsim=100, breeding.design="half.sib.sire", output="full")
#' 
#' data(hsdam)
#' data(aflp)
#' QstFstComp(aflp, hsdam, numpops=15, nsim=100, AFLP=TRUE, breeding.design="half.sib.dam", output="concise")
#' 
#' 
#' @export read.output

read.output <- function(y){
	
	
}
