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
#' @title Compare Qst to Fst
#'
#' @param fst.dat  A data frame with the first column indicating population of origin and
#'          the following columns representing genotypes at loci; see the 
#'          README \url{https://github.com/kjgilbert/QstFstComp/blob/master/README.md} for further description.
#'          If using AFLPs, this is a data frame of q-hat values, with pops in columns, loci in rows
#'          and the corresponding q-hat variances in the following columns, and \code{AFLP=TRUE} must be designated.
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
#' @references Whitlock MC and F Guillaume (2009) \href{http://www.genetics.org/content/183/3/1055}{Testing for spatially divergent selection: Comparing \emph{Qst} to \emph{Fst}.} \emph{Genetics}, 183:1055-1063.
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



make.input <- function(){
	
}