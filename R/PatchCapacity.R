# make the array of patch capacities



#' @title Create array of patch capacities
#'
#' @description Make the linear array pf patch capacities to be used in the Nemo .init file.
#'
#'  @param K1 The carrying capacity for patches in the array. If there will be two different carrying capacities, this is for the first set of patches in the array.
#'
#'  @param num.K1 The number of patches having carrying capacity, K1.
#' 
#'  @param K2 The carrying capacity for the second set of patches in the array. Only needed if there is more than one K for all patches.
#'
#'  @param num.K2 The number of patches having carrying capacity, K2. Only needed if there is more than one K for all patches.
#'
#'  @return
#'
#' Returns the array of patch capacities that can then be provided to the 'make.input' function for the final input file.
#'
#' @author Kimberly J Gilbert
#'
#' @references Gilbert KJ 
#'
#' @examples
#'
#' patch.cap(K1=10, num.K1=12, K2=0, num.K2=4)
#' 
#' 
#' @export patch.cap


patch.cap <- function(K1, num.K1, K2=NULL, num.K2=NULL){
	first.values <- paste(rep(K1, num.K1), collapse=",")
	if(!is.null(K2)){
		second.values <- paste(rep(K2, num.K2), collapse=",")
		capacities <- paste(c(first.values, second.values), collapse=",")
		cap.array <- paste("{{", capacities, "}}")
	}else{
		cap.array <- paste("{{", first.values, "}}")
	}
	return(cap.array)
}