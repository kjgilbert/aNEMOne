# CREATE THE LANDSCAPE OF PHENOTYPIC OPTIMA FOR SELECTION


#' @title Create the landscape of optimal phenotypic values
#'
#' @description Create the landscape of optimal phenotypic values over which the patches will exist. Size should match the number of cells in the x and y directions on the landscape. A figure is created to visualize the landscape as well as the matrix of values output to a file for use in the make.input function.
#'
#'  @param horizontal.patches Number of patches in the horizontal, x direction.
#'
#'  @param vertical.patches Number of patches in the vertical, y direction.
#' 
#'  @param range The distance at which covariance between patches goes to zero, i.e. how similar is each grain to the next. A larger value makes the landscape smoother.
#'
#'  @param sill The scale of the variance. See ?vgm, this is the sill of the variogram model component. A larger value makes the landscape more patchy. Default value is 1.
#'
#'  @param magnitude The first beta parameter for simple kriging in the gstat function. Default is zero.
#'
#'  @param directionality The second beta parameter for simple kriging in the gstat function. Default is 1. Larger values smooth the change in the direction of change.
#'
#'  @param nearest.obs Used for local kriging, see ?gstat: the number of nearest observations that should be used for a kriging prediction or simulation. Default is 20. 
#'
#'  @param scale.to Value to scale down the magnitude of change in optima over the total landscape. Default is 1, which scales the total change to approximately -1 to 1 across the full landscape, whereas, e.g. a value of 5 allows a higher degree of change across the landscape, approximately -5 to 5.
#'
#'  @param cylinder If the landscape is meant to be simulated as a cylinder, where the horizontal edges match up for individuals to disperse across, this should be set to true. This will halve the vertical length and mirror the landscape so that to horizontal edges match where they meet. Currently this capability is not yet coded in.
#'
#'  @return
#'
#' Returns the mean value of the first column of patches at the leftmost end of the landscape. Also prints to a file the matrix of values for the landscape, to be uses as "selection_local_optima" in Nemo, and additionally visualizes the landscape.
#'
#' @author Kimberly J Gilbert
#'
#' @import fields gstat
#'
#' @examples
#'
#' make.landscape(horizontal.patches=50, vertical.patches=25, range=5, sill=1)
#' 
#' 
#' @export make.landscape


make.landscape <- function(horizontal.patches, vertical.patches, range, sill=1, magnitude=0, directionality=1, nearest.obs=20, scale.to=1, cylinder=FALSE){
	
	# make the color palette for the visualization
	blues <- function(n){
		hsv(h=0.65,	#blue
		s=c(seq(1,0,length.out=(n))),
		v=c(seq(0,1,length.out=(n/2)), rep(1,(n/2)))
		)
	}

	if(cylinder == TRUE){
		print("This capability not yet enabled.")
		return(NULL)
	}

	xy <- expand.grid(1: horizontal.patches, 1: vertical.patches)
	names(xy) <- c("x", "y")

	# first number magnitude, second: directionality/step-ness, third: leave at 0 (angle of direction)
	trend <- matrix(c(magnitude,directionality,0))
		# Katie tested these at 1, 0.01, 0.1, 0 and combos of those

	#sill <- sill	# leave at 1 (higher also makes it more patchy) this is the scale of variance
	#range <- range	# larger makes less patchy - distance at which covariance goes to zero - how similar is each grain 		to the next
		# Katie tested ranges of 10, 30, 60, 120, 200
	#nearest.obs <- 20
	
	reps <- 1
	g.dummy <- gstat(formula=z~1+x+y, locations=~x+y, dummy=TRUE, beta=trend, model=vgm(psill=sill, model="Sph", 	range=range), nmax=nearest.obs)
	yy <- predict(g.dummy, newdata=xy, nsim=reps)
	beta.char <- paste(trend, collapse="_")
	env <- as.numeric(yy[,2+1])-mean(yy[,2+1])
	mat.env <- matrix(env, ncol= vertical.patches, nrow= horizontal.patches, byrow=FALSE)

	first.column.mean <- mean(mat.env[1,]) # yes it is still actually the mean of the column because of how they got put into rows

	mat.env.scaled <- mat.env/(abs(first.column.mean)/scale.to)
	scaled.first.column.mean <- mean(mat.env.scaled[1,])
	# visualize 
	# UNSCALED PLOT # image.plot(1: horizontal.patches, 1: vertical.patches, mat.env, col=blues(100))
	image.plot(1: horizontal.patches, 1: vertical.patches, mat.env.scaled, col=blues(100), ylab="Vertical Patches", xlab="Horizontal Patches")
	
	# transpose it so that what is actually a column (currently a row) becomes a column fo the nemo file
	t.mat.env <- t(mat.env.scaled)
	
	if(scaled.first.column.mean != mean(t.mat.env[,1])){
		#then the transpose didn't work properly because what was a row has not now become a column
		print("Transpose did not work properly.")
		stop("Rows did not properly become columns.")
	}
	
	# make it a data frame so can paste to file for other analyses
	frame.mat.env <- data.frame(t.mat.env)
	
	write.table(frame.mat.env, file="Landscape_Matrix.txt", sep=",", col.names=FALSE, row.names=FALSE)

	array.mat.env <- NULL
	for(i in 1:(dim(t.mat.env)[2])){
		temp.array.mat.env <- unlist(t.mat.env[,i])
		array.mat.env <- c(array.mat.env, temp.array.mat.env)
	}
	
	with.brackets <- paste(array.mat.env, collapse="}{")
	final.landscape.array <- paste(c("{{", with.brackets, "}}"), collapse="")	
	write.table(final.landscape.array, file=paste(c(getwd(), "/", "Landscape.txt"), collapse=""), col.names=FALSE, row.names=FALSE, quote=FALSE)
		
	
	return(scaled.first.column.mean)
}







#' @title Create the landscape of optimal phenotypic values that change at an equal and constant pace over space
#'
#' @description Create the landscape of optimal phenotypic values over which the patches will exist. Size should match the number of cells in the x and y directions on the landscape. A figure is created to visualize the landscape as well as the matrix of values output to a file for use in the make.input function.
#'
#'  @param horizontal.patches Number of patches in the horizontal, x direction.
#'
#'  @param vertical.patches Number of patches in the vertical, y direction.
#'
#' @param step.width The number of patches that horizontally make up the width of "steps" of equal optimal phenotype to create the landscape. At most can equal the number of vertical columns in the landscape to create a uniform landscape, or at minimum can equal 1 for the landscape to change every single column of patches. This number must equally divide into the total number of horizontal patches.
#'
#'  @param scale.to Value to scale down the magnitude of change in optima over the total landscape. Default is 1, which scales the total change to approximately -1 to 1 across the full landscape, whereas, e.g. a value of 5 allows a higher degree of change across the landscape, approximately -5 to 5.
#'
#'
#'
#'
#'  @return
#'
#' Returns the value of the first column of patches at the leftmost end of the landscape. Also prints to a file the matrix of values for the landscape, to be uses as "selection_local_optima" in Nemo, and additionally visualizes the landscape.
#'
#' @author Kimberly J Gilbert
#'
#' @import fields
#'
#' @examples
#'
#' step.landscape(horizontal.patches=100, vertical.patches=25, step.width=10)
#' 
#' 
#' @export step.landscape


step.landscape <- function(horizontal.patches, vertical.patches, step.width, scale.to=1){
	# make the color palette for the visualization
	blues <- function(n){
		hsv(h=0.65,	#blue
		s=c(seq(1,0,length.out=(n))),
		v=c(seq(0,1,length.out=(n/2)), rep(1,(n/2)))
		)
	}

	first.col <- -scale.to
	last.col <- scale.to
	
	num.steps <- horizontal.patches/step.width
	
	step.values <- seq(first.col, last.col, length.out= num.steps)
	
	total.num.patches <- horizontal.patches*vertical.patches
	
	patches.per.step <- total.num.patches/num.steps
	
	env <- NULL
	for(i in 1:num.steps){
		one.step <- rep(step.values[i], patches.per.step)
		env <- c(env, one.step)
	}
	
	mat.env <- matrix(env, ncol= vertical.patches, byrow=TRUE)
	image.plot(1: horizontal.patches, 1: vertical.patches, mat.env, col=blues(100), ylab="Vertical Patches", xlab="Horizontal Patches")

	# transpose it to have rows and columns right way round for looking at text
	t.mat.env <- t(mat.env)
	
	if(first.col != mean(t.mat.env[,1])){
		#then the transpose didn't work properly because what was a row has not now become a column
		print("Transpose did not work properly.")
		stop("Rows did not properly become columns.")
	}
	
	# make it a data frame so can paste to file for other analyses
	frame.mat.env <- data.frame(t.mat.env)

	write.table(frame.mat.env, file="Landscape_Matrix.txt", sep=",", col.names=FALSE, row.names=FALSE)


	array.env <- paste(env, collapse="}{")
	final.landscape.array <- paste(c("{{", array.env, "}}"), collapse="")	
	write.table(final.landscape.array, file=paste(c(getwd(), "/", "Landscape.txt"), collapse=""), col.names=FALSE, row.names=FALSE, quote=FALSE)
		
	
	return(first.col)
}
