# library(devtools)
# install_github("kjgilbert/aNEMOne")

library(aNEMOne)

setwd("~/Desktop")


horiz.patches <- 500
vert.patches <- 10
cell.size <- 10		# units, used to scale sigma
sigma.breed <- 5	# units, scaled relative to cell size
sigma.disp <- 30	# units, scaled relative to cell size



variance.seln <- 5	
fecundity <- 5
carrying.cap <- 10

# quanti loci
quanti.mut <- "0.0001"
quant.mut.var <- 0.05
q.env.var <- 1

# deleterious loci
delet.mut <- "0.00005"
gamma.mean <- 0.01
gamma.shape <- 0.3
del.dom.mean <- 0.3

# steepness of gradient, b
# set to 0 for no gradient
b <- 0.1

name.of.file <- "Example1"

# number of replicates
nreps <- 1
# number of generations
num.gens <- 5000
# at what generation to open the rest of the habitat for expansion
expand.at <- 1000

# how often to save outputs
save.delet.at <- 1000
save.quanti.at <- 1000
save.stats.at <- 50
rand.seed <- "12345"
root.directory <- "Example1_Directory"
stats.to.save <- c("demography", "extrate", "fecundity", "migrants", "delet", "quanti", "quanti.mean.patch", "fitness", "adlt.fitness.patch")

# simulation components of life cycle
life.cycles <- c("breed", "disperse", "viability_selection", "aging", "save_stats", "save_files")




#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#

num.patches <- (horiz.patches*vert.patches) + 2	# plus two for the ghost and absorbing patches
land.x <- cell.size*horiz.patches  # square units
land.y <- cell.size* vert.patches # square units

# make the landscape
if(b==0){
	landscape.init.optima <- 0
	landscape <- "{{0}}"
}else{
	scale.value <- (b*horiz.patches)/2
	landscape.init.optima <- step.landscape(horizontal.patches= horiz.patches, vertical.patches= vert.patches, step.width=1, scale.to=scale.value)
	landscape.file <- "Landscape.txt"
	landscape <- readChar(landscape.file, file.info(landscape.file)$size)
	landscape.init.optima <- -(scale.value-(vert.patches*b)) 	# sets pops to initialize at the mid-value of the landscape over which the pop is initialized
}

# make the dispersal kernel
disp.kernel <- make.kernel.and.matrix(cell.size= cell.size, horizontal.land= land.x, vertical.land= land.y, dist.mean=0, dist.sd= sigma.disp, breed.window=FALSE, two.kernels=FALSE, second.dist.mean=0, second.dist.sd=NULL)
disp.file <- "Dispersal_ConnMatrix.txt"
dispersal.connectivity.matrix <- readChar(disp.file, file.info(disp.file)$size)

# make mating window
breed.kernel <- make.kernel.and.matrix(cell.size= cell.size, horizontal.land= land.x, vertical.land= land.y, dist.mean=0, dist.sd= sigma.breed, breed.window=TRUE, two.kernels=FALSE, second.dist.mean=0, second.dist.sd=NULL)
breed.file <- "Breeding_ConnMatrix.txt"
breeding.connectivity.matrix <- readChar(breed.file, file.info(breed.file)$size)

# population -- because hermaphrodites, have to multiply desired carrying cap by 2!
cap1 <- patch.cap((2*carrying.cap),(vert.patches^2), 0,(num.patches-(vert.patches^2)), temp.gen=0)
cap2 <- patch.cap((2*carrying.cap),(num.patches-2), 0,2, temp.gen=expand.at)


make.delet.input(
	run.mode="overwrite", random.seed=rand.seed, log.file="logfile.log", root.dir=root.directory,
	filename=name.of.file,
	reps=nreps,	gens=num.gens, num.patches= num.patches,
	patch.capacity=paste(c("(", cap1, ",", cap2, ")"), collapse="\n"),
	patch.capacity.fem=NULL, patch.capacity.mal=NULL, cap.temp=FALSE, LCE.order= life.cycles,
	mating.system=4, mating.proportion=0,
	mean.fec=fecundity,	self.if.alone=FALSE, always.breed.window=TRUE, large.kernels=TRUE,
	breeding.connectivity.matrix= breeding.connectivity.matrix, breeding.kernel= breed.kernel,
	dispersal.connectivity.matrix= dispersal.connectivity.matrix, dispersal.kernel= disp.kernel,
	seln.trait="(delet, quant)",
	seln.model="(direct, gaussian)",
	seln.fitness.model="absolute",
	seln.var=variance.seln, seln.trait.dim=1, seln.local.optima= landscape, quanti.init= landscape.init.optima, num.quanti.traits=1,
	num.quanti.loci=100,
	quanti.mut.rate= quanti.mut, quanti.mut.var= quant.mut.var, quanti.recomb.rate= "{{100,100,100,100,100,100,100,100,100,100}}",
	quanti.init.model=1, quanti.env.var= q.env.var,
	num.delet.loci=1000, delet.mut.rate= delet.mut,
	delet.mut.model=1, delet.recomb.rate= "{{100,100,100,100,100,100,100,100,100,100}}",
	delet.init.freq=0, delet.effects.dist= "gamma", delet.effects.mean= gamma.mean, delet.effects.dist.param1= gamma.shape,
	delet.effects.dist.param2=NULL, delet.dominance.mean= del.dom.mean,
	delet.fitness.model= 1, save.delet=save.delet.at, save.quanti=save.quanti.at, save.stats=save.stats.at, stats=stats.to.save
)
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#


# Make a second input with a different fecundity but all else the same
name.of.file <- "Example2_HigherFecund"
fecundity <- 10


make.delet.input(
	run.mode="overwrite", random.seed=rand.seed, log.file="logfile.log", root.dir=root.directory,
	filename=name.of.file,
	reps=nreps,	gens=num.gens, num.patches= num.patches,
	patch.capacity=paste(c("(", cap1, ",", cap2, ")"), collapse="\n"),
	patch.capacity.fem=NULL, patch.capacity.mal=NULL, cap.temp=FALSE, LCE.order= life.cycles,
	mating.system=4, mating.proportion=0,
	mean.fec=fecundity,	self.if.alone=FALSE, always.breed.window=TRUE, large.kernels=TRUE,
	breeding.connectivity.matrix= breeding.connectivity.matrix, breeding.kernel= breed.kernel,
	dispersal.connectivity.matrix= dispersal.connectivity.matrix, dispersal.kernel= disp.kernel,
	seln.trait="(delet, quant)",
	seln.model="(direct, gaussian)",
	seln.fitness.model="absolute",
	seln.var=variance.seln, seln.trait.dim=1, seln.local.optima= landscape, quanti.init= landscape.init.optima, num.quanti.traits=1,
	num.quanti.loci=100,
	quanti.mut.rate= quanti.mut, quanti.mut.var= quant.mut.var, quanti.recomb.rate= "{{100,100,100,100,100,100,100,100,100,100}}",
	quanti.init.model=1, quanti.env.var= q.env.var,
	num.delet.loci=1000, delet.mut.rate= delet.mut,
	delet.mut.model=1, delet.recomb.rate= "{{100,100,100,100,100,100,100,100,100,100}}",
	delet.init.freq=0, delet.effects.dist= "gamma", delet.effects.mean= gamma.mean, delet.effects.dist.param1= gamma.shape,
	delet.effects.dist.param2=NULL, delet.dominance.mean= del.dom.mean,
	delet.fitness.model= 1, save.delet=save.delet.at, save.quanti=save.quanti.at, save.stats=save.stats.at, stats=stats.to.save
)
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#


# make the PBS files to run these Nemo inputs on Westgrid's bugaboo cluster:

# **MUST HAVE THE EXTRA FILES "LastLine.txt" and "MiddleLine.txt" from Github here: https://github.com/kjgilbert/aNEMOne/tree/master/extra in the same directory folder!

ini.dir <- "~/Desktop"
multi.pbs(ini.file.directory= ini.dir, RAM="6gb", server="bugaboo")

