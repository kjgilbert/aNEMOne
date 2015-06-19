#'
#' Take an individual Nemo .ini file and make a Westgrid PBS file for either orcinus, bugaboo, or grex servers.
#'
#' @title Make Westgrid PBS files
#'
#'
#'  @param filename The full name of the .ini file for which the PBS file is being created.
#'
#'  @param input.directory The directory on the Westgrid server which holds the nemo .ini file. Defaults are set for each server.
#'
#'  @param nemo.directory The directory on the Westgrid server from which the nemo executable is called. Defaults are set for each server.
#'
#'  @param RAM The amount of RAM the run will be allowed to use, there is no default, and must include gb, mb, or kb in specification, e.g. "6gb".
#'
#'  @param walltime Default is set at max for specified server (except bugaboo, which is set to default at 500 hours), but other lengths can be specified here in units "100:00:00" hours:minutes:seconds.
#'
#'  @param westgrid.server Which westgrid server the PBS file is meant for. Options are "grex", "orcinus" or "bugaboo" currently.
#'
#'  @return
#'
#'  Creates PBS files from a .ini file for running Nemo on Westgrid's Orcinus, Bugaboo, or Grex servers. *IMPORTANT NOTE* Must have files "MiddleLine.txt", "GrexLine.txt", and "LastLine.txt" in the same directory for make.pbs to work properly, as it concatenates text from those files into the final PBS file. 
#'
#' @author Kimberly J Gilbert
#'
#' @references \href{http://nemo2.sourceforge.net/index.html}{Nemo} is created and maintained by Fred Guillaume. The manual and source files are available online.
#'
#' @export make.pbs

make.pbs <- function(filename, input.directory=NULL, nemo.directory=NULL, RAM, walltime=NULL, westgrid.server=NULL){
	# currently have Nemo installed on: orcinus, bugaboo, and grex (though no compile abilities on grex, just copy over from orcinus)
	#
	# orcinus, walltime limit = 240:00:00  (10 days)
	#		orcinus RAM can be pretty slow
	#		/home/kgilbert/ fits 250 GB spae, /global/scratch fits 500 GB space
	# bugaboo, walltime limit = 122 days, so as much as needed
	#		seems slightly better about RAM and starting runs than orcinus, but maybe b/c I'm not sharing it w/ Remi and Jeremy
	#		/home/kgilbert/ is almost full, but global/scratch has 1TB memory
	# grex, walltime limit = 168:00:00  (7 days)
	#		much better for more ram, can use up to 48GB per node!
	#		must do all runs in /global/scratch/kgilbert/, so no backup and only 800GB memory
	if(is.null(westgrid.server)){
		print("Must specify 'grex', 'bugaboo', or 'orcinus' as the westgrid server being used.")
	}
	if(is.null(nemo.directory)){
		if(westgrid.server=="grex"){
			nemo.directory <- "/global/scratch/kgilbert/Nemo/bin/nemo2.3.38"
			if(is.null(input.directory)) input.directory <- "/global/scratch/kgilbert/RangeExpansion/"
		}
		if(westgrid.server=="bugaboo"){
			nemo.directory <- "/home/kgilbert/Nemo/bin/nemo2.3.38"
			if(is.null(input.directory)) input.directory <- "/home/kgilbert/RangeExpansion/My_HetLand_Project/"
		}	
		if(westgrid.server=="orcinus"){
			nemo.directory <- "/home/kgilbert/Nemo2.3/bin/nemo2.3.38"
			if(is.null(input.directory)) input.directory <- "/home/kgilbert/NemoRuns/MyRangeExpansionProject/"
		}
	}
	if(westgrid.server=="grex"){
		if(is.null(walltime)) walltime <- "168:00:00"
		first.lines1 <- "#!/bin/bash
#PBS -S /bin/bash
#PBS -l procs=1
#PBS -m bea
#PBS -M kgilbert@zoology.ubc.ca"
		first.lines2 <- paste(c("#PBS -l walltime=", walltime), collapse="")	
		first.lines3 <- paste(c("#PBS -l pmem=", RAM), collapse="")	
		
		first.lines <- paste(c(first.lines1, first.lines2, first.lines3), collapse="\n")
		write(first.lines, file="FirstLinesTemp.txt")
		system(paste(c("cat FirstLinesTemp.txt GrexLine.txt > FirstLines.txt")))
	}
	if(westgrid.server=="bugaboo"){
		if(is.null(walltime)) walltime <- "500:00:00"
		first.lines1 <- "#!/bin/bash
#PBS -l procs=1
#PBS -m bea
#PBS -M kgilbert@zoology.ubc.ca"
		first.lines2 <- paste(c("#PBS -l walltime=", walltime), collapse="")	
		first.lines3 <- paste(c("#PBS -l pmem=", RAM), collapse="")	

		first.lines <- paste(c(first.lines1, first.lines2, first.lines3), collapse="\n")
		write(first.lines, file="FirstLines.txt")
	}
	if(westgrid.server=="orcinus"){
		if(is.null(walltime)) walltime <- "240:00:00"
		first.lines1 <- "#!/bin/bash
#PBS -l procs=1
#PBS -m bea
#PBS -M kgilbert@zoology.ubc.ca"
		first.lines2 <- paste(c("#PBS -l walltime=", walltime), collapse="")	
		first.lines3 <- paste(c("#PBS -l pmem=", RAM), collapse="")	

		first.lines <- paste(c(first.lines1, first.lines2, first.lines3), collapse="\n")
		write(first.lines, file="FirstLines.txt")
	}
	
	input.file.path <- paste(c(input.directory, filename), collapse="")
	directory.lines <- paste(c(nemo.directory, input.file.path), collapse=" ")
	write(directory.lines, file="DirectoryLines.txt")
		
	system(paste(c("cat FirstLines.txt MiddleLine.txt DirectoryLines.txt LastLine.txt > PBS_", filename, ".pbs"), collapse=""))
}



#'
#' Take an individual Nemo .ini file and make a Westgrid PBS file for either orcinus, bugaboo, or grex servers.
#'
#' @title Make Westgrid PBS files
#'
#'
#'  @param ini.file.directory The directory containing .ini files for which PBS scripts will be made. Any file in this directory with the expression ".ini" in its name will be used to make a PBS file.
#'
#'  @param RAM The amount of RAM the run will be allowed to use, there is no default, and must include gb, mb, or kb in specification, e.g. "6gb".
#'
#'  @param walltime Default is set at max for specified server (except bugaboo, which is set to default at 500 hours), but other lengths can be specified here in units "100:00:00" hours:minutes:seconds.
#'
#'  @param server Which westgrid server the PBS file is meant for. Options are "grex", "orcinus" or "bugaboo" currently.
#'
#'  @param input.directory The directory on the Westgrid server which holds the nemo .ini file. Defaults are set for each server.
#'
#'  @param nemo.directory The directory on the Westgrid server from which the nemo executable is called. Defaults are set for each server.
#'
#'  @return
#'
#'  Creates PBS files from a .ini file for running Nemo on Westgrid's Orcinus, Bugaboo, or Grex servers.
#'
#' @author Kimberly J Gilbert
#'
#' @references \href{http://nemo2.sourceforge.net/index.html}{Nemo} is created and maintained by Fred Guillaume. The manual and source files are available online.
#'
#' @export multi.pbs


multi.pbs <- function(ini.file.directory, RAM, walltime, server, input.directory=NULL, nemo.directory=NULL){
	
	setwd(ini.file.directory)
	files.in.folder <- system("ls", intern=TRUE)
	# take only files with .ini in the name
	ini.files <- files.in.folder[grep(".ini", files.in.folder)]
	
	for(i in 1:length(ini.files)){
		make.pbs(filename=ini.files[i], 
		RAM=RAM, 
		walltime=walltime, 
		westgrid.server=server, 
		input.directory=input.directory, 
		nemo.directory=nemo.directory)
	}
}


