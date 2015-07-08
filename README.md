# aNEMOne

[![Build Status](https://travis-ci.org/kjgilbert/aNEMOne.png?branch=master)](https://travis-ci.org/kjgilbert/aNEMOne) 

An R package to create input files for and analyze output files from [Nemo](http://nemo2.sourceforge.net/). 

![aNEMOne logo](https://github.com/kjgilbert/aNEMOne/raw/master/extra/SeaAnemones.jpg)

## Documentation

Nemo is an individual-based forward-time simulation program created and maintained by Fred Guillaume. The latest release is available for download [here](http://nemo2.sourceforge.net/).

This package is meant to aid in writing input files for Nemo, particularly when dealing with large landscapes consisting of many patches that require large range-wide matrices for specifying dispersal distances or distances within which breeding can happen from a focal patch. See additions to the Nemo code [documented here](https://github.com/kjgilbert/NemoDispersalKernel) for some of the capabilities that this packages is meant to deal with. Future updates will include functions for analyzing stat file outputs and other output files.

The package can be installed using [devtools](https://github.com/hadley/devtools), which itself can be installed from CRAN with

```
install.packages("devtools")
```

Once devtools is installed, run

```
library(devtools)
install_github("kjgilbert/aNEMOne")
library(aNEMOne)
```
and the package will be installed and open.

#### Functions

An [example script](https://github.com/kjgilbert/aNEMOne/blob/master/extra/Example_MakeInputs.R) shows how this package's functions may be used to produce a .ini input file for Nemo.

The main function `make.input` takes all of the parameters that would normally appear in your input file with some defaults already set. See `?make.input` for details. This function really provides no advantage unless you use the additional functions to create large matrices and kernels that can now be automatically inserted into the input file without worry about copy and paste or bracket errors.

`make.delet.input` does the same as above except has deleterious loci instead of neutral loci.

##### Miscellaneous functions used for making inputs

See `?patch.cap` for details on making an array of patch capacities. Currently only two different *K*'s may be set across the landscape, but more may easily be added, so feel free to request so. This function also works for male and female specific carrying capacities as well as temporally set capacities.

See `?make.kernel.and.matrix` for details on making connectivity matrices and probability kernels for dispersal and breeding windows. The function returns the array for the kernel and prints the connectivity matrix to file.

##### Selection landscape

See `?make.landscape` for details on making the landscape of phenotypic optima for Nemo's `selection_local_optima` input parameter. The function returns the mean optimum value for the first column of landscape patches, visualizes the landscape, and write to file the matrix to be fed in for input.

See `?step.landscape` for details on making a stepwise gradient landscape of phenotypic optima for Nemo's `selection_local_optima` input parameter. The function returns the optimum value for the first column of landscape patches, visualizes the landscape, and write to file the matrix to be fed in for input.

See `?changing.landscape` for details on how to combine landscapes of two types or steepnesses.

##### Useful for running analyses on [Westgrid](https://www.westgrid.ca/)

`make.pbs` creates a PBS file for the provided .ini file. It is important that [GrexLine.txt, MiddleLine.txt, and LastLine.txt](https://github.com/kjgilbert/aNEMOne/tree/master/extra) are in the same directory as the .ini file (and the same holds for multi.pbs below).

`multi.pbs` creates PBS files for all .ini files in a directory.

##### Handle outputs

`dist.delet.effects` Plots histograms of the effect sizes of deleterious mutations from a .del file.

`sim.results.pop` IN PROGRESS (SimResults_Pop.R)
