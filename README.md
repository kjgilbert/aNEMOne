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

**PACKAGE CREATION IS CURRENTLY STILL IN PROGRESS** The functions outlined here are those currently tested to work properly.

An [example script] shows how this package's functions may be used to produce a .ini input file for Nemo.

The main function `make.input` takes all of the parameters that would normally appear in your input file with some defaults already set. See `?make.input` for details. This function really provides no advantage unless you use the additional functions to create large matrices and kernels that can now be automatically inserted into the input file without worry about copy and paste or bracket errors.

See `?patch.cap` for details on making an array of patch capacities. Currently only two different *K*'s may be set across the landscape, but more may easily be added, so feel free to request so. This function also works for male and female specific carrying capacities as well as temporally set capacities.

See `?make.kernel.and.matrix` for details on making connectivity matrices and probability kernels for dispersal and breeding windows. The function returns the array for the kernel and prints the connectivity matrix to file.

See `?make.landscape` for details on making the landscape of phenotypic optima for Nemo's `selection_local_optima` input parameter. The function returns the mean optimum value for the first column of landscape patches, visualizes the landscape, and write to file the matrix to be fed in for input.