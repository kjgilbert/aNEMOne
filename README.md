# aNEMOne

![aNEMOne logo](https://github.com/kjgilbert/aNEMOne/raw/master/extra/SeaAnemone.png)

An R package to create input files for and analyze output files from [Nemo](http://nemo2.sourceforge.net/). 

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