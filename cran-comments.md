## Submission

* First version 0.5


## Comment

* Win-dev returned one note:

> checking CRAN incoming feasibility ... NOTE
Maintainer: 'Denis Cousineau <denis.cousineau@uottawa.ca>'

* Corrected error with debian compiler:   array over-run in .Fortran("sublprimepdf") in integer argument 6

* R-hub returned three notes:

> checking CRAN incoming feasibility ... NOTE
Maintainer: 'Denis Cousineau <denis.cousineau@uottawa.ca>'
I do not know what this message means or how it can be resolved.

> checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: 'Rdpack'
    All declared Imports should be used.
Rdpack is for building documentation but is not required in the present package

> checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'
This directory is created by r-hub but is not part of the present package.

> checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found

## Test environments

* local WIN-64x install, R 4.1.0

* win-builder devel 

* win-builder release

* r-hub.io


## R CMD check results

* There were no ERRORs and no WARNINGs.

* There was 0 NOTE: 


## Downstream dependencies

* This package does not rely on any other packages

* It requires Rtools4

