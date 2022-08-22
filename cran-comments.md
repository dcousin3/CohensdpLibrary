## Submission

* Updated version 0.5.8

## Comment

* Man, thanks for seeing this. Corrected ier type to integer.

## Test environments

* local WIN-64x install, R 4.1.0

* win-builder devel 

* r-hub.io


## R CMD check results

* There were no ERRORs and no WARNINGs.

* There was a few NOTEs on remote servers: 

* Win-dev returned one note:

> checking CRAN incoming feasibility ... NOTE
Maintainer: 'Denis Cousineau <denis.cousineau@uottawa.ca>'

* R-hub returned four notes:

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
I do not know what this message means or how it can be resolved.



## Downstream dependencies

* This package does not rely on any other packages

* It requires Rtools4

