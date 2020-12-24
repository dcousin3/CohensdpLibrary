# Cohen's D library
Getting the Cohen's D and its confidence interval

This repository contains source code and a package suitable for R 3.0 and above. It computes 
the Cohen's D and its confidence intervals in within-subject and between-subject designs.

See :
Cousineau, D. (under preparation) n.d. 
for the formal mathematical derivations of the results.

Note that for between-subject designs, the library MBESS already has a potent function to get confidence intervals (Kelley, 2019).

You can install this library on you computer if the library devtools is installed with:

	devtools::install_github("dcousin3/CohenDLibrary")<br>
	library(CohenDLibrary)

Check the <ToCome>.pdf document for more on the functions. The main function is CohensD so help(CohenD) will get you started.


