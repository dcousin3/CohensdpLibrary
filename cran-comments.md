## Submission

* Updated version 0.6.0

## Comment

* adjusted pooled standard deviation in between design

## Test environments

* local WIN-64x install, R 4.4.2
* win-builder devel 


## R CMD check results

* Local
    -There were no ERRORs and no WARNINGs.

* Win-dev returned:
  - Status: 1 NOTE
  - indicated that all the URL were possibly missing: I re-checked them all.
  - indicated that all the DOI were possibly missing: I re-checked them all.
  - indicated that title case was not respected because I wrote d_p. This is a mathematical symbol.
  
## Downstream dependencies

* This package does not rely on any other packages
* It requires Rtools4.3

