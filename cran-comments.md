## Submission

* Updated version 0.5.9

## Comment

* Just re-vamp fortran code for standard fortran 2008, following Prof Ripley's message.

## Test environments

* local WIN-64x install, R 4.3.0

* win-builder devel 

* r-hub.io


## R CMD check results

* There were no ERRORs and no WARNINGs.

* There was a few NOTEs on remote servers: 

* Win-dev returned zero note:

* R-hub returned two notes on some servers:

  - nothing on Debian Linux, R-devel, GCC ASAN/UBSAN
  
  - Ubuntu Linux 20.04.1 LTS, R-release, GCC says
      - * checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
Skipping checking math rendering: package 'V8' unavailable

  - Fedora Linux, R-devel, clang, gfortran says the same as Ubuntu

## Downstream dependencies

* This package does not rely on any other packages

* It requires Rtools4.3

