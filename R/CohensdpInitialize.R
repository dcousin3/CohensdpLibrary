# things to initialize
.onLoad <- function(libname, pkgname) {    
    # set the default arguments for the iterative functions:
    options(
        "CohensdpLibrary.MAXITER" = 32500,     # this is the maximum in short integer; 
                                               # that should be more than enough!
                                               # The nbre of iterations only rarely exceeds 2000
        "CohensdpLibrary.TOLERAN" = 0.0000001, # less than 7 decimals in the additional steps...
        "CohensdpLibrary.FORMAT"  = "%5.3f",   # printing results with 3 decimals 
                                               # Should be more than enough!
                                               # See Cousineau, 2020, JMP, for the number of decimals
        "CohensdpLibrary.SHOWWARNINGS" = TRUE  # use to inhibit messages and warnings
    )

    # load the external dynamically-link library Cohensdp
    # dyn.load("CohensdpLibrary") #loaded automatically by R

}

# things to clean?
.onUnload <- function(libpath) {
    # unload the dll
    # dyn.unload("CohensdpLibrary") #unloaded automatically by R

}


# define the error messages here as they are used in many functions
messageSnotL <- function() {"Argument `statistics` is missing or is not a list. Exiting..."}
messageSnotE <- function() {"Argument `statistics` is empty. Exiting..." }
messageDnotG <- function(design) {paste("Not a known `design` \"",design,"\". Use within, between, single. Exiting...", sep="") }
messageDempt <- function() {"Argument `design` not given. Exiting..."}
messageGnotG <- function(gamma) {paste("The confidence level `gamma` \"",gamma,"\" is not between 0 and 1. Exiting...", sep="")}
messageSinct <- function(statname) {
                    paste("The list of statistics is incomplete. Are needed: ", 
                        paste(statname, collapse = " "), sep = "")}
messageWier  <- function( fctname, ier) {paste("The subroutine ",fctname," signals convergence problems ", ier, sep="" ) }
messageNtsm  <- function(n) {paste("Sample size ",n," too small. Exiting...",sep="")}
messageSneg  <- function(s) {paste("Sample standard deviation ",s," cannot be negative. Exiting...",sep="")}
messageRwrg  <- function(r) {paste("Correlation ",r," must be between -1 and +1. Exiting...",sep="")}


# verify that the named statistics are in the list or else issue an error message
vfyStat <- function(statlist, statname) {
    # check that the required statistics were provided
    if (!(all(statname %in% names(statlist))))
        stop( messageSinct(statname) )

    # check the domain of the statistics
    if ("n" %in% names(statlist) )  {if (statlist$n  <2) stop(messageNtsm(statlist$n)) }
    if ("n1" %in% names(statlist) ) {if (statlist$n1 <2) stop(messageNtsm(statlist$n1)) }
    if ("n2" %in% names(statlist) ) {if (statlist$n2 <2) stop(messageNtsm(statlist$n2)) }
    if ("s" %in% names(statlist) )  {if (statlist$s  <0) stop(messageSneg(statlist$s)) }
    if ("s1" %in% names(statlist) ) {if (statlist$s1 <0) stop(messageSneg(statlist$s1)) }
    if ("s2" %in% names(statlist) ) {if (statlist$s2 <0) stop(messageSneg(statlist$s2)) }
    if ("rho" %in% names(statlist) ){if ((statlist$rho < -1) | (statlist$rho > +1)) stop(messageRwrg(statlist$rho)) }
    if ("r" %in% names(statlist) )  {if (  (statlist$r < -1) |   (statlist$r > +1)) stop(messageRwrg(statlist$r)) }

    statlist
}

####################################################################
## new METHODS : here only router function and generic
####################################################################


#' @title explain 
#' 
#' @description
#' explain provides a human-readable, exhaustive, description of
#' the results. It also provides references to the key results.
#' 
#' @param x   an object to explain
#' @param ... ignored
#' 
#' @export
explain <- function(x, ...) {  UseMethod("explain") }

#' @export 
explain.default <- function(x, ...) { print(x) } 



#' @title summarize 
#' 
#' @description
#' Summarize provides a human-readable output of a dpObject. it is 
#' synonyme of summary (but as actions are verbs, I used a verb).
#' 
#' @param x   an object to summarize
#' @param ... ignored
#' 
#' @export 
summarize <- function(x, ...) {  UseMethod("summarize") }

#' @method summarize default 
#' @export 
summarize.default <- function(x, ...) { print(x) }



