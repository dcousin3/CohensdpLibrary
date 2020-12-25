.onLoad <- function(libname, pkgname) {
    
    # set the default method:
    options(CohenD.MAXITER = 10000, CohenD.TOLERAN = 0.000001)

    # load the external dynamically-link library CohenDLibrary
    #dyn.load("CohenDLibrary") #loaded automatically by R?

}

.onUnload <- function(libpath) {

    # unload the dll
    #dyn.unload("CohenDLibrary")

    # other cleaning to do?

}