#' @title lambdap
#'
#' @details
#' lambdap are functions that compute the Lambda prime distribution. 
#'
#' @description
#' plambdap(x, nu, ncp) cumulative probability of the lambda prime distribution with parameters nu, ncp
#' dlambdap(x, nu, ncp) density of the lambda prime distribution
#' qlambdap(x, nu, ncp) quantile of the lambda prime distribution
#'  
#' These functions are implemented from the FORTRAN source of Poitevineau & Lecoutre, 2010.
#' Note that the library sadists also has implementations for these two distributions. However,
#' the sadists::kprime distribution is innacurate for small nu1 or small nu2.
#'
#' @examples
#' 
#' dlambdap(6.0, 99, 6.19) # 0.436873
#' plambdap(6.0, 99, 6.19) # 0.3607524
#' qlambdap(0.55, 99, 6.19) # 6.31081
#' 

#' @export
plambdap <- function( x, nu, ncp) {
    res <- .Fortran("sublprimepdf",
        as.double(x), 
        as.double(nu), 
        as.double(ncp), 
        as.double(getOption("CohenD.TOLERAN")), 
        as.integer(getOption("CohenD.MAXITER")),
        as.integer(-99),
        as.double(0.00)  )

    if (res[[6]] != 0) {
        warning(paste("The subroutine signals convergence problems ", res[[6]] ))
    }
    return( res[[7]] )
}

#' @export
qlambdap <- function( x, nu, ncp) {
    res <- .Fortran("sublprimeidf",
        as.double(x), 
        as.double(nu), 
        as.double(ncp), 
        as.double(getOption("CohenD.TOLERAN")), 
        as.integer(getOption("CohenD.MAXITER")),
        as.integer(-99),
        as.double(0.00)  )

    if (res[[6]] != 0) {
        warning(paste("The subroutine signals convergence problems ", res[[6]] ))
    }
    return( res[[7]] )
}

#' @export
dlambdap <- function( x, nu, ncp) {
    res <- .Fortran("sublprimecdf",
        as.double(x), 
        as.double(nu), 
        as.double(ncp), 
        as.double(getOption("CohenD.TOLERAN")), 
        as.integer(getOption("CohenD.MAXITER")),
        as.integer(-99),
        as.double(0.00)  )

    if (res[[6]] != 0) {
        warning(paste("The subroutine signals convergence problems ", res[[6]] ))
    }
    return( res[[7]] )
}
