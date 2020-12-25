#' @useDynLib CohenDLibrary
#
#' @details
#' CohenDLibrary is a library that computes the Cohen's d along with confidence intervals. 
#' Its main contribution is to allow confidence intervals in the repeated-measure designs.
#' The main function is 
#' 
#' CohenD(statistics = list(...), method, design, conf.level, verbose = FALSE)  
#' 
#' see help(CohenD) for more
#' 
#' Subsidiary functions are
#' plambdap(x, nu, ncp) cumulative probability of the lambda prime distribution with parameters nu, ncp
#' dlambdap(x, nu, ncp) density
#' qlambdap(x, nu, ncp) quantile
#'  
#' pkprime(x, nu1, nu2, ncp) cumulative probability of the K prime distribution with parameters nu1, nu2, ncp
#' dkprime(x, nu1, nu2, ncp) density
#' qkprime(x, nu1, nu2, ncp) quantile
#'  
#' These functions are implemented from the FORTRAN source of Poitevineau & Lecoutre, 2010.
#' Note that the library sadists also has implementations for these two distributions. However,
#' the sadists::kprime distribution is innacurate for small nu1 or small nu2.
#' 
#' @keywords internal
"_PACKAGE"
#> [1] "_PACKAGE"

#' @title CohenDLibrary
#'
#' @description
#' CohenDLibrary provides the main command CohenD to compute Cohen's d (noted d_p) and its 
#' confidence interval in either within-subject, between-subject design
#' and single-group design. For the between-subject design, MBESS already 
#' has an implementation based on the "pivotal" method. This method is exact for between-subject  
#' design but slower to compute than the "lambdaprime" method.
#'  
#'
#' 
#'
#' @author Denis Cousineau, \email{denis.cousineau@@uottawa.ca}
#' @references \url{https://.../...}
#' @keywords Measurement precision; rounding
#'
#' @examples
#' # define a vector (it could be a 1-colum matrix or a one-column data.frame)
#' x1 <- c(3,4,5)
#' x2 <- c(6,7,8,9)
#'
#' # get the cohen's d with the confidence interval
#' CohenD(statistics=list(m1=mean(x1), m2=mean(x2), sd1=sd(x1), sd2=sd(x2), n1=len(x1), n2=len(x2)))
#' CohenD(statistics=list(m1=15, m2=20, sd1=4, sd2=4, n1=25, n2=25), design="between", conf.level=.80)
#'
#'
#' # By default, the functions assumes maximum iterations and tolerance values.
#' # You can change these defaults with the options:
#' options(CohenD.MAXITER = 10000)    # this is the default maximum iterations
#' options(CohenD.TOLERAN = 0.000001) # this is the default tolerance value
#' # Increasing MAXITER or decreasing TOLERAN slows computations.
#'
