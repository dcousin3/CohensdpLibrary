#' @title Cohen's standardized mean difference.
#'
#' @aliases Cohensdp
#'
#' @description
#' Cohensdp computes the Cohen's d (noted d_p) and its confidence interval in 
#' either within-subject, between-subject design and single-group design. For 
#' the between-subject design, MBESS already has an implementation based on the 
#' "pivotal" method. This method is exact for between-subject design but slower 
#' to compute than the currently implemented  method based on the Lambda prime
#' distribution. See \insertCite{h81,c22a,c22b,gc18;textual}{CohensdpLibrary}.
#'
#' @usage
#' Cohensdp(statistics, design, gamma )
#'
#' @param statistics    a list of pre-computed statistics. The statistics to provide 
#'                      depend on the design:
#'                      - for "between": m1, m2 the means of the two groups, s1, s2 
#'                              the standard deviation of the two groups, and n1, n2, 
#'                              the sample sizes of the two groups;
#'                        - for "within": m1, m2, s1, s2, n, and r or rho the 
#'                          correlation between the measure; 
#'                        - for "single": m, s, n and mu the reference mean from which 
#'                          m is standardized).
#' @param design        the design of the measures ("within", "between", or "single");
#' @param gamma         the confidence level of the confidence interval (default 0.95) 
#'
#' @return            the Cohen's d_p statistic and its confidence interval
#'
#' @details
#' This function is currently unavailable in within-subject design when the population 
#' rho is unknown.
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @examples
#'
#' # example in which the means are 114 vs. 101 with sds of 14.3 and 12.5 respectively
#' Cohensdp( statistics = list( m1= 101, m2= 114, s1= 12.5, s2= 14.3, n1= 12, n2= 12 ), 
#'           design     = "between")
#'
#' # example in a repeated-measure design
#' Cohensdp( statistics = list( m1= 101, m2= 114, s1= 12.5, s2= 14.3, n= 12, r= 0.53 ), 
#'           design     = "within")
#'
#' # example with a single-group design where mu is a population parameter
#' Cohensdp( statistics = list( m = 101, m0 = 114, s = 12.5, n = 10 ), 
#'           design     = "single")
#'
#' # The results can be displayed in three modes
#' res <- Cohensdp( statistics = list( m = 101, m0 = 114, s = 12.5, n = 10), 
#'                  design     = "single")
#' res              # a raw result of the Cohen's d_p and its confidence interval
#' summarize( res ) # a human-readable output
#' explain( res )   # a human-readable ouptut with additional explanations.
#'
#' 

#' @export 
Cohensdp <- function( 
            statistics = NULL, 
            design     = NULL, 
            gamma      = 0.95
    ) {

    ##############################################################################
    # STEP 1: Input validation
    ##############################################################################
    # 1.1: check that statistics is a non-empty list
    if( !(is.list(statistics)) )   stop( messageSnotL() ) 
    if( length(statistics) == 0 )  stop( messageSnotE() ) 

    # 1.2: check that the designs are legitimate
    if (is.null(design)) stop(messageDempt)
    if (!(design %in% c("within","between","single"))) 
        stop( messageDnotG(design) )

    # 1.3: check that the confidence level gamma is legitimate
    if (gamma < 0 | gamma > 1) 
        stop( messageGnotG(gamma)  )


    ##############################################################################
    # STEP 2: let's do the computation and return everything
    ##############################################################################
    fct = paste("Cohensdp", design, sep=".")
    res = lapply(list(statistics), fct, gamma = gamma)[[1]]

    # preserve everything in an object of class CohensdpObject
    dpObject = list(
        type       = "dp",
        estimate   = res[1],
        interval   = c(res[2], res[3]),
        statistics = statistics, 
        design     = design, 
        gamma      = gamma
    )
    class(dpObject) <- c("CohensdpObject", class(dpObject) )
    return( dpObject )
}





##############################################################################
# DEFINITIONS of all the designs subfunctions
# there are 3 functions in this package:
#
# within
#    using exact = within using lsecond mixture (Cousineau, 2022)
#    NOT IMPLEMENTED
#        Algina & Kesselman, 2003 -- see MBESS implementation (slow)
#        Morris, 2000, Goulet-Pelletier & Cousineau, 2018, Adjusted Lprime;
#    See Cousineau & Goulet-Pelletier, 2020, TQMP, for a review.
#
# between 
#    using exact = between using lprime(Lecoutre, 1999, 2007)
#    NOT IMPLEMENTED
#        Steiger & Fouladi, 1997 -- see MBESS implementation (slow)
#        Bayes, 1763, Hedges, 1981; 
#    See Cousineau & Goulet-Pelletier, 2020, PsyArXiv, for a review.)
#
# single
#    using lprime    (Cousineau, 2022)
#    NOT IMPLEMENTED
#        using pivotal -- see MBESS implementation (slow)
#
##############################################################################

##############################################################################
##### single #################################################################
##############################################################################
Cohensdp.single <- function(statistics, gamma = .95) {
    sts  <- vfyStat(statistics, c("m","m0","s","n"))

    #get pairwise statistics Delta means and pooled SD
    dmn  <- sts$m - sts$m0
    #compute biased Cohen's d_p
    dp   <- dmn / sts$s  

    dlow <- qlprime(1/2-gamma/2, nu = sts$n-1, ncp = dp * sqrt(sts$n) ) 
    dhig <- qlprime(1/2+gamma/2, nu = sts$n-1, ncp = dp * sqrt(sts$n) ) 

    limits <- c(dlow, dhig) / sqrt(sts$n)
    c(dp, limits)

}


##############################################################################
##### between ################################################################
##############################################################################
Cohensdp.between <- function(statistics, gamma = .95) {
    sts  <- vfyStat(statistics, c("m1","s1","n1","m2","s2","n2"))

    #get pairwise statistics Delta means and pooled SD
    dmn  <- sts$m1 - sts$m2
    sdp  <- sqrt((sts$s1^2 + sts$s2^2)/2)
    n    <- 2 / (1/sts$n1 + 1/sts$n2)           #harmonic mean
    #compute biased Cohen's d_p
    dp   <- dmn / sdp  

    dlow <- qlprime(1/2-gamma/2, nu = 2*(n-1), ncp = dp * sqrt(n/2) ) 
    dhig <- qlprime(1/2+gamma/2, nu = 2*(n-1), ncp = dp * sqrt(n/2) ) 

    limits <- c(dlow, dhig) / sqrt(n/2)
    c(dp, limits)

}


##############################################################################
##### within #################################################################
##############################################################################
Cohensdp.within <- function(statistics, gamma = .95) {
    if ("r" %in% names(statistics)) 
        res <- Cohensdp.within.rhounknown( statistics, gamma )
    else
        res <- Cohensdp.within.rhoknown( statistics, gamma )
    res
}

Cohensdp.within.rhounknown <- function(statistics, gamma = .95) {
    sts  <- vfyStat(statistics, c("m1","s1","m2","s2","n", "r"))

    #get pairwise statistics Delta means and pooled SD
    dmn  <- sts$m1 - sts$m2
    sdp  <- sqrt((sts$s1^2 + sts$s2^2)/2)
    W    <- (sts$s1 * sts$s2) / ((sts$s1^2 + sts$s2^2)/2)

    #compute biased Cohen's d_p 
    dp   <- dmn / sdp  

    if( getOption("CohensdpLibrary.SHOWWARNINGS") )
        warning("No known confidence interval at this time when rho is unknown. Using Adjusted lambda' method...")

    W    <- (sts$s1 * sts$s2) / ((sts$s1^2 + sts$s2^2)/2)
    rW <- sts$r * W

    lambda <- dp * J.single( sts ) * sqrt(sts$n/(2*(1-rW)))

    #quantile of the noncentral t distribution 
    dlow = qlprime(1/2-gamma/2, nu = 2/(1+sts$r^2)*(sts$n-1), ncp = lambda )
    dhig = qlprime(1/2+gamma/2, nu = 2/(1+sts$r^2)*(sts$n-1), ncp = lambda )

    limits <- c(dlow, dhig) / sqrt(sts$n/(2*(1-rW))) / J.single( list( n=2*(sts$n-1)/(1+sts$r^2)+1 ))
    c(dp, limits)

}

Cohensdp.within.rhoknown <- function(statistics, gamma = .95) {
    sts  <- vfyStat(statistics, c("m1","s1","m2","s2","n", "rho"))

    #get pairwise statistics Delta means and pooled SD
    dmn  <- sts$m1 - sts$m2
    sdp  <- sqrt((sts$s1^2 + sts$s2^2)/2)

    #compute biased Cohen's d_p 
    dp   <- dmn / sdp  

    #quantile of the (noncentral, nonstandard) lambda'' distribution 
    dlow = qlsecond(1/2-gamma/2., n = sts$n, d = dp, rho = sts$rho )
    dhig = qlsecond(1/2+gamma/2., n = sts$n, d = dp, rho = sts$rho )

    limits <- c(dlow, dhig) 
    c(dp, limits)

}


##############################################################################
##### This is it!#############################################################
##############################################################################
