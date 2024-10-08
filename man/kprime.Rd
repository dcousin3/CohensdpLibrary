% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kprime.R
\name{kprime}
\alias{kprime}
\alias{dkprime}
\alias{pkprime}
\alias{qkprime}
\title{The K' or non-central K distribution.}
\usage{
pkprime(x, nu1, nu2, ncp) 
dkprime(x, nu1, nu2, ncp) 
qkprime(p, nu1, nu2, ncp)
}
\arguments{
\item{x}{the value from which a probability is sought;}

\item{nu1}{the first degree of freedom;}

\item{nu2}{the second degree of freedom;}

\item{ncp}{the noncentrality parameter;}

\item{p}{the probability from which a quantile is requested;}
}
\value{
\if{html}{\out{<div class="sourceCode">}}\preformatted{     The probability or quantile of a K' distribution.
}\if{html}{\out{</div>}}
}
\description{
The K' distribution was created to solve various problems in linear model statistics.
pkprime returns the cumulative probability of the lambda prime distribution with
parameters nu1, nu2, ncp; dkprime returns its density and qkprime, a quantile.
\insertCite{l99,pl10;textual}{CohensdpLibrary}.
}
\details{
kprime is a (p,d,q) set of functions that compute the K-prime distribution. This distribution has many applications,
including to obtain the sampling distribution of \emph{r} given a population rho and the predictive distributions
of rho given a sample \emph{r}. See \insertCite{l99,pl10;textual}{CohensdpLibrary}.

These functions are herein implemented from the FORTRAN source code of \insertCite{pl10b;textual}{CohensdpLibrary}.
Note that the library \emph{sadists} also implements this distribution \insertCite{p20}{CohensdpLibrary}.
However, the sadists::kprime distribution is inaccurate for small nu1 or small nu2.
}
\examples{

dkprime(11.1, 9, 8, 10.0)  # 0.09410193
pkprime(11.1, 9, 8, 10.0)  # 0.606652
qkprime(0.01, 9, 8, 10.0)  # 3.875234

}
\references{
\insertAllCited{}
}
