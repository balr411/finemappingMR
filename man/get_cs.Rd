% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_cs.R
\name{get_cs}
\alias{get_cs}
\title{Get credible sets for exposure or outcome}
\usage{
get_cs(V, alpha, G_t_G)
}
\arguments{
\item{V}{The L x M matrix of prior variance estimates}

\item{alpha}{The list of length M of L x p matrices of pips for each region}

\item{G_t_G}{The list of length M of covariance matrices for the genotypes}
}
\value{
A list of lists containing the credible sets for each region
}
\description{
Function that gets credible sets for the exposure or outcome.
Note that much of the code is directly copied from the susieR package, but
modified so that a SuSiE object is not needed
}
