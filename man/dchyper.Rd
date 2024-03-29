% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chyper.R
\name{dchyper}
\alias{dchyper}
\title{Probability mass function for conditional hypergeometric distributions}
\usage{
dchyper(k, s, n, m, verbose = T)
}
\arguments{
\item{k}{an integer or vector of integers representing the overlap size}

\item{s}{an integer representing the size of the intersecting population}

\item{n}{a vector of integers representing the sizes of each non-intersecting population}

\item{m}{a vector of integers representing the sample sizes}

\item{verbose}{T/F should intermediate messages be printed?}
}
\value{
The probability of sampling \code{k} of the same items in all samples
}
\description{
Calculates the PMF of a conditional hypergeometric distribution: the distribution of how many items are in the overlap of all samples when samples of arbitrary size are each taken without replacement from populations of arbitrary size.
}
\examples{
dchyper(c(3,5), 10, c(12,13,14), c(7,8,9))

}
