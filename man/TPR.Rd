% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{TPR}
\alias{TPR}
\title{Description of the function}
\usage{
TPR(b_hat, b)
}
\arguments{
\item{b_hat}{an estimation of this vector}

\item{b}{a vector with non null values}
}
\value{
The proportion of non nulle values of B found non null in B_hat
}
\description{
Description of the function
}
\examples{
B <- c(1, -1, 1.5, 1.5, rep(0, 6), 2, 0, 2, 0)
B_hat <- c(0, 0, 1, 1, rep(0,4) , rep(1 , 0))
TPR(b_hat = B_hat , b =B)
}
