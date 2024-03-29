% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlOSP_utils.R
\name{treeDivide.BW.1d}
\alias{treeDivide.BW.1d}
\title{1-d version of \code{\link{treeDivide.BW}} that stores all the local fits}
\usage{
treeDivide.BW.1d(grid, curDim, model, test)
}
\arguments{
\item{grid}{dataset of x-values}

\item{curDim}{dimension of the grid}

\item{model}{a list containing all model parameters. In particular must have 
\code{model$nBins} defined}
}
\description{
1-d version of \code{\link{treeDivide.BW}} that stores all the local fits
}
