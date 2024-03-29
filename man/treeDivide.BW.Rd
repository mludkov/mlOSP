% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlOSP_utils.R
\name{treeDivide.BW}
\alias{treeDivide.BW}
\title{Create a Bouchard-Warin equi-probable grid of regression sub-domains}
\usage{
treeDivide.BW(grid, curDim, model, test.paths)
}
\arguments{
\item{grid}{dataset of x-values}

\item{curDim}{dimension of the grid}

\item{model}{a list containing all model parameters. In particular must have 
\code{model$nBins} defined}

\item{test.paths}{testing paths to predict along as well}
}
\description{
Create a Bouchard-Warin equi-probable grid of regression sub-domains
}
\details{
Recursively sort along each of the d-coordinates
At the end do local linear regression at each leaf
This is a recursive algorithm!
first column is reserved for the y-coordinate (timingValue)
It's safest if nrows(grid) is divisible by \code{model$nBinsmodel$dim}
}
