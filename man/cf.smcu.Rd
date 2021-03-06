% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlOSP_utils.R
\name{cf.smcu}
\alias{cf.smcu}
\title{Straddle Maximum Contour Uncertainty criterion}
\usage{
cf.smcu(objMean, objSd, gamma = 1.96)
}
\arguments{
\item{objMean}{predicted mean response}

\item{objSd}{posterior standard deviation of the response}

\item{gamma}{weight on the variance}
}
\description{
straddle MCU with a specified variance weight
}
\details{
compute the UCB criterion with constant weight: gamma*s(x) - |f(x)|
}
\references{
Mike Ludkovski, Kriging Metamodels and Experimental Design for Bermudan Option Pricing
 Journal of Computational Finance, 22(1), 37-77, 2018
 
 X.Lyu, M Binois, M. Ludkovski (2020+) Evaluating Gaussian Process Metamodels and Sequential Designs for 
 Noisy Level Set Estimation <https://arxiv.org/abs/1807.06712>
}
\seealso{
[osp.seq.design]
}
\author{
Mike Ludkovski
}
