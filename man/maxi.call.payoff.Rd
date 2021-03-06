% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/payoffs.R
\name{maxi.call.payoff}
\alias{maxi.call.payoff}
\title{Max Call payoff}
\usage{
maxi.call.payoff(x, model)
}
\arguments{
\item{x}{Matrix of asset prices with \code{model$dim} columns and N rows}

\item{model}{list containing model params. Uses \code{model$K} as the Put strike.}
}
\description{
Multivariate max call
}
