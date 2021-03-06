% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/payoffs.R
\name{digital.put.payoff}
\alias{digital.put.payoff}
\title{geometric digital Put payoff}
\usage{
digital.put.payoff(x, model)
}
\arguments{
\item{x}{is a matrix of asset prices of dimension N * \code{model$dim}}

\item{model}{list containing model params. Uses \code{model$K} as the Put strike.}
}
\description{
Digital Put
}
