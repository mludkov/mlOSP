% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlOSP_utils.R
\name{swing.policy}
\alias{swing.policy}
\title{Forward simulation of a swing payoff based on a sequence of emulators}
\usage{
swing.policy(
  x,
  M,
  fit,
  model,
  offset = 1,
  use.qv = FALSE,
  n.swing = 1,
  verbose = FALSE
)
}
\arguments{
\item{x}{a matrix of starting values (N x \code{model$dim}).
If input \code{x} is a list, then use the grids specified by x}

\item{M}{number of time steps to forward simulate}

\item{fit}{a list of fitted emulators that determine the stopping classifiers to be used}

\item{model}{List containing all model parameters. In particular uses \code{model$dt,model$r} 
for discounting and \code{model$swing.payoff} to compute payoffs}

\item{offset}{deprecated}

\item{use.qv}{experimental, do not use}

\item{n.swing}{number of swing rights (integer, at least 1)}

\item{verbose}{for debugging purposes}
}
\value{
a list containing:
\itemize{
 \item \code{payoff}: a vector of length `nrow(x)` containing the resulting payoffs NPV from $t=0$
 \item \code{tau} matrix of the times when stopped. Columns represent the rights exercised
 \item  \code{nsims} number of total 1-step simulations performed
}
}
\description{
Simulate \eqn{\sum_k h(X_{tau_k})} using \code{fit} emulators
}
\details{
Should be used in conjuction with the \code{\link[mlOSP]{swing.fixed.design}} function that builds the emulators.
}
