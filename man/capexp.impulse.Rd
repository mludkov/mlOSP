% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mlOSP_utils.R
\name{capexp.impulse}
\alias{capexp.impulse}
\title{Compute intervention function for Capacity Expansion impulse problems}
\usage{
capexp.impulse(cur_x, model, fit, ext = FALSE)
}
\arguments{
\item{cur_x}{Set of inputs where to compute the intervention function
Should be a n x 2 matrix, with first column for prices and second column
for capacities. Impulse affects second column only.}

\item{model}{a list containing all model parameters. 
In particular must have \code{model$imp.cost.capexp} to compute cost of impulses}

\item{fit}{Object containing the one-step-ahead functional approximator for V(k,x)}

\item{ext}{logical flag (default is FALSE) whether to return extended information}
}
\description{
Compute intervention function for Capacity Expansion impulse problems
}
\details{
Calculates the intervention operator for a 2-D capacity
expansion problem. This is done by running \code{optimize} on the
cost-to-go based on \code{fit}. Calls \code{ospPredict}
}
