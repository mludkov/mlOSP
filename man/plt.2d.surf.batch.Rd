% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batchDesign_utils.R
\name{plt.2d.surf.batch}
\alias{plt.2d.surf.batch}
\title{Visualize and compare 2D emulator + stopping region for two fits}
\usage{
plt.2d.surf.batch(
  fit1,
  fit2,
  r1,
  r2,
  batch1,
  batch2,
  x = seq(31, 43, len = 201),
  y = seq(31, 43, len = 201)
)
}
\arguments{
\item{fit1, fit2}{can be any of the types supported by \code{\link{forward.sim.policy}}}

\item{r1, r2}{batch vectors for the two fits}

\item{batch1, batch2}{batch heristics for two fits; Passed to \code{ggplot()}
This only works for \code{km} and \code{het/homGP} objects}

\item{x, y}{locations to use for the \code{predict()} functions. Default is a 200x200 fine grid.
Passed to \code{expand.grid}}
}
\description{
two-dimensional image of contour + site + batch plot for two fits
}
