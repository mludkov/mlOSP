##########################################
#'  Arithmetic basket Put for d-dim x
#'  @title American Put payoff
#'
#' @param K is the strike
#' @param x is a vector of asset prices
#' @details in more than 1D, the prices are averaged and maxed with zero.
#' @export
put.payoff <- function(x,model)
{
   return( pmax(model$K-apply(x,1,mean),0))
}

########################
#' Call
#' @title Arithmetic average Call payoff
#' @details arithmetic basket for d-dim x)
#' @export
#' @inheritParams put.payoff
call.payoff <- function(x,model)
{
   return( pmax(apply(x,1,mean)-model$K,0))
}


########################
#' Digital Put
#' @title geometric digital Put payoff
#'
#' @export
#' @inheritParams put.payoff
digital.put.payoff <- function(x,model)
{
   return ( as.numeric(model$K > apply(x,1,prod)))
}

######################
#' Geometric Put
#' @title geometric basket Put
#' @export
#' @inheritParams put.payoff
geom.put.payoff <- function(x,model)
{
   return( pmax(model$K-apply(x,1,prod),0))
}


######################
#' Multivariate min Put
#' @title Min Put payoff
#'
#' @param x: asset prices
#' @export
#' @inheritParams put.payoff
mini.put.payoff <- function(x,model)
{
   return( pmax(model$K-apply(x,1,min),0))
}

######################
#' Multivariate max call
#' @title Max Call payoff
#'
#' @param x: asset prices
#' @export
#' @inheritParams put.payoff
maxi.call.payoff <- function(x, model)
{
   return( pmax(apply(x,1,max)-model$K,0))
}


######################
#' Put payoff for a stoch vol model
#' @title Basket Put with stochastic volatility
#' @param x: First coordinate is the asset price, the rest are vol factors
#' @export
#' @inheritParams put.payoff
#' @details Uses K-x[,1], ignores other coordinates
sv.put.payoff <- function(x, model)
{
    return (pmax(model$K-x[,1],0))
}
