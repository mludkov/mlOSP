##########################################
#'  Arithmetic basket Put for d-dim x
#'  @ title American Put payoff
#'
#' @param K is the strike
#' @param x is a vector of asset prices
#' @details in more than 1D, the prices are averaged and maxed with zero.
#' @export
put.payoff <- function(x,K=40)
{
   return( pmax(K-apply(x,1,mean),0))
}

#' @title American Call payoff
#' @details arithmetic basket for d-dim x)
#' @export
#' @inheritParams put.payoff
call.payoff <- function(x,K=40)
{
   return( pmax(apply(x,1,mean)-K,0))
}


########################
#' geometric digital Put payoff
#'
#' @export
#' @inheritParams put.payoff
digital.put.payoff <- function(x,K=10)
{
   return ( as.numeric(K > apply(x,1,prod)))
}

######################
#' geometric basket Put
#' @export
#' @inheritParams put.payoff
geom.put.payoff <- function(x,K=40)
{
   return( pmax(K-apply(x,1,prod),0))
}


######################
#' Min Put payoff
#'
#' @param x: asset prices
#' @export
#' @inheritParams put.payoff
minPut <- function(x,K=40)
{
   return( pmax(K-apply(x,1,min),0))
}

######################
#' @title Max Call payoff
#'
#' @param x: asset prices
#' @export
#' @inheritParams put.payoff
maxCall <- function(x, K=100)
{
   return( pmax(apply(x,1,max)-K,0))
}


######################
#' Put payoff for a stoch vol model
#' @title Basket Put with stochastic volatility
#' @param x: First coordinate is the asset price, the rest are vol factors
#' @export
#' @inheritParams put.payoff
#' @details Uses K-x[,1], ignores other coordinates
sv.put <- function(x, K=40)
{
    return (pmax(K-x[,1],0))
}
