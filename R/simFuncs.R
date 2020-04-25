
####################################
#' Simulate paths of Geometric Brownian Motion with constant parameters
#'
#' Simulate from \eqn{p(X_t|X_{t-1})}
#' Use log-normal transition density specified by the model
#' @param x0 is the starting values (vector)
#' @param dt is the step size
#' @param model contains all the other parameters, including volatility \code{model$sigma}
#' interest rate \code{model$r} and continuous dividend yield \code{model$div}
#' @export
#' @md
sim.gbm <- function( x0, model, dt=model$dt)
{
    len <- nrow(x0)

    newX <- x0
    for (j in 1:model$dim)
       newX[,j] <- x0[,j,drop=F]*exp( rnorm(len)*model$sigma[j]*sqrt(dt) +
            (model$r- model$div- model$sigma[j]^2/2)*dt)

    return (newX)
}

####################################
#' Simulate from exponential Ornstein Uhlenbeck process
#' @inheritParams sim.gbm
#' @export
sim.ouExp <- function( x0, model, dt=model$dt)
{
    len <- nrow(x0)
    newX <- log(x0)
    for (j in 1:model$dim) {
       newX[,j] <- newX[,j]*exp( - model$alpha*dt) + (1-exp(-model$alpha*dt))*(log(model$meanRev) -
       model$sigma^2/2/model$alpha)
       newX[,j] <- exp( newX[,j] + rnorm(len)*model$sigma[j]*sqrt( (1-exp(-2*model$alpha*dt))/2/model$alpha)
       )
    }
    return (newX)
}

####################################
#' Simulate from 1-D discretized exponential Ornstein Uhlenbeck process
#' See Bender (SIFIN 2011). Only works in one dimension
#' @inheritParams sim.gbm
#' @export
#' @details Requires
#' \itemize{
#' \item  \code{model$rho} -- similar to mean-reversion rate, should be close to 1
#' \item \code{model$mu} -- mean-reversion level
#' \item \code{model$sigma} -- volatility
#' }
sim.logOU_Discrete <- function( x0, model, dt=model$dt)
{
  len <- nrow(x0)
  newX <- log(x0)
  newX <- (1-model$rho)*(newX - model$mu) + model$mu + model$sigma*rnorm(len)
  newX <- exp(newX)
  
  return (newX)
}




########################################
#' Simulate an exp-OU stoch volatility model with 1 or 2 vol factors
#' @inheritParams sim.gbm
#' @export
#' @param x0 should have 2 or 3 columns
#' @section Usage:  Need the following fields in model: \code{svMean} (mean-rev level), \code{svAlpha} (men-rev strength),
#' svEpsY (fast scaling param), svVol (vol-vol), svRho (corr)
#' for 2-factor also need: \code{svMeanZ} (slow scale mean-reversion level), \code{svAlphaZ} (mean-reversion strength), \code{svDeltaZ} (slow
#' scaling parameter)
#' \code{svVolZ} (Z volatility), \code{svRhoZ} (correlation between Z and S), \code{svRhoYZ} ( correlation between the 2 factors)
sim.expOU.sv <- function(x0, model, dt=model$dt,useEuler=F)
{
    len <- nrow(x0)
    newX <- x0
    rhoBar <- sqrt(1- model$svRho^2)

    t.seq <- seq(0,dt,by=model$eulerDt)
    if (max(t.seq) < dt)
       t.seq <- c(t.seq,dt)
    step.seq <- diff(t.seq)
    nSteps <- length(step.seq)

    #if (ncol(x0) == 2)
    #  epsY <- 1
    epsY <- model$svEpsY
    alphaY <- model$svAlpha/epsY
    cur.vol <- exp(newX[,2])

    if (ncol(x0) == 3)  {
      alphaZ <- model$svAlphaZ/model$svDeltaZ
      rhoBarZ <- sqrt( 1-model$svRhoZ^2 - model$svRhoYZ^2)
      cur.vol <- exp(newX[,2] + newX[,3])
    }

    for (j in 1:nSteps) {
        w1 <- rnorm(len)
        w2 <- rnorm(len)
        newX[,1] <- newX[,1]*exp( (model$r - cur.vol^2/2)*step.seq[j]  + cur.vol*sqrt(step.seq[j])*w1 )

        # new SV factor
        if (useEuler == T)
           newX[,2] <- newX[,2] + model$svAlpha/epsY*(model$svMean - newX[,2])*step.seq[j] +
                  model$svVol/sqrt(epsY)*(model$svRho*w1 + rhoBar*w2)*sqrt(step.seq[j])
        else
           newX[,2] <- model$svMean + exp(-alphaY*step.seq[j])*(newX[,2]-model$svMean) +
                 model$svVol/sqrt(epsY)*sqrt( 1-exp(-2*alphaY*step.seq[j]))/sqrt(2*alphaY)*(model$svRho*w1 +
                 rhoBar*w2)
        # new asset price
         cur.vol <- exp(newX[,2])

        if (ncol(newX) == 3) {
            w3 <- rnorm(len)
            newX[,3] <- model$svMeanZ + exp(-alphaZ*step.seq[j])*(newX[,3]-model$svMeanZ) +
             model$svVolZ/sqrt(model$svDeltaZ)*sqrt(
             1-exp(-2*alphaZ*step.seq[j]))/sqrt(2*alphaZ)*(model$svRhoZ*w1 + model$svRhoYZ*w2 + rhoBarZ*w3)

            cur.vol <- exp( newX[,2] + newX[,3])

        }




    }
    return (newX)
}

####################################
#' Simulate paths of correlated GBM
#'
#' Simulate correlated multivariate Geometric Brownian motion
#' with a given \code{model$rho} and **identical** \code{model$sigma}'s
#' 
#'
#' @param x0 is the starting values (vector)
#' @param dt is the step size
#' @param model contains all the other parameters.
#' In particular, need \code{model$r, model$rho, model$sigma, model$div}
#' Note that \code{model$sigma} is the **volatility** parameter
#' @export
#' @md
sim.gbm.cor <- function( x0, model, dt=model$dt)
{
    len <- nrow(x0)
    sigm <- diag(model$sigma^2) + model$sigma[1]^2*model$rho*(1- diag(ncol(x0)))

    newX <- x0*exp( rmvnorm(len, sig=sigm)*sqrt(dt) +
            (model$r- model$div- model$sigma[1]^2/2)*dt)

    return (newX)
}

