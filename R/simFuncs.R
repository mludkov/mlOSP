
####################################
#' Simulate paths of Geometric Brownian Motion with constant parameters
#'
#' @details Simulate from \eqn{p(X_t|X_{t-1})}.
#' Use log-normal transition density specified by the \code{model} parameters
#' @param x0 is the starting values (matrix of size N x model$dim)
#' @param dt is the step size in time. Defaults to \code{model$dt}
#' @param model a list containing all the other parameters, including volatility \code{model$sigma},
#' interest rate \code{model$r} and continuous dividend yield \code{model$div}.
#' @return a vector of same dimensions as x0
#' @export
#'
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
#' @details Uses \code{model$alpha} for the mean-reversion strength, 
#' \code{model$meanrev} for the mean-reversion level, and \code{model$sigma} for the volatility.
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
#' @param useEuler flag, whether to use the exact transition for the StochVol factor, or its 
#' Euler approximation
#' @details:  Need the following fields in model: \code{svMean} (mean-reversion level), 
#' \code{svAlpha} (mean-reversion strength), \code{svEpsY} (fast scaling parameter), 
#' \code{svVol} (volatility of volatility), \code{svRho} (correlation with asset S).
#' For 2-factor also need: \code{svMeanZ} (slow scale mean-reversion level), \code{svAlphaZ} 
#' (mean-reversion strength), \code{svDeltaZ} (slow scaling parameter),
#' \code{svVolZ} (Z volatility), \code{svRhoZ} (correlation between Z and S), \code{svRhoYZ} 
#' ( correlation between the fast and slow SV factors)
sim.expOU.sv <- function(x0, model, dt=model$dt,useEuler=FALSE)
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
#' Simulate paths of correlated GBM with a constant correlation
#'
#' @details Simulate correlated multivariate Geometric Brownian motion
#' with a given \code{model$rho} and **identical** \code{model$sigma}'s.
#' Calls \code{rmvnorm} from \pkg{mvtnorm}
#' 
#' 
#'
#' @param x0 is the starting values (vector)
#' @param dt is the step size
#' @param model contains all the other parameters.
#' In particular, need \code{model$r, model$rho, model$sigma, model$div, model$dim}
#' Note that \code{model$sigma} is the **volatility** parameter (scalar)
#' @export
#' @md
sim.gbm.cor <- function( x0, model, dt=model$dt)
{
    len <- nrow(x0)
    sig <- rep(model$sigma[1], model$dim)
    sigm <- diag(sig^2) + kronecker(sig,t(sig))*model$rho*(1- diag(model$dim))

    newX <- x0*exp( rmvnorm(len, sig=sigm)*sqrt(dt) +
            (model$r- model$div- model$sigma[1]^2/2)*dt)

    return (newX)
}

####################################
#' Simulate paths of correlated GBM
#'
#' @details Simulate correlated multivariate Geometric Brownian motion
#' with a given \code{model$rho} and arbitrary \code{model$sigma}'s. Calls \code{rmvnorm}
#' from \pkg{mvtnorm}
#' 
#'
#' @param x0 is the starting values. Should be a matrix of size N x \code{model$dim})
#' @param dt is the step size (Defaults to \code{model$dt})
#' @param model contains all the other parameters.
#' In particular, need \code{model$r, model$rho, model$sigma, model$div, model$dim}
#' Note that \code{model$sigma} is the **volatility vector**
#' @return a vector of the new states (same dimension as x0)
#' @export
#' @md
sim.gbm.matrix <- function( x0, model, dt=model$dt)
{
  len <- nrow(x0)
  
  newX <- x0*exp( mvtnorm::rmvnorm(len, sig=model$sigma)*sqrt(dt) +
                    (model$r- model$div- diag(model$sigma)/2)*dt)
  
  return (newX)
}

####################################
#' Simulate 1D Brownian Motion for Asian Options
#' @details first column is t, second column is S_t
#' third column is A_t (arithmetic average)
#' fourth column is tilde{A}_t (geometric average)
#' @inheritParams sim.gbm
#' @export
sim.gbm.asian <- function( x0, model, dt=model$dt)
{
  len <- nrow(x0)
  newX <- x0
  
  newX[,1] <- x0[,1]*exp( (model$r - 0.5*model$sigma^2)*dt + model$sigma*sqrt(dt)*rnorm(len))
  newX[,2] <- x0[,2] + dt # time
  newX[,3] <- x0[,3]*x0[,2]/newX[,2] + newX[,1]/newX[,2]  # (t-1)/t*A_{t-1} + S_t/t
  newX[,4] <- exp( (x0[,2]/newX[,2]*log(x0[,4])) + log(newX[,1])/newX[,2])
  
  return (newX)
}

####################################
#' Simulate 1D Brownian Motion for Moving Average Asian Options
#' @details first column is S_t, other columns are lagged S_t's
#' the lags are in terms of dt
#' @inheritParams sim.gbm
#' @export
sim.gbm.moving.ave <- function( x0, model, dt=model$dt)
{
  len <- nrow(x0)
  newX <- x0
  
  for (j in 2:ncol(x0))
    newX[,j] <- x0[,j-1]
  
  newX[,1] <- x0[,1]*exp( (model$r - 0.5*model$sigma^2)*dt + model$sigma*sqrt(dt)*rnorm(len))
  
  return (newX)
}

