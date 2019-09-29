##########################################
#' Longstaff Schwartz Algorithm using Bouchard-Warin method
#'
#' Uses the Bouchard-Warin recursive partitioning to created N-d trees
#' for local linear regression fits
#' @param N     is the number of paths
#' @param model must contain the following fields: \code{T, dt, dim, nChildren},
#'        \code{sim.func, x0, r}
#'   @return a list with the following fields:
#'  \code{price} is the scalar optimal reward;
#'  \code{tau} is a vector of stopping times over in-sample paths;
#'  \code{test} is a vector of out-of-sample pathwise rewards;
#'  \code{val} is a vector of in-sample pathwise rewards
#'  \code{timeElapsed} timing information based on \code{Sys.time}
#'
#' @examples
#' option.payoff <-sv.put
#' set.seed(1)
#' modelSV5 <- list(K=100,x0=c(90, log(0.35)),r=0.0225,div=0,sigma=1,
#'    T=50/252,dt=1/252,svAlpha=0.015,svEpsY=1,svVol=3,svRho=-0.03,svMean=2.95,
#'    eulerDt=1/2520, dim=2,sim.func=sim.expOU.sv,nChildren=10)
#' putPr <- osp.probDesign.piecewisebw(20000,modelSV5)
#' putPr$price
#'   \# get [1] 17.30111
#' @export
###########################################
osp.probDesign.piecewisebw <- function(N,model,test.paths=NULL, verb=0)
{
  t.start <- Sys.time()
  M <- model$T/model$dt
  grids <- list()
  preds <- list()
  meanPrice <- rep(0,M)
  bnd <- array(0, dim=c(M,3))  
  # used for 1-D case, saves the stopping boundary for each time step
  bnd[M,] <- c(model$K,model$K,N)
  
  # in 1-d save all the models to analyze the fits
  if (model$dim == 1) {
    all.models <- array(list(NULL), dim=c(M,model$nChildren))
    all.bounds <- array(0, dim=c(M,model$nChildren))
  }

  # Build the grids from a global simulation of X_{1:T}
  grids[[1]] <- model$sim.func( matrix(rep(model$x0, N), nrow=N,byrow=T), model, model$dt)
  for (i in 2:M)
    grids[[i]] <- model$sim.func( grids[[i-1]], model, model$dt)
  
  # make sure to have something for out-of-sample
  if (is.null(test.paths)) {
    test.paths <- list()
    for (i in 1:M)
      test.paths[[i]] <- grids[[i]][1:min(N,100),,drop=F]
  }
  
  # initialize the continuation values/stopping times
  contValue <- exp(-model$r*model$dt)*option.payoff( grids[[M]], model$K)
  tau <- rep(model$T, N)
  meanPrice[M] <- mean(contValue)  # average price along paths (for debugging purpose)
  test.value <- exp(-model$r*model$dt)*option.payoff(test.paths[[M]], model$K)
  
  ###### Main loop: Backward step in time
  # Estimate T(t,x)
  for (i in (M-1):1)
  {
    # forward predict
    if (model$dim == 1) {
      preds[[i]] <- treeDivide.BW.1d( data.frame(cont=contValue,grid=grids[[i]]), 2, 
                                      model, test=test.paths[[i]])
      all.models[i,] <- preds[[i]]$lm.model
      all.bounds[i,] <- preds[[i]]$boundaries
      #preds[[i]]$out.sample <- predict( all.models[i,], test.paths[[i]])
    }
    else
      preds[[i]] <- treeDivide.BW( data.frame(cont=contValue,grid=grids[[i]]), 2, 
                                   model, test=test.paths[[i]])
    
    # compute T(t,x) = C(t,x) - h(x)
    immPayoff <- option.payoff(grids[[i]],model$K)
    timingValue <- preds[[i]]$in.sample - immPayoff
    
    # same on the out-of-sample paths
    test.payoff <- option.payoff(test.paths[[i]],model$K)
    test.tv <- preds[[i]]$out.sample - test.payoff
    
    # figure out the boundary in 1d
    if (model$dim == 1) {
      stop.ndx <- which( timingValue < 0 & grids[[i]][,1] < model$K)
      bnd[i,] <- c(max( grids[[i]][stop.ndx,1]),quantile( grids[[i]][stop.ndx,1], 0.98), length(stop.ndx)/N )
      if (verb > 0 & i > 1)
        if (i %% verb == 1)  {
          plot(grids[[i]][,1] ,timingValue,cex=0.5,col="red",xlim=c(model$K*0.6,model$K*1.6),ylim=c(-0.5,2),
               xlab="S_t", ylab="Timing Value",cex.lab=1.1,main=i)
          abline(h=0, lty=2, lwd=2)
          rug( grids[[i]][1:1000,1], col="blue")
        }
    }
    
    # paths on which stop
    stopNdx <- which( timingValue <= 0 & immPayoff > 0)
    contValue[stopNdx] <- immPayoff[stopNdx]
    contValue <- exp(-model$r*model$dt)*contValue
    tau[stopNdx] <- i*model$dt
    meanPrice[i] <- mean(contValue)
    
    # paths on which stop out-of-sample
    stop.test <- which( test.tv <= 0 & test.payoff > 0)
    test.value[stop.test] <- test.payoff[ stop.test]
    test.value <- exp(-model$r*model$dt)*test.value
  }
  
  # final answer
  price <- mean(contValue)
  print(paste(price, " and out-of-sample ", mean(test.value) ))
  if (verb >0) {
    plot(seq(model$dt,model$T,by=model$dt),bnd[,2], lwd=2, type="l", ylim=c(33,41),
         xlab="Time t", ylab="Stopping Boundary",col="red", cex.lab=1.2)
    points(0,40, cex=2, pch= 18,col='black')
    text(0.5,model$K-2, 'Continue',cex=1.7)
    text(0.75,model$K-6, 'Stop', cex=1.7)
    
  }
  
  # Return a list with the following fields:
  #  price is the scalar optimal reward
  #  tau is a vector of stopping times over in-sample paths
  #  test is a vector of out-of-sample pathwise rewards
  #  val is a vector of in-sample pathwise rewards
  return( list(price=price,tau=tau,test=test.value,val=contValue, timeElapsed=Sys.time()-t.start))
}
