##########################################
#' @title Longstaff Schwartz Algorithm using the Bouchard-Warin method
#'
#' Uses the Bouchard-Warin recursive partitioning to create N-d trees
#' for local linear regression fits. Each tree node contains N/model$nBins^model$dim inputs.
#' @details Calls \link{treeDivide.BW} to create the equi-probable partitions.
#' Must have N/model$nBins^model$dim as an integer.
#' 
#' @param N     the number of forward training paths
#' @param model a list defining all model parameters. Must contain the following fields:
#' \cr \code{T, dt, dim, nBins},
#'        \code{sim.func, x0, r, payoff.func}
#' @param verb if specified, produces plots of the 1-dim fit every \code{verb} time-steps
#' [default is zero, no plotting]
#' @param tst.paths (optional) a list containing out-of-sample paths to obtain a price estimate
#'       
#' @return a list with the following fields:
#' \itemize{
#' \item \code{price} is the scalar optimal reward;
#' \item \code{tau} is a vector of stopping times over in-sample paths;
#' \item \code{test} is a vector of out-of-sample pathwise rewards;
#' \item \code{val} is a vector of in-sample pathwise rewards
#' \item \code{timeElapsed} total running time based on \code{Sys.time}
#' }
#' 
#' @references 
#'  Bruno Bouchard and Xavier Warin. Monte-Carlo valorisation of American options: facts and new
#' algorithms to improve existing methods. In R. Carmona, P. Del Moral, P. Hu, and N. Oudjane, editors,
#' Numerical Methods in Finance, volume 12 of Springer Proceedings in Mathematics. Springer, 2011.
#' 
#' @examples
#' set.seed(1)
#' modelSV5 <- list(K=100,x0=c(90, log(0.35)),r=0.0225,div=0,sigma=1,
#'    T=50/252,dt=1/252,svAlpha=0.015,svEpsY=1,svVol=3,svRho=-0.03,svMean=2.95,
#'    eulerDt=1/2520, dim=2,sim.func=sim.expOU.sv,nBins=10,payoff.func=sv.put.payoff)
#' putPr <- osp.probDesign.piecewisebw(20000,modelSV5)
#' putPr$price
#'   # get [1] 17.30111
#' @export
###########################################
osp.probDesign.piecewisebw <- function(N,model,tst.paths=NULL, verb=0)
{
  t.start <- Sys.time()
  M <- as.integer(round(model$T/model$dt))
  grids <- list()
  preds <- list()
  
  if (is.null(model$nBins) )
    stop("Missing model parameters: must specify nBins (number of partitions per dimension)")
  
  
  # in 1-d save all the models to analyze the fits
  if (model$dim == 1) {
    all.models <- array(list(NULL), dim=c(M,model$nBins))
    all.bounds <- array(0, dim=c(M,model$nBins))
    bnd <- array(0, dim=c(M,3))  
    # used for 1-D case, saves the stopping boundary for each time step
    bnd[M,] <- c(model$K,model$K,N)
  }

  # Build the grids from a global simulation of X_{1:T}
  grids[[1]] <- model$sim.func( matrix(rep(model$x0, N), nrow=N,byrow=T), model, model$dt)
  for (i in 2:M)
    grids[[i]] <- model$sim.func( grids[[i-1]], model, model$dt)
  
  # make sure to have something for out-of-sample
  if (is.null(tst.paths)) {
    test.paths <- list()
    for (i in 1:M)
      test.paths[[i]] <- grids[[i]][1:min(N,100),,drop=F]
  }
  else 
    test.paths <- tst.paths
  
  # initialize the continuation values/stopping times
  contValue <- exp(-model$r*model$dt)*model$payoff.func( grids[[M]], model)
  tau <- rep(model$T, N)
  test.value <- exp(-model$r*model$dt)*model$payoff.func(test.paths[[M]], model)
  
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
    immPayoff <- model$payoff.func(grids[[i]],model)
    timingValue <- preds[[i]]$in.sample - immPayoff
    
    # same on the out-of-sample paths
    test.payoff <- model$payoff.func(test.paths[[i]],model)
    test.tv <- preds[[i]]$out.sample - test.payoff
    
    # figure out the Put boundary in 1d
    if (model$dim == 1) {
      stop.ndx <- which( timingValue < 0 & model$payoff.func(grids[[i]], model) > 0)
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
    
    # paths on which stop out-of-sample
    stop.test <- which( test.tv <= 0 & test.payoff > 0)
    test.value[stop.test] <- test.payoff[ stop.test]
    test.value <- exp(-model$r*model$dt)*test.value
  }
  
  # final answer
  price <- mean(contValue) 
  out.message <- paste("In-sample estimated price: ", round(price, digits=3))
  if (is.null(tst.paths) == FALSE)
    out.message <- paste(out.message, " and out-of-sample ", round(mean(test.value), digits=3 ))
  cat(out.message)
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
  # time elapsed (for benchmarking)
  return( list(price=price,tau=tau,test=test.value,val=contValue, timeElapsed=Sys.time()-t.start))
}
