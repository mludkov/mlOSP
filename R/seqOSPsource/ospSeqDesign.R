#################
#' @title Sequential design for optimal stopping
#'
#' @details Implements the EI strategy defined in model/al.heuristic. Calls lhs (library \pkg{tgp}).
#' Empirical losses are computed using \code{cf.el} function
#' @param method: either \code{km} or \code{hetgp} to select the GP emulator to apply
#' @export
#' @return a list containing:
#' \code{price}: v(0,x_0); \code{fit} a list of fitted response surfaces. 
#' \code{timeElapsed}, \code{nsims} total number of 1-step sim.func calls
#' \code{empLoss} --vector of empirical losses 
osp.seq.design <- function(model,method="km")
{
  M <- model$T/model$dt
  t.start <- Sys.time()
  cur.sim <- 0
  
  fits <- list()   # list of emulator objects at each step
  pilot.paths <- list()
  emp.loss <- array(0, dim=c(M,model$adaptive.grid.loop-model$init.size))
  update.kernel.iters <- seq(0,model$adaptive.grid.loop,by=20)   # when to refit the whole GP
  
  # set-up a skeleton to understand the distribution of X
  pilot.paths[[1]] <- model$sim.func( matrix(rep(model$x0[1:model$dim], 5*model$init.size), 
                                             nrow=5*model$init.size, byrow=T), model, model$dt)
  for (i in 2:(M-1))
    pilot.paths[[i]] <- model$sim.func( pilot.paths[[i-1]], model, model$dt)
  pilot.paths[[1]] <- pilot.paths[[3]]
  init.grid <- pilot.paths[[M-1]]
  
  ############ step back in time
  for (i in (M-1):1)
  {
    all.X <- matrix( rep(0, (model$dim+2)*model$adaptive.grid.loop), ncol=model$dim+2)
    
    # construct the input domain where candidates will be looked for
    if (is.null(model$lhs.rect))
      model$lhs.rect <- 0.02
    if (length(model$lhs.rect) == 1)
    {
      lhs.rect <- matrix( rep(0, 2*model$dim), ncol=2)
      # create a box using empirical quantiles of the init.grid cloud
      for (jj in 1:model$dim)
        lhs.rect[jj,] <- quantile( init.grid[,jj], c(model$lhs.rect, 1-model$lhs.rect) )
    }
    else   # already specified
      lhs.rect <- model$lhs.rect
    
    # Candidate grid of potential NEW sites to add. Will be ranked using the AL acquisition function
    # only keep in-the-money sites
    ei.cands <- lhs( model$cand.len, lhs.rect )  # from tgp package
    ei.cands <- ei.cands[ option.payoff( ei.cands,model$K) > 0,,drop=F]
    
    # initial design
    if (is.null(model$init.grid))
      init.grid <- lhs( model$init.size, lhs.rect)
    else
      init.grid <- model$init.grid
    K0 <- dim(init.grid)[1]
    
    # initial conditions for all the forward paths: replicated design with km.batch
    big.grid <- init.grid[ rep(1:K0, model$km.batch),]
    
    fsim <- forward.sim.policy( big.grid, M-i,fits[i:M],model, compact=T,offset=0)
    cur.sim <- cur.sim + fsim$nsims
    
    # payoff at t
    immPayoff <- option.payoff( init.grid, model$K)
    
    # batched mean and variance
    for (jj in 1:K0) {
      all.X[jj,model$dim+1] <- mean( fsim$payoff[ jj + seq(from=0,len=model$km.batch,by=K0)]) - immPayoff[ jj]
      all.X[jj,model$dim+2] <- var( fsim$payoff[ jj + seq(from=0,len=model$km.batch,by=K0)])
    }
    all.X[1:K0,1:model$dim] <- init.grid  # use  first dim+1 columns for batched GP regression
    k <- K0
    
    # create the km object
    if (method == "km")
      fits[[i]] <- km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[1:k,model$dim+1]),
                        noise.var=all.X[1:k,model$dim+2]/model$km.batch,
                        control=list(trace=F), lower=rep(0.1,model$dim), covtype=model$covfamily,
                        coef.trend=0, coef.cov=model$km.cov, coef.var=model$km.var)
    else {
      big.payoff <- option.payoff(big.grid,model$K)
      hetData <- find_reps(big.grid, fsim$payoff-big.payoff)
      fits[[i]] <- hetGP::mleHetGP(X = list(X0=hetData$X0, Z0=hetData$Z0,mult=hetData$mult), Z= hetData$Z,
                                   lower = rep(0.1,model$dim), upper = model$km.upper,
                                   covtype=model$covfamily)
    } 
    
    # active learning loop:
    # to do it longer for later points to minimize error build-up use: *(1 + i/M)
    for (k in (K0+1):round(model$adaptive.grid.loop))
    {
      # predict on the candidate grid. Need predictive GP mean, posterior GP variance and similation Stdev
      if (method == "km") {
         pred.cands <- predict(fits[[i]],data.frame(x=ei.cands), type="UK")
         cand.mean <- pred.cands$mean
         cand.sd <- pred.cands$sd
         nug <- sqrt(mean(all.X[1:(k-1),model$dim+2])/model$km.batch)
      }
      else {   # hetGP
         pred.cands <- predict(x=ei.cands, object=fits[[i]])
         cand.mean <- pred.cands$mean
         cand.sd <- sqrt(pred.cands$sd2)
         nug <- sqrt(pred.cands$nugs/model$km.batch)
      }
      losses <- cf.el(cand.mean, cand.sd)
      
      # use active learning measure to select new sites
      #if (model$al.heuristic == 'ucb')
      #  al.weights <- qEI.ucb(pred.cands, model$ucb.gamma*sqrt(log(k)))
      if (model$al.heuristic == 'sur')  # the only one tested right now
        al.weights <- cf.sur(cand.mean, cand.sd, nugget = nug)
      
      # multiply by weights based on the distribution of pilot paths
      dX_1 <- pilot.paths[[i]]; dX_2 <- ei.cands
      for (dd in 1:model$dim) {
        dX_1[,dd] <- dX_1[,dd]/(lhs.rect[dd,2]-lhs.rect[dd,1])
        dX_2[,dd] <- dX_2[,dd]/(lhs.rect[dd,2]-lhs.rect[dd,1])
      }
      ddx <- distance( dX_1, dX_2) #ddx <- distance( pilot.paths[[i]], ei.cands)
      #x.dens <- apply( exp(-ddx/2/(lhs.rect[1,2]-lhs.rect[1,1])), 2, sum)
      x.dens <- apply( exp(-ddx*10), 2, sum)
      
      #x.dens <- dnorm( (log( ei.cands/model$x0[1]) - (model$r - 0.5*model$sigma^2)*i*model$dt)/model$sigma/sqrt(i*model$dt))
      
      ei.weights <- x.dens*al.weights
      emp.loss[i,k-model$init.size] <- sum(losses*x.dens)/sum(x.dens)
      
      #add.grid <- ei.cands[sample(dim(ei.cands)[1],model$n.propose,replace=T, prob=ei.weights),,drop=F]
      # select site with highest EI score
      add.grid <- ei.cands[ which.max(ei.weights),,drop=F]
      add.grid <- matrix(rep(add.grid[1,,drop=F], model$km.batch), nrow=model$km.batch, byrow=T)  #batch 
      
      # compute corresponding y-values
      fsim <- forward.sim.policy( add.grid,M-i,fits[i:M],model,offset=0)
      cur.sim <- cur.sim + fsim$nsims
      
      immPayoff <- option.payoff(add.grid, model$K)
      add.mean <- mean(fsim$payoff - immPayoff)
      add.var <- var(fsim$payoff - immPayoff)
      all.X[k,] <- c(add.grid[1,],add.mean,add.var)
      
      if (k %in% update.kernel.iters) {
        if (method == "km")         
           fits[[i]] <- km(y~0, design=data.frame(x=all.X[1:k,1:model$dim]), response=data.frame(y=all.X[1:k,model$dim+1]),
                          noise.var=all.X[1:k,model$dim+2]/model$km.batch, coef.trend=0, coef.cov=model$km.cov, 
                          coef.var=model$km.var, control=list(trace=F), lower=rep(0.1, model$dim),upper=model$km.upper)
        if (method == "hetgp")
           fits[[i]] <- update(object=fits[[i]], Xnew=add.grid,Znew=fsim$payoff-immPayoff,method="mixed",
                               lower = rep(0.1,model$dim), upper=model$km.upper)
      }
      else {
        if (method == "km")
           fits[[i]] <- update(fits[[i]],newX=add.grid[1,,drop=F],newy=add.mean,
                               newnoise=add.var/model$km.batch,  cov.re=F)
        if (method == "hetgp")
           fits[[i]] <- update(object=fits[[i]], Xnew=add.grid, Znew=fsim$payoff-immPayoff, maxit = 0)
      }
      
      # resample the candidate set
      ei.cands <- lhs( model$cand.len, lhs.rect )
      ei.cands <- ei.cands[ option.payoff( ei.cands,model$K) > 0,,drop=F]
    }
  }
  
  return (list(fit=fits,timeElapsed=Sys.time()-t.start,nsims=cur.sim,empLoss=emp.loss))
}