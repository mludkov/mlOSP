#################
#' Regression Monte Carlo via sequential experimental design. The experimental design is augmented 
#' one input at a time, using an Expected Improvement (EI) acquisition function. This is repeated at
#' each time step. The method is likely to be somewhat slow, but highly efficient in its use of underlying
#' simulations. See Gramacy & Ludkovski (2013), Ludkovski (2018) for details.
#' 
#' EI criteria are based on posterior and/or predictive variance and therefore require the use of a 
#' Gaussian-process based surrogate (currently from \pkg{DiceKriging} or \pkg{hetGP}).
#' 
#' @title Sequential design for optimal stopping
#'
#' @details Implements the EI strategy defined in \code{model$ei.func}. Calls \code{lhs} from library \pkg{tgp}.
#' Empirical losses are computed using \code{cf.el} function. The acquisition function is specified via
#' \code{ei.func} which can be \code{csur} (Default), \code{sur}, \code{smcu}, \code{amcu},
#' \code{tmse} and \code{icu}.
#' 
#' The experimental design is initialized via \code{init.size}/\code{init.grid} parameters and then is grown 
#' one input-at-a-time until it is of size \code{model$seq.design.size}. Thus, there are a total of 
#' seq.design.size-init.size sequential iterations.
#' 
#' @param model a list containing all the model parameters. 
#' 
#' The following model parameters are used:
#' \itemize{
#' \item \code{init.size}: size of starting grid (will be generated via lhs sampling if \code{init.grid} is not given)
#' \item \code{pilot.nsims}: number of pilot simulations to create the search space where new inputs will be added
#' (Default is 5*\code{model$init.size})
#' \item \code{cand.len}: number of candidate new inputs to be proposed. The next input is chosen greedily as
#' the candidate that maximizes the EI criterion. Candidate inputs are selected via \code{tgp::lhs}
#' (Default is 500*\code{model$dim})
#' \item \code{lhs.rect}: specification of the bounding hyper-rectangle where search is conducted
#' (Default: construct based on 0.02/0.98 quantiles of the pilot paths in each dimension)
#' \item \code{update.freq}: how often to re-fit the entire GP surrogate as new inputs are added (Default is 10)
#' \item \code{batch.nrep} (REQUIRED): number of replicates at each unique input
#' \item \code{min.lengthscale}: vector with minimum lengthscales of the surrogate (Default: 1% of lhs.rect dimensions)
#' \item \code{max.lengthscale}: vector with maximum lengthscale of the surrogate (Default: 10x each of lhs.rect dimensions)
#' \item \code{ei.func}: acquisition function (cSUR by Default)
#' }
#' 
#' @param method one of \code{km}, \code{trainkm}, \code{homtp} or \code{hetgp} to select the GP emulator to apply
#' @export
#' @return a list containing:
#' \itemize{
#' \item \code{fit} a list of fitted response surfaces.
#' \item \code{timeElapsed},
#' \item \code{nsims} total number of 1-step sim.func calls
#' \item \code{budget} -- number of sequential iterations per time-step
#' \item \code{empLoss} --matrix of empirical losses (rows for time-steps, columns for iterations)
#' \item \code{theta.fit} -- 3d array of estimated lengthscales (sorted by time-steps,iterations,dimensions-of-x)
#' }
#' 
#' @references 
#' Mike Ludkovski, Kriging Metamodels and Experimental Design for Bermudan Option Pricing
#'  Journal of Computational Finance, 22(1), 37-77, 2018
#' @author Mike Ludkovski
osp.seq.design <- function(model,method="km")
{
  M <- model$T/model$dt
  t.start <- Sys.time()
  cur.sim <- 0
  if (is.null(model$ei.func)) 
     model$ei.func <- "csur"
  if (is.null(model$update.freq)) 
    model$update.freq <- 10
  if (is.null(model$pilot.nsims))
    model$pilot.nsims <- 5*model$init.size
  if (is.null(model$cand.len))
    model$cand.len <- 500*model$dim
  if(is.null(model$ucb.gamma))
    model$ucb.gamma <- 1.96

  fits <- list()   # list of emulator objects at each step
  pilot.paths <- list()
  emp.loss <- array(0, dim=c(M,model$seq.design.size-model$init.size))
  update.kernel.iters <- seq(0,model$seq.design.size,by=model$update.freq)   # when to refit the whole GP

  # set-up a skeleton to understand the distribution of X
  pilot.paths[[1]] <- model$sim.func( matrix(rep(model$x0[1:model$dim], model$pilot.nsims),
                                             nrow=model$pilot.nsims, byrow=T), model, model$dt)
  for (i in 2:(M-1))
    pilot.paths[[i]] <- model$sim.func( pilot.paths[[i-1]], model, model$dt)
  pilot.paths[[1]] <- pilot.paths[[3]]
  budget.used <- rep(0,M-1)
  theta.fit <- array(0, dim=c(M,model$seq.design.size-model$init.size+1,model$dim))

  #------- step back in time
  for (i in (M-1):1)
  {
    all.X <- matrix( rep(0, (model$dim+2)*model$seq.design.size), ncol=model$dim+2)

    # construct the input domain where candidates will be looked for
    if (is.null(model$lhs.rect))
      model$lhs.rect <- 0.02
    if (length(model$lhs.rect) == 1)
    {
      lhs.rect <- matrix( rep(0, 2*model$dim), ncol=2)
      # create a box using empirical quantiles of the pilot.paths cloud
      for (jj in 1:model$dim)
        lhs.rect[jj,] <- quantile( pilot.paths[[i]][,jj], c(model$lhs.rect, 1-model$lhs.rect) )
    }
    else   # already specified
      lhs.rect <- model$lhs.rect
      
        

    # Candidate grid of potential NEW sites to add. Will be ranked using the EI acquisition function
    # only keep in-the-money sites
    ei.cands <- tgp::lhs( model$cand.len, lhs.rect )  # from tgp package
    ei.cands <- ei.cands[ model$payoff.func( ei.cands,model) > 0,,drop=F]
    
    if (is.null(model$min.lengthscale)) 
      model$min.lengthscale <- 0.01*diff(t(lhs.rect))
    if (is.null(model$max.lengthscale)) 
      model$max.lengthscale <- 10*diff(t(lhs.rect)) 

    # initial design
    if (is.null(model$init.grid))
      init.grid <- tgp::lhs( model$init.size, lhs.rect)
    else
      init.grid <- model$init.grid
    if (model$dim > 1) {
       K0 <- dim(init.grid)[1]

       # initial conditions for all the forward paths: replicated design with batch.nrep
       big.grid <- init.grid[ rep(1:K0, model$batch.nrep),]
    }
    else {
      K0 <- length(init.grid)
      big.grid <- as.matrix(init.grid[ rep(1:K0, model$batch.nrep)])
    }

    fsim <- forward.sim.policy( big.grid, M-i,fits[i:M],model, compact=T,offset=0)
    cur.sim <- cur.sim + fsim$nsims

    # payoff at t
    immPayoff <- model$payoff.func( init.grid, model)

    # batched mean and variance
    for (jj in 1:K0) {
      all.X[jj,model$dim+1] <- mean( fsim$payoff[ jj + seq(from=0,len=model$batch.nrep,by=K0)]) - immPayoff[ jj]
      all.X[jj,model$dim+2] <- var( fsim$payoff[ jj + seq(from=0,len=model$batch.nrep,by=K0)])
    }
    all.X[1:K0,1:model$dim] <- init.grid  # use  first dim+1 columns for batched GP regression
    k <- K0

    # create the km object
    if (method == "km")
      fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[1:k,model$dim+1]),
                        noise.var=all.X[1:k,model$dim+2]/model$batch.nrep,
                        control=list(trace=F), lower=model$min.lengthscale, covtype=model$kernel.family, 
                        coef.trend=0, coef.cov=model$km.cov, coef.var=model$km.var)
    if (method == "trainkm") {
      fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[1:k,model$dim+1]),
                      noise.var=all.X[1:k,model$dim+2]/model$batch.nrep, covtype=model$kernel.family,
                      control=list(trace=F), lower=model$min.lengthscale, upper=model$max.lengthscale)
      theta.fit[i,k-model$init.size+1,] <- coef(fits[[i]])$range
    }
    if (method == "homtp") {
      fits[[i]] <- hetGP::mleHomTP(X=init.grid, Z=all.X[1:k,model$dim+1],
                                   lower=model$min.lengthscale,upper=model$max.lengthscale,
                      covtype=model$kernel.family, noiseControl=list(nu_bounds=c(2.01,10),sigma2_bounds = c(5e-2, 4)))
      theta.fit[i,k-model$init.size+1,] <- fits[[i]]$theta
    }
    
    if (method== "hetgp") {
      big.payoff <- model$payoff.func(big.grid,model)
      hetData <- find_reps(big.grid, fsim$payoff-big.payoff)
      fits[[i]] <- hetGP::mleHetGP(X = list(X0=hetData$X0, Z0=hetData$Z0,mult=hetData$mult), Z= hetData$Z,
                                   lower = model$min.lengthscale, upper = model$max.lengthscale,
                                   covtype=model$kernel.family)
      theta.fit[i,k-model$init.size+1,] <- fits[[i]]$theta
    }

    #----- active learning loop:
    add.more.sites <- TRUE; k <- K0 + 1
    # to do it longer for later points to minimize error build-up use: *(1 + i/M)
    while (add.more.sites == TRUE)
    {
      # predict on the candidate grid. Need predictive GP mean, posterior GP variance and similation Stdev
      if (method == "km" | method == "trainkm") {
         pred.cands <- predict(fits[[i]],data.frame(x=ei.cands), type="UK")
         cand.mean <- pred.cands$mean
         cand.sd <- pred.cands$sd
         nug <- sqrt(mean(all.X[1:(k-1),model$dim+2])/model$batch.nrep)
      }
      else {   # hetGP/homTP
         pred.cands <- predict(x=ei.cands, object=fits[[i]])
         cand.mean <- pred.cands$mean
         cand.sd <- sqrt(pred.cands$sd2)
         nug <- sqrt(pred.cands$nugs/model$batch.nrep)
      }
      losses <- cf.el(cand.mean, cand.sd)
      
     

      #---------- use active learning measure to select new sites
      #if (model$ei.func == 'ucb')
      #  al.weights <- qEI.ucb(pred.cands, model$ucb.gamma*sqrt(log(k)))
      if (model$ei.func == 'sur')  
        al.weights <- cf.sur(cand.mean, cand.sd, nugget = nug)
      if (model$ei.func == 'tmse')
        al.weights <- cf.tMSE(cand.mean, cand.sd, seps = model$tmse.eps)
      if (model$ei.func == 'mcu')
        al.weights <- cf.mcu(cand.mean, cand.sd)
      if (model$ei.func == 'smcu')
        al.weights <- cf.smcu(cand.mean, cand.sd, model$ucb.gamma)
      if (model$ei.func == 'amcu') {  # adaptive gamma from paper with Xiong
        adaptive.gamma <- (quantile(cand.mean, 0.75) - quantile(cand.mean, 0.25))/mean(cand.sd)
        al.weights <- cf.smcu(cand.mean, cand.sd, adaptive.gamma)
      }
      if (model$ei.func == 'csur')
        al.weights <- cf.csur(cand.mean, cand.sd,nugget=nug)
      if (model$ei.func == 'icu') {
      # Integrated contour uncertainty with weights based on *Hard-coded* log-normal density
        x.dens2 <- dlnorm( ei.cands[,1], meanlog=log(model$x0[1])+(model$r-model$div - 0.5*model$sigma[1]^2)*i*model$dt,
                           sdlog = model$sigma[1]*sqrt(i*model$dt) )
        jdim <- 2
        while (jdim <= model$dim ) {
          
          x.dens2 <- x.dens2*dlnorm( ei.cands[,jdim], meanlog=log(model$x0[jdim])+(model$r-model$div-0.5*model$sigma[jdim]^2)*i*model$dt,
                              sdlog = model$sigma[jdim]*sqrt(i*model$dt) )
          jdim <- jdim+1
        }
        #plot(ei.cands[,1], ei.cands[,2], cex=cf.mcu(cand.mean, cand.sd)*x.dens2*4000,pch=19)
      
        kxprime <- cov_gen(X1 = fits[[i]]$X0, X2 = ei.cands, theta = fits[[i]]$theta, type = fits[[i]]$covtype)
        al.weights <- apply(ei.cands,1, crit_ICU, model=fits[[i]], thres = 0, Xref=ei.cands, 
                            w = x.dens2, preds = pred.cands, kxprime = kxprime)

      }
      # multiply by weights based on the distribution of pilot paths
      dX_1 <- pilot.paths[[i]]; dX_2 <- ei.cands
      for (dd in 1:model$dim) {
        dX_1[,dd] <- dX_1[,dd]/(lhs.rect[dd,2]-lhs.rect[dd,1])
        dX_2[,dd] <- dX_2[,dd]/(lhs.rect[dd,2]-lhs.rect[dd,1])
      }
      # from package lagp
      ddx <- laGP::distance( dX_1, dX_2) #ddx <- distance( pilot.paths[[i]], ei.cands)
      #x.dens <- apply( exp(-ddx/2/(lhs.rect[1,2]-lhs.rect[1,1])), 2, sum)
      x.dens <- apply( exp(-ddx*dim(ei.cands)[1]*0.01), 2, sum)

      
      if (model$ei.func != 'icu') 
          ei.weights <- x.dens*al.weights
      else 
           ei.weights<- al.weights # those are negative for ICU
      emp.loss[i,k-model$init.size] <- sum(losses*x.dens)/sum(x.dens)
      if (is.null(model$el.thresh) == FALSE)  # stop if below the Emp Loss threshold
        if (emp.loss[i,k-model$init.size] < model$el.thresh)
             add.more.sites <- FALSE


      #add.grid <- ei.cands[sample(dim(ei.cands)[1],model$n.propose,replace=T, prob=ei.weights),,drop=F]
      # select site with highest EI score
      add.grid <- ei.cands[ which.max(ei.weights),,drop=F]
      add.grid <- matrix(rep(add.grid[1,,drop=F], model$batch.nrep), nrow=model$batch.nrep, byrow=T)  #batch

      # compute corresponding y-values
      fsim <- forward.sim.policy( add.grid,M-i,fits[i:M],model,offset=0)
      cur.sim <- cur.sim + fsim$nsims

      immPayoff <- model$payoff.func(add.grid, model)
      add.mean <- mean(fsim$payoff - immPayoff)
      add.var <- var(fsim$payoff - immPayoff)
      all.X[k,] <- c(add.grid[1,],add.mean,add.var)

      if (k %in% update.kernel.iters) {
        if (method == "km")
           fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=all.X[1:k,1:model$dim]), response=data.frame(y=all.X[1:k,model$dim+1]),
                          noise.var=all.X[1:k,model$dim+2]/model$batch.nrep, coef.trend=0, coef.cov=model$km.cov,
                          coef.var=model$km.var, control=list(trace=F), 
                          lower=model$min.lengthscale,upper=model$max.lengthscale)
        if (method == "trainkm") {
          fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=all.X[1:k,1:model$dim]), response=data.frame(y=all.X[1:k,model$dim+1]),
                          noise.var=all.X[1:k,model$dim+2]/model$batch.nrep, control=list(trace=F), 
                          lower=model$min.lengthscale, upper=model$max.lengthscale)
          theta.fit[i,k-model$init.size+1,] <- coef(fits[[i]])$range
        }
        if (method == "hetgp") {
           fits[[i]] <- update(object=fits[[i]], Xnew=add.grid, Znew=fsim$payoff-immPayoff,method="mixed",
                               lower = model$min.lengthscale, upper=model$max.lengthscale)
           theta.fit[i,k-model$init.size+1,] <- fits[[i]]$theta
        }
        if (method == "homtp") {
          fits[[i]] <- hetGP::mleHomTP(X=all.X[1:k,1:model$dim], Z=all.X[1:k,model$dim+1],noiseControl=list(nu_bounds=c(2.01,10),sigma2_bounds=c(5e-2,4)),
                                 lower=model$min.lengthscale,upper=model$max.lengthscale, covtype=model$kernel.family)
          theta.fit[i,k-model$init.size+1,] <- fits[[i]]$theta
        }
      }
      else {
        if (method == "km" | method == "trainkm")
           fits[[i]] <- update(fits[[i]],newX=add.grid[1,,drop=F],newy=add.mean,
                               newnoise=add.var/model$batch.nrep,  cov.re=F)
        if (method == "hetgp" | method == "homtp")
           fits[[i]] <- update(object=fits[[i]], Xnew=add.grid, Znew=fsim$payoff-immPayoff, maxit = 0)
      }

      # resample the candidate set
      ei.cands <- tgp::lhs( model$cand.len, lhs.rect )
      ei.cands <- ei.cands[ model$payoff.func( ei.cands,model) > 0,,drop=F]
      if (k >= model$seq.design.size)
         add.more.sites <- FALSE
      k <- k+1
    }
    budget.used[i] <- k
  }

  return (list(fit=fits,timeElapsed=Sys.time()-t.start,nsims=cur.sim,empLoss=emp.loss,
               budget=budget.used,theta=theta.fit))
}


########################################
#### Adaptive batching
osp.batch.design <- function(model,input.domain=NULL, method ="km",inTheMoney.thresh = 0)
{
  M <- model$T/model$dt
  t.start <- Sys.time()

  fits <- list()   # list of fits at each step
  grids <- list()

  cur.sim <- 0

  # set-up pilot design using a forward simulation of X
  if (model$pilot.nsims > 0) {
    grids[[1]] <- model$sim.func( matrix(rep(model$x0[1:model$dim], model$pilot.nsims),
                                         nrow=model$pilot.nsims, byrow=T), model, model$dt)
    for (i in 2:(M-1))
      grids[[i]] <- model$sim.func( grids[[i-1]], model, model$dt)
    grids[[1]] <- grids[[3]]
    cur.sim <- model$pilot.nsims
  }

  ############ step back in time
  design.size <- rep(0,M)

  for (i in (M-1):1) {
    # figure out design size -- this is the number of unique sites, before restricting to in-the-money
    if (length(model$N) == 1)
      design.size[i] <- model$N
    else
      design.size[i] <- model$N[i]

    #figure out batch size
    if (is.null(model$batch.nrep))
      n.reps <- 1
    else if (length(model$batch.nrep) == 1)
      n.reps <- model$batch.nrep
    else
      n.reps <- model$batch.nrep[i]

    if (is.null(input.domain))  {   # empirical design using simulated pilot paths
      init.grid <- grids[[i]]

      init.grid <- init.grid[sample(1:min(design.size[i],dim(init.grid)[1]), design.size[i], replace=F),,drop=F]
    }
    else if (length(input.domain)==2*model$dim | length(input.domain)==1) {
      # space-filling design over a rectangle
      if (length(input.domain) == 1){
        # specifies the box as empirical quantiles, should be 0.01. If zero then use the full range
        my.domain <- matrix( rep(0, 2*model$dim), ncol=2)
        if (input.domain > 0) {
          for (jj in 1:model$dim)
            my.domain[jj,] <- quantile( grids[[i]][,jj], c(input.domain, 1-input.domain))
        }
        else {
          for (jj in 1:model$dim)
            my.domain[jj,] <- range( grids[[i]][,jj] )
        }
      }
      else my.domain <- input.domain  #  user-specified box

      # now choose how to space-fill
      if (is.null(model$qmc.method)) {
        init.grid <- tgp::lhs( design.size[i], my.domain)
      }
      else {
        init.grid <- model$qmc.method( design.size[i], dim=model$dim)
        # rescale to the correct rectangle
        for (jj in 1:model$dim)
          init.grid[,jj] <- my.domain[jj,1] + (my.domain[jj,2]-my.domain[jj,1])*init.grid[,jj]

      }
    }
    else  {   # fixed pre-specified design
      init.grid <- matrix(input.domain,nrow=length(input.domain)/model$dim)
      design.size[i] <- nrow(init.grid)
    }

    init.grid <- init.grid[ model$payoff.func(init.grid, model) > inTheMoney.thresh,,drop=F]

    design.size[i] <- dim(init.grid)[1]
    all.X <- matrix( rep(0, (model$dim+2)*design.size[i]), ncol=model$dim+2)

    # construct replicated design
    if (i == (M-1))
      my.reps <- rep(n.reps,design.size[i])
    else {
      pred <- predict(fits[[i+1]], x=init.grid)
      my.reps <- n.reps*model$batch.func( pred$mean)
      my.reps <- ceiling( model$N*my.reps/sum(my.reps))
      print(sum(my.reps))
    }

    big.grid <- matrix(nrow=sum(my.reps), ncol=model$dim)
    for (jj in 1:model$dim)
      big.grid[,jj] <- rep(init.grid[,jj], my.reps)
    #big.grid <- matrix(rep(t(init.grid), my.reps), ncol = ncol(init.grid), byrow = TRUE)

    fsim <- forward.sim.policy( big.grid, M-i,fits[i:M],model,offset=0)
    immPayoff <- model$payoff.func( init.grid, model)
    cur.sim <- cur.sim + fsim$nsims

    # pre-averaged mean/variance
    #for (jj in 1:design.size[i]) {
    #  all.X[jj,model$dim+1] <- mean( fsim$payoff[ jj + seq(from=0,len=my.reps[jj],by=design.size[i])]) - immPayoff[ jj]
    #  all.X[jj,model$dim+2] <- var( fsim$payoff[ jj + seq(from=0,len=n.reps,by=design.size[i])])
    #}

    all.X[,1:model$dim] <- init.grid  # use the first dim+1 columns for the batched GP regression.

    # create the km object
    if (n.reps > 10 & method == "km")
      fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[,model$dim+1]),
                                   noise.var=all.X[,model$dim+2]/n.reps, control=list(trace=F),
                                   lower=model$min.lengthscale, upper=model$max.lengthscale,
                                   coef.trend=0,coef.cov=model$km.cov, coef.var=model$km.var, covtype=model$kernel.family)
    else if (method == "km")  # manually estimate the nugget for small batches
      fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[,model$dim+1]),
                                   control=list(trace=F), lower=model$min.lengthscale, upper=model$max.lengthscale,
                                   coef.trend=0, coef.cov=model$km.cov, coef.var=model$km.var,
                                   nugget.estim=TRUE, nugget=sqrt(mean(all.X[,model$dim+2])), covtype=model$kernel.family)
    else if (method =="trainkm")
      fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[,model$dim+1]),
                                   control=list(trace=F), lower=rep(0.1,model$dim), upper=model$max.lengthscale,
                                   noise.var=all.X[,model$dim+2]/n.reps, covtype=model$kernel.family)

    else if (n.reps < 10 & method == "lagp")  # laGP library implementation
      fits[[i]]  <- laGP::newGP(X=init.grid, Z=all.X[,model$dim+1],
                                d=list(mle=FALSE, start=model$km.cov), g=list(start=1, mle=TRUE))
    else if(method =="hetgp") {
      big.payoff <- model$payoff.func(big.grid,model)
      hetData <- find_reps(big.grid, fsim$payoff-big.payoff)
      fits[[i]] <- hetGP::mleHetGP(X = list(X0=hetData$X0, Z0=hetData$Z0,mult=hetData$mult), Z= hetData$Z,
                                   lower = model$min.lengthscale, upper = model$max.lengthscale, covtype=model$kernel.family)
      #ehtPred <- predict(x=check.x, object=hetModel)
    }
    else if (method =="homgp") {
      big.payoff <- model$payoff.func(big.grid,model)
      hetData <- hetGP::find_reps(big.grid, fsim$payoff-big.payoff)
      fits[[i]] <- hetGP::mleHomGP(X = list(X0=hetData$X0, Z0=hetData$Z0,mult=hetData$mult), Z= hetData$Z,
                                   lower = model$min.lengthscale, upper = model$max.lengthscale, covtype=model$kernel.family)
    }
    else if (model$dim == 1 & method=="spline")  # only possible in 1D
      fits[[i]] <- smooth.spline(x=init.grid,y=all.X[,2],nknots=model$nk)
    else if (method == "rvm") {
      if (is.null(model$rvm.kernel))
        rvmk <- "rbfdot"
      else
        rvmk <- model$rvm.kernel
      fits[[i]] <- rvm(x=init.grid, y=all.X[,model$dim+1],kernel=rvmk)
    }
    else if (method == "npreg") {
      if (is.null(model$np.kertype))
        kertype = "gaussian"
      else
        kertype <- model$np.kertype
      if (is.null(model$np.regtype))
        regtype <- "lc"
      else
        regtype <- model$np.regtype
      if (is.null(model$np.kerorder))
        kerorder <- 2
      else
        kerorder <- model$np.kerorder

      fits[[i]] <- npreg(txdat = init.grid, tydat = all.X[,model$dim+1],
                         regtype=regtype, ckertype=kertype,ckerorder=kerorder)
    }

  }  # end of loop over time-steps

  return (list(fit=fits,timeElapsed=Sys.time()-t.start,nsims=cur.sim))
}
