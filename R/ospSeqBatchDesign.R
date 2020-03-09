#################
#' @title Adaptive Batch design for optimal stopping
#'
#' @details Implements the EI strategy defined in model/al.heuristic. Calls lhs (library \pkg{tgp}).
#' Empirical losses are computed using \code{cf.el} function
#' @param method: either \code{km} or \code{hetgp} to select the GP emulator to apply
#' @export
#' @return a list containing:
#' \code{price}: v(0,x_0); \code{fit} a list of fitted response surfaces. 
#' \code{timeElapsed}, \code{nsims} total number of 1-step sim.func calls
#' \code{empLoss} --vector of empirical losses 
#' \code{ndesigns}: number of different designs
#' \code{batches}: batch size
osp.seq.batch.design <- function(model, method="km", nt = 1000, t0 = 0.01)
{
  M <- model$T/model$dt
  t.start <- Sys.time()
  cur.sim <- 0
  if (is.null(model$ei.func)) {
    model$ei.func <- "csur"
  }
  if (is.null(model$update.freq)) {
    model$update.freq <- 10
  }
  if (is.null(model$batch.heuristic)) {
    model$batch.heuristic <- 'fb'
  }
  
  # parameters in absur
  if (is.null(model$total.budget)) {
    model$total.budget = model$seq.design.size * model$batch.nrep
  }
  if (is.null(model$r.cand)) {
    model$r.cand = c(model$batch.nrep, model$batch.nrep)
  }
  r_lower = model$r.cand[1]
  r_upper = min(model$r.cand[length(model$r.cand)], 0.1 * model$total.budget)
  r_interval = seq(r_lower, r_upper, length = 1000)
  theta_for_optim = c(0.1371, 0.000815, 1.9871E-6)
  batch_matrix <- matrix(rep(0, M*model$seq.design.size), ncol=M)
  
  # parameter in adsa and ddsa
  c_batch = 20 / model$dim
  
  fits <- list()   # list of emulator objects at each step
  pilot.paths <- list()
  emp.loss <- array(0, dim=c(M,model$seq.design.size-model$init.size))
  update.kernel.iters <- seq(0,model$seq.design.size,by=model$update.freq)   # when to refit the whole GP
  
  # set-up a skeleton to understand the distribution of X
  pilot.paths[[1]] <- model$sim.func( matrix(rep(model$x0[1:model$dim], 5*model$init.size),
                                             nrow=5*model$init.size, byrow=T), model, model$dt)
  for (i in 2:(M-1)) {
    pilot.paths[[i]] <- model$sim.func( pilot.paths[[i-1]], model, model$dt)
  }
  pilot.paths[[1]] <- pilot.paths[[3]]
  init.grid <- pilot.paths[[M-1]]
  budget.used <- rep(0,M-1)
  theta.fit <- array(0, dim=c(M,model$seq.design.size-model$init.size+1,model$dim))
  
  ############ step back in time
  for (i in (M-1):1)
  {
    all.X <- matrix( rep(0, (model$dim+2)*model$seq.design.size), ncol=model$dim+2)
    
    # construct the input domain where candidates will be looked for
    if (is.null(model$lhs.rect)) {
      model$lhs.rect <- 0.02
    }
    if (length(model$lhs.rect) == 1) {
      lhs.rect <- matrix( rep(0, 2*model$dim), ncol=2)
      # create a box using empirical quantiles of the init.grid cloud
      for (jj in 1:model$dim)
        lhs.rect[jj,] <- quantile( init.grid[,jj], c(model$lhs.rect, 1-model$lhs.rect) )
    } else { # already specified 
      lhs.rect <- model$lhs.rect
    }
    
    if (is.null(model$min.lengthscale)) {
      model$min.lengthscale <- rep(0.1, model$dim) 
    }
    
    # Candidate grid of potential NEW sites to add. Will be ranked using the EI acquisition function
    # only keep in-the-money sites
    ei.cands <- lhs( model$cand.len, lhs.rect )  # from tgp package
    ei.cands <- ei.cands[ model$payoff.func( ei.cands,model) > 0,,drop=F]
    
    
    # initial design
    if (is.null(model$init.grid)) {
      init.grid <- lhs( model$init.size, lhs.rect)
    } else {
      init.grid <- model$init.grid
    }
    if (model$dim > 1) {
      K0 <- dim(init.grid)[1]
      
      # initial conditions for all the forward paths: replicated design with batch.nrep
      big.grid <- init.grid[ rep(1:K0, model$batch.nrep),]
    } else {
      K0 <- length(init.grid)
      big.grid <- as.matrix(init.grid[ rep(1:K0, model$batch.nrep)])
    }
    
    fsim <- forward.sim.policy( big.grid, M-i,fits[i:M],model, compact=T, offset=0)
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
    batch_matrix[1:K0, i] <- model$batch.nrep
    
    # create the km object
    if (method == "km") {
      fits[[i]] <- km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[1:k,model$dim+1]),
                      noise.var=all.X[1:k,model$dim+2]/model$batch.nrep,
                      control=list(trace=F), lower=model$min.lengthscale, covtype=model$kernel.family, 
                      nugget.estim= TRUE,
                      coef.trend=0, coef.cov=model$km.cov, coef.var=model$km.var)
    }
    if (method == "trainkm") {
      fits[[i]] <- km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[1:k,model$dim+1]),
                      noise.var=all.X[1:k,model$dim+2]/model$batch.nrep, covtype=model$kernel.family,
                      control=list(trace=F), lower=model$min.lengthscale, upper=model$max.lengthscale)
      theta.fit[i,k-model$init.size+1,] <- coef(fits[[i]])$range
    }
    if (method == "homtp") {
      fits[[i]] <- hetGP::mleHomTP(X=all.X[1:k,1:model$dim], Z=all.X[1:k,model$dim+1],
                                   lower=model$min.lengthscale,upper=model$max.lengthscale,
                                   covtype=model$kernel.family, noiseControl=list(nu_bounds=c(2.01,5)))
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
    
    # initialize gamma for RB and MLB
    gamma <- sqrt(mean(all.X[1:K0,model$dim+2])/model$batch.nrep) / 2
    r_batch = model$batch.nrep
    
    # active learning loop
    add.more.sites <- TRUE
    k <- K0 + 1
    running_budget = model$batch.nrep*K0
    n <- K0 + 1
    
    # active learning loop:
    # to do it longer for later points to minimize error build-up use: *(1 + i/M)
    while(add.more.sites)
    {
      # predict on the candidate grid. Need predictive GP mean, posterior GP variance and similation Stdev
      if (method == "km" | method == "trainkm") {
        pred.cands <- predict(fits[[i]],data.frame(x=ei.cands), type="UK")
        cand.mean <- pred.cands$mean
        cand.sd <- pred.cands$sd
        nug <- sqrt(mean(all.X[1:(k-1), model$dim+2]))
        nug <- rep(nug, length(cand.mean))
      } else {   # hetGP/homTP
        pred.cands <- predict(x=ei.cands, object=fits[[i]])
        cand.mean <- pred.cands$mean
        cand.sd <- sqrt(pred.cands$sd2)
        nug <- sqrt(pred.cands$nugs)
      }
      losses <- cf.el(cand.mean, cand.sd)

      # multiply by weights based on the distribution of pilot paths
      dX_1 <- pilot.paths[[i]]; dX_2 <- ei.cands
      for (dd in 1:model$dim) {
        dX_1[,dd] <- dX_1[,dd]/(lhs.rect[dd,2]-lhs.rect[dd,1])
        dX_2[,dd] <- dX_2[,dd]/(lhs.rect[dd,2]-lhs.rect[dd,1])
      }
      # from package lagp
      ddx <- distance( dX_1, dX_2) #ddx <- distance( pilot.paths[[i]], ei.cands)
      x.dens <- apply( exp(-ddx*dim(ei.cands)[1]*0.01), 2, sum)
      emp.loss[i,k-model$init.size] <- sum(losses*x.dens)/sum(x.dens)
      if (is.null(model$el.thresh) == FALSE) { # stop if below the Emp Loss threshold
        if (emp.loss[i,k-model$init.size] < model$el.thresh) {
          add.more.sites <- FALSE
        }
      }
      
      # use active learning measure to select new sites and associated batch size
      if (model$ei.func == 'absur') {
          overhead = 3 * model$dim * CalcOverhead(theta_for_optim[1], theta_for_optim[2], theta_for_optim[3], k + 1)
          al.weights <- cf.absur(cand.mean, cand.sd, nug, r_interval, overhead, t0)
          
          # select site with highest EI score
          x.dens.matrix <- matrix(x.dens, nrow=length(x.dens), ncol=length(r_interval))
          ei.weights <- x.dens.matrix * al.weights
          
          # select next input location
          max_index <- which.max(ei.weights)
          x_index <- max_index %% length(x.dens)
          x_index <- ifelse(x_index, x_index, length(x.dens))
          add.grid <- ei.cands[x_index,,drop=F]
          
          # select associated batch size
          r_index <- (max_index - 1) / model$cand.len + 1
          r_batch <- min(round(r_interval[r_index]), model$total.budget - running_budget)
      } else {
        #---------- use active learning measure to select new sites
        #if (model$ei.func == 'ucb')
        #  al.weights <- qEI.ucb(pred.cands, model$ucb.gamma*sqrt(log(k)))
        if (model$ei.func == 'sur')  
          al.weights <- cf.sur(cand.mean, cand.sd, nugget = nug / sqrt(model$batch.nrep))
        if (model$ei.func == 'tmse')
          al.weights <- cf.tMSE(cand.mean, cand.sd, seps = model$tmse.eps)
        if (model$ei.func == 'mcu')
          al.weights <- cf.mcu(cand.mean, cand.sd)
        if (model$ei.func == 'smcu')
          al.weights <- cf.smcu(cand.mean, cand.sd, model$ucb.gamma)
        if (model$ei.func == 'amcu') {  # adaptive gamma from paper with Xiong
          adaptive.gamma <- (quantile(cand.mean, 0.75, na.rm = T) - quantile(cand.mean, 0.25, na.rm = T))/mean(cand.sd)
          al.weights <- cf.smcu(cand.mean, cand.sd, adaptive.gamma)
        }
        if (model$ei.func == 'csur')
          al.weights <- cf.csur(cand.mean, cand.sd,nugget=nug / sqrt(model$batch.nrep))
        if (model$ei.func == 'icu' || model$batch.heuristic == 'adsa' || model$batch.heuristic == 'ddsa') {
          # Integrated contour uncertainty with weights based on *Hard-coded* log-normal density
          if (model$dim >= 2) {
            x.dens2 <- dlnorm( ei.cands[,1], meanlog=log(model$x0[1])+(model$r-model$div - 0.5*model$sigma[1]^2)*i*model$dt,
                               sdlog = model$sigma[1]*sqrt(i*model$dt) )
            x.dens2 <- x.dens2*dlnorm( ei.cands[,2], meanlog=log(model$x0[2])+(model$r-model$div-0.5*model$sigma[2]^2)*i*model$dt,
                                       sdlog = model$sigma[2]*sqrt(i*model$dt) )
          }
          if (model$dim == 3) {
            x.dens2 <- x.dens2*dlnorm( ei.cands[,3], meanlog=log(model$x0[3])+(model$r-model$div-0.5*model$sigma[3]^2)*i*model$dt,
                                       sdlog = model$sigma[3]*sqrt(i*model$dt) )
          }
          #plot(ei.cands[,1], ei.cands[,2], cex=cf.mcu(cand.mean, cand.sd)*x.dens2*4000,pch=19)
          
          if (model$ei.func == 'icu') {
            kxprime <- cov_gen(X1 = fits[[i]]$X0, X2 = ei.cands, theta = fits[[i]]$theta, type = fits[[i]]$covtype)
            al.weights <- apply(ei.cands,1, crit_ICU, model=fits[[i]], thres = 0, Xref=ei.cands, 
                                w = x.dens2, preds = pred.cands, kxprime = kxprime)
          }
          
        }
      
        # select site with highest EI score
        if (model$ei.func != 'icu') {
          ei.weights <- x.dens*al.weights
        } else {
          ei.weights<- al.weights # those are negative for ICU
        }
        
        x_index <- which.max(ei.weights)
        add.grid <- ei.cands[x_index,,drop=F]
      }
      
      # use active batching algorithms to select batch size
      if (model$batch.heuristic == 'fb') {
        r_batch = model$batch.nrep
      } else {
        r0 = min(model$total.budget - running_budget, round(c_batch*sqrt(k)))
        if (model$batch.heuristic == 'rb') {
          rb_batch <- batch.rb(cand.sd[x_index], model$r.cand, r_batch, nug[x_index], gamma)
          r_batch <- min(rb_batch$roptim, model$total.budget - running_budget)
          gamma <- rb_batch$gamma
        } 
        if (model$batch.heuristic == 'mlb') {
          mlb_batch <- batch.mlb(cand.sd[x_index], model$r.cand, nug[x_index], gamma)
          r_batch <- min(mlb_batch$roptim, model$total.budget - running_budget)
          gamma <- mlb_batch$gamma
        }
        if (model$batch.heuristic == 'adsa') {
          adsa_batch <- batch.adsa(fits[[i]], batch_matrix[1:(n - 1), i], ei.cands, x.dens2, add.grid, r0, nug, method)
          add.grid <- adsa_batch$x_optim
          r_batch <- adsa_batch$r_optim
        }
        if (model$batch.heuristic == 'ddsa') {
          if (k%%2) {
            r_batch <- batch.ddsa(fits[[i]], batch_matrix[1:(n - 1), i], ei.cands, x.dens2, r0, method)$r_new
          } else {
            r_batch <- r0
          }
        }
      }
      
      # if using up all budget, move on to next time step
      if (is.numeric(r_batch) && r_batch == 0) {
        add.more.sites <- FALSE
        next
      }
      
      if(is.null(add.grid) || model$batch.heuristic == 'ddsa' && k%%2) { # Reallocation
        
        # Indexes for inputs which receives further allocation
        idx_diff = which(r_batch != batch_matrix[1:(n - 1), i])
        r_seq_diff = r_batch[idx_diff] - batch_matrix[idx_diff, i]
        add.mean = rep(0, length(idx_diff))
        add.var = rep(0, length(idx_diff))
        
        if (method == 'hetgp') {
          newX <- matrix(, nrow = 0, ncol = model$dim)
          newY <- vector()
        }
        
        # Reallocates new observations
        for(l in 1:length(idx_diff)) {
          idx = idx_diff[l]
          change.grid = matrix(all.X[idx,1:model$dim], nrow = r_seq_diff[l], ncol = model$dim, byrow = T)
          
          # compute corresponding y-values
          fsim <- forward.sim.policy( change.grid, M-i, fits[i:M], model, compact=T, offset=0)
          cur.sim <- cur.sim + fsim$nsims
          
          # payoff at t
          immPayoff <- option.payoff(change.grid, model)
          add.mean[l] = mean(fsim$payoff - immPayoff)
          
          if (length(immPayoff) == 1) {
            add.var[l] = 0.00001
          } else {
            add.var[l] = var(fsim$payoff - immPayoff) + 0.00001 # avoid unstable results %
          }
          
          if (method == 'hetgp') {
            newX = rbind(newX, change.grid)
            newY = c(newY, fsim$payoff - immPayoff)
          }
        }
        
        y_new = (all.X[idx_diff, model$dim+1] * batch_matrix[idx_diff, i] + add.mean * r_seq_diff) / r_batch[idx_diff]
        var_new = (all.X[idx_diff, model$dim+2] * (batch_matrix[idx_diff, i] - 1) + batch_matrix[idx_diff, i] * all.X[idx_diff, model$dim+1] ^ 2 + add.var * (r_seq_diff - 1) + r_seq_diff * add.mean ^ 2) / (batch_matrix[idx_diff, i] - 1)
        all.X[idx_diff, model$dim + 1] = y_new
        all.X[idx_diff, model$dim + 2] = var_new
        batch_matrix[1:(n - 1), i] = r_batch
        
        running_budget = running_budget + sum(r_seq_diff)
        
        if (k %in% update.kernel.iters) {
          if (method == "km")
            fits[[i]] <- km(y~0, design=data.frame(x=all.X[1:n - 1, 1:model$dim]), response=data.frame(y=all.X[1:n - 1, model$dim+1]),
                            noise.var=all.X[1:n - 1,model$dim+2]/r_batch, covtype=model$kernel.family, coef.trend=0, coef.cov=model$km.cov,
                            coef.var=model$km.var, control=list(trace=F), 
                            lower=model$min.lengthscale,upper=model$max.lengthscale)
          if (method == "trainkm") {
            fits[[i]] <- km(y~0, design=data.frame(x=all.X[1:n - 1, 1:model$dim]), response=data.frame(y=all.X[1:n - 1, model$dim+1]),
                            noise.var=all.X[1:n - 1,model$dim+2]/r_batch, covtype=model$kernel.family, control=list(trace=F), 
                            lower=model$min.lengthscale, upper=model$max.lengthscale)
            theta.fit[i,n-model$init.size+1,] <- coef(fits[[i]])$range
          }
          if (method == "hetgp") {
            fits[[i]] <- update(object=fits[[i]], Xnew=newX, Znew=newY, method="mixed",
                                lower = model$min.lengthscale, upper=model$max.lengthscale)
            theta.fit[i,n-model$init.size+1,] <- fits[[i]]$theta
          }
          if (method == "homtp") {
            fits[[i]] <- mleHomTP(X=all.X[1:n - 1, 1:model$dim], Z=all.X[1:n - 1, model$dim+1],noiseControl=list(nu_bounds=c(2.01,5)),
                                  lower=model$min.lengthscale,upper=model$max.lengthscale, covtype=model$kernel.family)
            theta.fit[i,n-model$init.size+1,] <- fits[[i]]$theta
          }
        } else {
          if (method == "km" | method == "trainkm")
            fits[[i]] <- update(fits[[i]], newX=all.X[idx_diff, 1:model$dim, drop = F], newy=y_new,
                                newnoise=var_new / r_batch, newX.alreadyExist=T, cov.re=F)
          if (method == "hetgp")
            fits[[i]] <- update(object=fits[[i]], Xnew=newX, Znew=newY, maxit = 0, method = "quick")
          if (method == "homtp")
            fits[[i]] <- update(object=fits[[i]], Xnew=all.X[idx_diff, 1:model$dim, drop = F], Znew=y_new, maxit = 0)
        }
      } else { # add a new input
        #add.grid <- ei.cands[sample(dim(ei.cands)[1],model$n.propose,replace=T, prob=ei.weights),,drop=F]
        # select site with highest EI score
        
        add.grid <- matrix(rep(add.grid[1, ,drop=F], r_batch), nrow = r_batch, byrow=T)  #batch
        
        # compute corresponding y-values
        fsim <- forward.sim.policy( add.grid,M-i,fits[i:M],model,offset=0)
        cur.sim <- cur.sim + fsim$nsims
        
        immPayoff <- model$payoff.func(add.grid, model)
        add.mean <- mean(fsim$payoff - immPayoff)
        if (r_batch == 1) {
          add.var <- 0.00001
        } else {
          add.var <- var(fsim$payoff - immPayoff)
        }
        all.X[n,] <- c(add.grid[1,],add.mean,add.var)
        batch_matrix[n, i] <- r_batch
        
        if (k %in% update.kernel.iters) {
          if (method == "km")
            fits[[i]] <- km(y~0, design=data.frame(x=all.X[1:n,1:model$dim]), response=data.frame(y=all.X[1:n,model$dim+1]),
                            noise.var=all.X[1:n,model$dim+2]/r_batch, covtype=model$kernel.family, coef.trend=0, coef.cov=model$km.cov,
                            coef.var=model$km.var, control=list(trace=F), 
                            lower=model$min.lengthscale,upper=model$max.lengthscale)
          if (method == "trainkm") {
            fits[[i]] <- km(y~0, design=data.frame(x=all.X[1:n,1:model$dim]), response=data.frame(y=all.X[1:n,model$dim+1]),
                            noise.var=all.X[1:n,model$dim+2]/r_batch, covtype=model$kernel.family, control=list(trace=F), 
                            lower=model$min.lengthscale, upper=model$max.lengthscale)
            theta.fit[i,n-model$init.size+1,] <- coef(fits[[i]])$range
          }
          if (method == "hetgp") {
            fits[[i]] <- update(object=fits[[i]], Xnew=matrix(rep(add.grid[1,,drop=F], r_batch), nrow = r_batch, byrow = T), Znew=fsim$payoff-immPayoff,method="mixed",
                                lower = model$min.lengthscale, upper=model$max.lengthscale)
            theta.fit[i,n-model$init.size+1,] <- fits[[i]]$theta
          }
          if (method == "homtp") {
            fits[[i]] <- mleHomTP(X=all.X[1:n,1:model$dim], Z=all.X[1:n,model$dim+1],noiseControl=list(nu_bounds=c(2.01,5)),
                                  lower=model$min.lengthscale,upper=model$max.lengthscale, covtype=model$kernel.family)
            theta.fit[i,n-model$init.size+1,] <- fits[[i]]$theta
          }
        } else {
          if (method == "km" | method == "trainkm")
            fits[[i]] <- update(fits[[i]], newX=add.grid[1,,drop=F], newy=add.mean,
                                newnoise=add.var / r_batch,  cov.re=F)
          if (method == "hetgp")
            fits[[i]] <- update(object=fits[[i]], Xnew=matrix(rep(add.grid[1,,drop=F], r_batch), nrow = r_batch, byrow = T), Znew=fsim$payoff-immPayoff, maxit = 0, method = "quick")
          if (method == "homtp")
            fits[[i]] <- update(object=fits[[i]], Xnew=add.grid[1,,drop=F], Znew=add.mean, maxit = 0)
        }
        
        running_budget = running_budget + sum(r_batch)
        n = n + 1
      }
      
      # resample the candidate set
      ei.cands <- lhs(model$cand.len, lhs.rect)
      ei.cands <- ei.cands[ model$payoff.func( ei.cands,model) > 0,,drop=F]
      if (n >= model$seq.design.size || running_budget >= model$total.budget) {
        add.more.sites <- FALSE
      }
      k <- k+1
    }
    budget.used[i] <- n
  }
  
  return (list(fit=fits,timeElapsed=Sys.time()-t.start,nsims=cur.sim,empLoss=emp.loss,
               ndesigns=budget.used,theta=theta.fit, batches = batch_matrix))
}
