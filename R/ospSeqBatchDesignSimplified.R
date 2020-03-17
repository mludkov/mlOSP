#################
#' @title Adaptive Batch design for optimal stopping (simplified version)
#'
#' @details Implements the adaptive batching strategy defined in batch.heuristic with model defined in method. 
#' @param method: \code{km} to select the GP emulator to apply
#' @export
#' @return a list containing:
#' \code{fit} a list of fitted response surfaces
#' \code{ndesigns}: number of design size k_T
#' \code{batches}: matrix of replications r_i
osp.seq.batch.design.simplified <- function(model, method="km")
{
  #################
  #' @title Adaptive Batch design for optimal stopping (simplified version)
  #'
  #' @details Implements the adaptive batching strategy defined in batch.heuristic with model defined in method. 
  #' @param method: \code{km} to select the GP emulator to apply
  #' @export
  #' @return a list containing:
  #' \code{fit} a list of fitted response surfaces
  #' \code{ndesigns}: number of design size k_T
  #' \code{batches}: matrix of replications r_i
  #' 
  M <- model$T/model$dt
  
  # replication range in absur
  r_lower <- model$r.cand[1]
  r_upper <- min(model$r.cand[length(model$r.cand)], 0.1 * model$total.budget)
  r_interval <- seq(r_lower, r_upper, length = 1000)
  theta_for_optim <- c(0.1371, 0.000815, 1.9871E-6)  # c_{ovh} in eq. (13)
  batch_matrix <- matrix(rep(0, M*model$seq.design.size), ncol=M)
  
  # parameter in adsa and ddsa
  c_batch <- 20 / model$dim
  
  fits <- list()   # list of emulator objects at each step
  pilot.paths <- list()
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
  
  ############ step back in time
  for (i in (M-1):1)
  {
    all.X <- matrix( rep(0, (model$dim+2)*model$seq.design.size), ncol=model$dim+2)
    
    # construct the input domain where candidates will be looked for
    if (is.null(model$lhs.rect)) {
      model$lhs.rect <- 0.02
    }
    if (length(model$lhs.rect) == 1) {
      lhs.rect <- matrix(rep(0, 2*model$dim), ncol=2)
      # create a box using empirical quantiles of the init.grid cloud
      for (jj in 1:model$dim)
        lhs.rect[jj,] <- quantile( init.grid[,jj], c(model$lhs.rect, 1-model$lhs.rect) )
    } else { # already specified
      lhs.rect <- model$lhs.rect
    }
    
    # Candidate grid of potential NEW sites to add. Will be ranked using the EI acquisition function
    # only keep in-the-money sites
    ei.cands <- lhs( model$cand.len, lhs.rect )  # from tgp package
    ei.cands <- ei.cands[ model$payoff.func( ei.cands,model) > 0,,drop=F]
    
    # initial design
    init.grid <- model$init.grid
    K0 <- dim(init.grid)[1]
    # initial conditions for all the forward paths: replicated design with batch.nrep
    big.grid <- init.grid[ rep(1:K0, model$batch.nrep),]
    fsim <- forward.sim.policy( big.grid, M-i,fits[i:M],model, compact=T, offset=0)
    
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
    fits[[i]] <- km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[1:k,model$dim+1]),
                    noise.var=all.X[1:k,model$dim+2]/model$batch.nrep, covtype=model$kernel.family,
                    control=list(trace=F), lower=model$min.lengthscale, upper=model$max.lengthscale)
    
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
      pred.cands <- predict(fits[[i]],data.frame(x=ei.cands), type="UK")
      cand.mean <- pred.cands$mean
      cand.sd <- pred.cands$sd
      nug <- sqrt(mean(all.X[1:(k-1), model$dim+2]))
      nug <- rep(nug, length(cand.mean))
      
      losses <- cf.el(cand.mean, cand.sd)
      
      # multiply by weights based on the distribution of pilot paths
      dX_1 <- pilot.paths[[i]]; dX_2 <- ei.cands
      for (dd in 1:model$dim) {
        dX_1[,dd] <- dX_1[,dd]/(lhs.rect[dd,2]-lhs.rect[dd,1])
        dX_2[,dd] <- dX_2[,dd]/(lhs.rect[dd,2]-lhs.rect[dd,1])
      }
      # from package lagp
      ddx <- distance( dX_1, dX_2)
      x.dens <- apply( exp(-ddx*dim(ei.cands)[1]*0.01), 2, sum)
      
      # use active learning measure to select new sites and associated replication
      if (model$ei.func == 'absur') {
        overhead = 3 * model$dim * CalcOverhead(theta_for_optim[1], theta_for_optim[2], theta_for_optim[3], k + 1)
        al.weights <- cf.absur(cand.mean, cand.sd, nug, r_interval, overhead, model$t0)
        
        # select site and replication with highest EI score
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
        # use active learning measure to select new sites
        adaptive.gamma <- (quantile(cand.mean, 0.75, na.rm = T) - quantile(cand.mean, 0.25, na.rm = T))/mean(cand.sd)
        al.weights <- cf.smcu(cand.mean, cand.sd, adaptive.gamma)
        x.dens2 <- dlnorm( ei.cands[,1], meanlog=log(model$x0[1])+(model$r-model$div - 0.5*model$sigma[1]^2)*i*model$dt,
                           sdlog = model$sigma[1]*sqrt(i*model$dt) )
        x.dens2 <- x.dens2*dlnorm( ei.cands[,2], meanlog=log(model$x0[2])+(model$r-model$div-0.5*model$sigma[2]^2)*i*model$dt,
                                   sdlog = model$sigma[2]*sqrt(i*model$dt) )
        
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
        if (model$batch.heuristic == 'adsa') {
          adsa_batch <- batch.adsa(fits[[i]], batch_matrix[1:(n - 1), i], ei.cands, x.dens2, add.grid, r0, nug, method)
          add.grid <- adsa_batch$x_optim
          r_batch <- adsa_batch$r_optim
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
        
        ids <- seq(1, length(r_seq_diff))
        ids_rep <- unlist(mapply(rep, ids, r_seq_diff))
        
        newX <- all.X[idx_diff, 1:model$dim]
        newX <- matrix(unlist(mapply(rep, newX, r_seq_diff)), ncol = model$dim)
        
        # compute corresponding y-values
        fsim <- forward.sim.policy(newX, M-i, fits[i:M], model, compact=T, offset=0)
        
        # payoff at t
        immPayoff <- model$payoff.func(newX, model)
        newY <- fsim$payoff - immPayoff
        
        add.mean <- tapply(newY, ids_rep, mean)
        add.var <- tapply(newY, ids_rep, var)
        add.var[is.na(add.var)] <- 0.0000001
        
        y_new = (all.X[idx_diff, model$dim+1] * batch_matrix[idx_diff, i] + add.mean * r_seq_diff) / r_batch[idx_diff]
        var_new = (all.X[idx_diff, model$dim+2] * (batch_matrix[idx_diff, i] - 1) + batch_matrix[idx_diff, i] * all.X[idx_diff, model$dim+1] ^ 2 + add.var * (r_seq_diff - 1) + r_seq_diff * add.mean ^ 2) / (batch_matrix[idx_diff, i] - 1)
        all.X[idx_diff, model$dim + 1] = y_new
        all.X[idx_diff, model$dim + 2] = var_new
        batch_matrix[1:(n - 1), i] = r_batch
        
        running_budget = running_budget + sum(r_seq_diff)
        
        if (k %in% update.kernel.iters) {
            fits[[i]] <- km(y~0, design=data.frame(x=all.X[1:n - 1, 1:model$dim]), response=data.frame(y=all.X[1:n - 1, model$dim+1]),
                            noise.var=all.X[1:n - 1,model$dim+2]/r_batch, covtype=model$kernel.family, control=list(trace=F),
                            lower=model$min.lengthscale, upper=model$max.lengthscale)
        } else {
            fits[[i]] <- update(fits[[i]], newX=all.X[idx_diff, 1:model$dim, drop = F], newy=y_new,
                                newnoise=var_new / r_batch, newX.alreadyExist=T, cov.re=F)
        }
      } else { # add a new input
        add.grid <- matrix(rep(add.grid[1, ,drop=F], r_batch), nrow = r_batch, byrow=T)  #batch
        
        # compute corresponding y-values
        fsim <- forward.sim.policy( add.grid,M-i,fits[i:M],model,offset=0)
        
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
          fits[[i]] <- km(y~0, design=data.frame(x=all.X[1:n,1:model$dim]), response=data.frame(y=all.X[1:n,model$dim+1]),
                          noise.var=all.X[1:n,model$dim+2]/r_batch, covtype=model$kernel.family, control=list(trace=F),
                          lower=model$min.lengthscale, upper=model$max.lengthscale)
        } else {
          fits[[i]] <- update(fits[[i]], newX=add.grid[1,,drop=F], newy=add.mean,
                              newnoise=add.var / r_batch,  cov.re=F)
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
  
  return (list(fit=fits, ndesigns=budget.used, batches = batch_matrix))
}
