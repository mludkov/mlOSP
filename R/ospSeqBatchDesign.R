#################
#' @title Adaptive Batching designs for optimal stopping 
#' 
#' @description Sequential experimental design for optimal stopping problems with several 
#' adaptive batching heuristics based on Lyu & Ludkovski (2020+)
#'
#' @details Implements the adaptive batching strategy defined in \code{mode$batch.heuristic}. 
#' Calls \code{lhs}  from library \pkg{tgp}. Possible batch heuristics are:
#' \itemize{
#'   \item \code{fb}: [Default] fixed batch amounts (essentially same as \link{osp.seq.design})
#'   \item \code{mlb}: Multi-level batching; relies on \code{model$r.cand}
#'   \item \code{rb}: Ratchet batching; relies on \code{model$r.cand}
#'   \item \code{absur}: Adaptive batching with Stepwise Uncertainty Reduction; relies on \code{model$t0}
#'   \item \code{adsa}: Adaptive Design with Sequential Allocation
#'   \item \code{ddsa}: Deterministic ADSA that alternates between adding a new input site and allocating
#'   to existing sites
#' }
#' 
#' All heuristics also require specifying the acquisition function for expected improvement criterion 
#' via \code{model$ei.func}, see \link{osp.seq.design}
#' 
#' @param model  a list containing all the model parameters.
#' 
#' @param method either \code{km} or \code{hetgp} to select the GP emulator to apply
#' @param t0 parameter \code{t0} for \code{ABSUR} heuristic [Default value is 0.01]
#' @param is.gbm flag to indicate whether the underlying simulator is independent log-normals (used 
#' as part of density computation for integrated EI criteria) [Default FALSE]
#' @export
#' @return a list containing:
#' \itemize{
#' \item \code{fit} a list of fitted response surfaces
#' \item \code{timeElapsed} vector of time costs for each round
#' \item \code{nsims} total number of 1-step \code{model$sim.func} calls
#' \item \code{empLoss} vector of empirical losses
#' \item \code{ndesigns}: number of unique designs k_T
#' \item \code{batches}: matrix of replications r_i, indexed by time-steps and by sequential rounds
#' }
#' @references 
#' M. Ludkovski, X. Lyu (2020+) Adaptive Batching for Gaussian Process Surrogates with Application 
#' in Noisy Level Set Estimation, <http://arxiv.org/abs/2003.08579>
#' 
#' @seealso [mlOSP::osp.seq.design]
#' 
#' @examples
#' sob30 <- randtoolbox::sobol(55, d=2)  # construct a space-filling initial design
#' sob30 <- sob30[ which( sob30[,1] + sob30[,2] <= 1) ,]  
#' sob30 <- 25+30*sob30 
#' model2d <- list(x0 = rep(40,2),K=40,sigma=rep(0.2,2),r=0.06,
#'  div=0,T=1,dt=0.04,dim=2,sim.func=sim.gbm, 
#'  payoff.func=put.payoff, look.ahead=1, pilot.nsims=1000,
#'  cand.len=1000,max.lengthscale=c(40,40),min.lengthscale=c(3,3),
#'  seq.design.size=50,batch.nrep=25,total.budget=1000,init.size=30,
#'  init.grid=sob30, kernel.family="gauss",update.freq=5,
#'  r.cand=c(20, 30,40,50,60, 80, 120, 160))
#' set.seed(11)
#' require(tgp)
#' require(DiceKriging)
#' require(laGP)
#' require(ks)
#' model2d$batch.heuristic <- 'adsa'
#' model2d$ei.func <- 'amcu'
#' oos.obj.adsa <- osp.seq.batch.design(model2d,method="trainkm")
#' # not run: plt.2d.surf.with.batch(oos.obj.adsa$fit[[10]],25)
osp.seq.batch.design <- function(model, method="km", t0 = 0.01, is.gbm=FALSE)
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
  if (is.null(model$c.batch))  # parameter for new replicates in adsa and ddsa
        model$c.batch = 20 / model$dim
  if (is.null(model$r.cand)) {
    model$r.cand = c(model$batch.nrep, model$batch.nrep)
  }
  r_lower = model$r.cand[1]
  r_upper = min(model$r.cand[length(model$r.cand)], 0.1 * model$total.budget)
  r_interval = seq(r_lower, r_upper, length = 1000)
  theta_for_optim = c(0.1371, 0.000815, 1.9871E-6)  # c_{oh} in eq. (13)
  batch_matrix <- matrix(rep(0, M*model$seq.design.size), ncol=M)


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
    ei.cands <- tgp::lhs( model$cand.len, lhs.rect )  # from tgp package
    ei.cands <- ei.cands[ model$payoff.func( ei.cands,model) > 0,,drop=F]


    # initial design
    if (is.null(model$init.grid)) {
      init.grid <- tgp::lhs(model$init.size, lhs.rect)
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
      fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[1:k,model$dim+1]),
                      noise.var=all.X[1:k,model$dim+2]/model$batch.nrep,
                      control=list(trace=F), lower=model$min.lengthscale, covtype=model$kernel.family,
                      #nugget.estim= TRUE,
                      coef.trend=0, coef.cov=model$km.cov, coef.var=model$km.var)
    }
    if (method == "trainkm") {
      fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), 
                      response=data.frame(y=all.X[1:k,model$dim+1]),
                      nugget.estim=TRUE, #
                      #noise.var=all.X[1:k,model$dim+2]/model$batch.nrep, 
                      covtype=model$kernel.family,
                      control=list(trace=F), lower=model$min.lengthscale, upper=model$max.lengthscale)
      #if (coef(fits[[i]])$sd2 < 1e-10)
      #  fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), 
      #                               response=data.frame(y=all.X[1:k,model$dim+1]),
      #                               nugget.estim=TRUE, #
      #                               #noise.var=all.X[1:k,model$dim+2]/model$batch.nrep, 
      #                               covtype=model$kernel.family,
      #                               control=list(trace=F), lower=model$min.lengthscale, upper=model$max.lengthscale)
      
      theta.fit[i,k-model$init.size+1,] <- coef(fits[[i]])$range
    }
    if (method == "homtp") {
      fits[[i]] <- hetGP::mleHomTP(X=all.X[1:k,1:model$dim], Z=all.X[1:k,model$dim+1],
                                   lower=model$min.lengthscale,upper=model$max.lengthscale,
                                   covtype=model$kernel.family, noiseControl=list(nu_bounds=c(2.01,5)))
      theta.fit[i,k-model$init.size+1,] <- fits[[i]]$theta
    }
    if (method == "homgp") {
      big.payoff <- model$payoff.func(big.grid,model)
      hetData <- hetGP::find_reps(big.grid, fsim$payoff-big.payoff)
      fits[[i]] <- hetGP::mleHomGP(X = list(X0=hetData$X0, Z0=hetData$Z0,mult=hetData$mult), Z= hetData$Z,
                                   lower = model$min.lengthscale, upper = model$max.lengthscale, 
                                   covtype=model$kernel.family)
      
      theta.fit[i,k-model$init.size+1,] <- fits[[i]]$theta
    }
    if (method== "hetgp") {
      big.payoff <- model$payoff.func(big.grid,model)
      hetData <- hetGP::find_reps(big.grid, fsim$payoff-big.payoff)
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
    
    if (i == 43)
       browser()

    # active learning loop:
    # to do it longer for later points to minimize error build-up use: *(1 + i/M)
    while(add.more.sites)
    {
      # predict on the candidate grid. Need predictive GP mean, posterior GP variance and similation Stdev
      if (method == "km" | method == "trainkm") {
        pred.cands <- predict(fits[[i]],data.frame(x=ei.cands), type="UK")
        cand.mean <- pred.cands$mean
        cand.sd <- pred.cands$sd
        nug <- sqrt(coef(fits[[i]])$nug) 
        #nug <- sqrt(mean(all.X[1:(k-1), model$dim+2]))
        nug <- rep(nug, length(cand.mean)) # noise.var=all.X[1:n,model$dim+2]/batch_matrix[1:n,i]
      } else {   # hetGP/homTP/homGP
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
      ddx <- laGP::distance( dX_1, dX_2)
      x.dens <- apply( exp(-ddx*dim(ei.cands)[1]*0.01), 2, sum)
      emp.loss[i,k-model$init.size] <- sum(losses*x.dens)/sum(x.dens)
      #if (is.null(model$el.thresh) == FALSE) { # stop if below the Emp Loss threshold
      #  if (emp.loss[i,k-model$init.size] < model$el.thresh) {
      #    add.more.sites <- FALSE
      #  }
      #}

      # use active learning measure to select new sites and associated replication
      if (model$ei.func == 'absur') {
          overhead = 3 * model$dim * CalcOverhead(theta_for_optim[1], theta_for_optim[2], theta_for_optim[3], k + 1)
          al.weights <- cf.absur(cand.mean, cand.sd, nug, r_interval, overhead, t0)

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
          # Integrated contour uncertainty with weights based on *Hard-coded* log-normal density and create testing points
          if (is.gbm) {
            x.dens2 <- dlnorm( ei.cands[,1], meanlog=log(model$x0[1])+(model$r-model$div - 0.5*model$sigma[1]^2)*i*model$dt,
                               sdlog = model$sigma[1]*sqrt(i*model$dt) )
            jdim <- 2
            while (jdim <= model$dim) {
              x.dens2 <- x.dens2*dlnorm( ei.cands[,jdim], 
                                         meanlog=log(model$x0[jdim])+(model$r-model$div-0.5*model$sigma[jdim]^2)*i*model$dt,
                                         sdlog = model$sigma[jdim]*sqrt(i*model$dt) )
              jdim <- jdim+1
            }
          } else {
            x.dens2 <- x.dens }

          if (model$ei.func == 'icu' & method == "hetgp") {  # only works with hetGP
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
        r0 = min(model$total.budget - running_budget, round(model$c.batch*sqrt(k)))
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
          #browser()
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
      if (length(r_batch) == 1 && is.numeric(r_batch) && r_batch == 0) {
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
        cur.sim <- cur.sim + fsim$nsims

        # payoff at t
        immPayoff <- model$payoff.func(newX, model)
        newY <- fsim$payoff - immPayoff

        add.mean <- tapply(newY, ids_rep, mean)
        add.var <- tapply(newY, ids_rep, var)
        add.var[is.na(add.var)] <- 0.0000001

        y_new = (all.X[idx_diff, model$dim+1] * batch_matrix[idx_diff, i] + add.mean * r_seq_diff) / r_batch[idx_diff]
        var_new = (all.X[idx_diff, model$dim+2] * (batch_matrix[idx_diff, i] - 1) + 
                     batch_matrix[idx_diff, i] * all.X[idx_diff, model$dim+1] ^ 2 +
                     add.var * (r_seq_diff - 1) + r_seq_diff * add.mean ^ 2) / (batch_matrix[idx_diff, i] - 1)
        all.X[idx_diff, model$dim + 1] = y_new
        all.X[idx_diff, model$dim + 2] = var_new
        batch_matrix[1:(n - 1), i] = r_batch

        running_budget = running_budget + sum(r_seq_diff)

        if (k %in% update.kernel.iters) {
          if (method == "km")
            fits[[i]] <- km(y~0, design=data.frame(x=all.X[1:n - 1, 1:model$dim]), 
                            response=data.frame(y=all.X[1:n - 1, model$dim+1]),
                            noise.var=all.X[1:n - 1,model$dim+2]/batch_matrix[1:(n - 1), i], 
                            covtype=model$kernel.family, coef.trend=0, coef.cov=model$km.cov,
                            coef.var=model$km.var, control=list(trace=F),
                            lower=model$min.lengthscale,upper=model$max.lengthscale)
          if (method == "trainkm") {
            fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=all.X[1:n - 1, 1:model$dim]), 
                            response=data.frame(y=all.X[1:n - 1, model$dim+1]),
                            nugget.estim=TRUE, 
                            #noise.var=all.X[1:n - 1,model$dim+2]/batch_matrix[1:(n - 1), i], 
                            covtype=model$kernel.family, control=list(trace=F),
                            lower=model$min.lengthscale, upper=model$max.lengthscale)
            theta.fit[i,n-model$init.size+1,] <- coef(fits[[i]])$range
          }
          if (method == "hetgp" | method == "homgp") {
            fits[[i]] <- update(object=fits[[i]], Xnew=newX, Znew=newY, method="mixed",
                                lower = model$min.lengthscale, upper=model$max.lengthscale)
            theta.fit[i,n-model$init.size+1,] <- fits[[i]]$theta
          }
          if (method == "homtp") {
            fits[[i]] <- hetGP::mleHomTP(X=all.X[1:n - 1, 1:model$dim], Z=all.X[1:n - 1, model$dim+1],
                                         noiseControl=list(nu_bounds=c(2.01,5)),
                                  lower=model$min.lengthscale,upper=model$max.lengthscale, covtype=model$kernel.family)
            theta.fit[i,n-model$init.size+1,] <- fits[[i]]$theta
          }
        } else {
          if (method == "km" | method == "trainkm")
            fits[[i]] <- update(fits[[i]], newX=all.X[idx_diff, 1:model$dim, drop = F], newy=y_new,
                                #newnoise=var_new / r_batch, 
                                newX.alreadyExist=T, cov.re=F)
          if (method == "hetgp" | method == "homgp")
            fits[[i]] <- update(object=fits[[i]], Xnew=newX, Znew=newY, maxit = 0, method = "quick")
          if (method == "homtp" )
            fits[[i]] <- update(object=fits[[i]], Xnew=all.X[idx_diff, 1:model$dim, drop = F], Znew=y_new, maxit = 0)
        }
      } else { # add a new input
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
            fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=all.X[1:n,1:model$dim]), 
                                         response=data.frame(y=all.X[1:n,model$dim+1]),
                            noise.var=all.X[1:n,model$dim+2]/batch_matrix[1:n,i], covtype=model$kernel.family, 
                            coef.trend=0, coef.cov=model$km.cov,
                            coef.var=model$km.var, control=list(trace=F),
                            lower=model$min.lengthscale,upper=model$max.lengthscale)
          if (method == "trainkm") {
            fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=all.X[1:n,1:model$dim]), 
                                         response=data.frame(y=all.X[1:n,model$dim+1]),
                            nugget.estim=TRUE,
                            #noise.var=all.X[1:n,model$dim+2]/batch_matrix[1:n,i], 
                            covtype=model$kernel.family, control=list(trace=F),
                            lower=model$min.lengthscale, upper=model$max.lengthscale)
            theta.fit[i,n-model$init.size+1,] <- coef(fits[[i]])$range
          }
          if (method == "hetgp" | method == "homgp") {
            fits[[i]] <- update(object=fits[[i]], Xnew=matrix(rep(add.grid[1,,drop=F], r_batch), nrow = r_batch, byrow = T), 
                                Znew=fsim$payoff-immPayoff,method="mixed",
                                lower = model$min.lengthscale, upper=model$max.lengthscale)
            theta.fit[i,n-model$init.size+1,] <- fits[[i]]$theta
          }
          if (method == "homtp") {
            fits[[i]] <- hetGP::mleHomTP(X=all.X[1:n,1:model$dim], Z=all.X[1:n,model$dim+1],
                                         noiseControl=list(nu_bounds=c(2.01,5)),
                                  lower=model$min.lengthscale,upper=model$max.lengthscale, covtype=model$kernel.family)
            theta.fit[i,n-model$init.size+1,] <- fits[[i]]$theta
          }
      
        } else {
          if (method == "km" | method == "trainkm")
            fits[[i]] <- update(fits[[i]], newX=add.grid[1,,drop=F], newy=add.mean,
                                #newnoise=add.var / r_batch,  
                                cov.re=F)
          if (method == "hetgp" | method == "homgp")
            fits[[i]] <- update(object=fits[[i]], Xnew=matrix(rep(add.grid[1,,drop=F], r_batch), nrow = r_batch, byrow = T), 
                                Znew=fsim$payoff-immPayoff, maxit = 0, method = "quick")
          if (method == "homtp" )
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
               ndesigns=budget.used, theta=theta.fit, batches = batch_matrix))
}
