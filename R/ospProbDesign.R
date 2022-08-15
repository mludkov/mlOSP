#############################
#' RMC using probabilistic design: backpropagation along fixed set of paths (a la Longstaff-Schwartz).
#' All designs are kept in memory. By default produces only an in-sample estimate. Use in conjuction
#' with \code{\link{forward.sim.policy}} to generate out-of-sample price estimates.
#' 
#' @title Longstaff-Schwartz RMC algorithm with a variety of regression methods.
#' 
#' @param model defines the simulator and reward model, with the two main model hooks being. 
#' Initial condition is \code{model$x0}. Can be either a vector of length \code{model$dim} 
#' or a vector of length \code{model$dim}*N 
#' option payoff  \code{payoff.func} (plus parameters) and stochastic simulator \code{sim.func} (plus parameters)
#' @param N is the number of paths
#' @param subset To have out-of-sample paths, specify \code{subset} (e.g 1:1000) to use for testing.
#' By default everything is in-sample
#'
#' @param method a string specifying regression method to use
#' \itemize{
#'  \item spline: Smoothing splines \code{smooth.spline} from \pkg{base}. Only works \emph{in 1D}.
#'  Requires number of knots \code{nk}.
#'  \item randomforest: (from \pkg{randomForest} package) requires \code{rf.maxnode}
#'  and \code{rf.ntree} (number of trees) model parameters
#'  \item loess: local polynomial regression. Only works in \emph{1D or 2D}, 
#'  requires \code{lo.span} model parameter
#'  \item earth: multivariate adaptive regression splines (MARS) using \pkg{earth} package.
#'  Requires \code{earth.deg} (interaction degree), \code{earth.nk} (max number of terms to keep),
#'  \code{earth.thresh} params
#'  \item rvm: relevance vector machine from \pkg{kernlab} package. Optional \code{rvm.kernel}
#'  model parameter to decide which kernel family to utilize. Default kernel is rbfdot
#'  \item npreg: kernel regression using \pkg{np} package. Can optionally provide \code{np.kertype} 
#'  (default is "gaussian"); \code{np.regtype} (default is "lc"); \code{np.kerorder} (default is 2)
#'  \item nnet: neural network using \pkg{nnet}. This is a single-layer neural net. Specify a scalar \code{nn.nodes} 
#'  to describe the number of nodes at the hidden layer
#'  \item lagp: local approximate Gaussian Process regression using \pkg{lagp} package. Can 
#'  optionally provide \code{lagp.type} (default is "alcray" which is fastest, other choices are "alc" 
#'  and "mspe") that determines how the local design is constructed, and \code{lagp.end} which determines
#'  how many inputs are in the above local design. 
#'  \item dynatree: dynamic trees using \pkg{dynaTree}. Requires \code{dt.type} ("constant" 
#'  or "linear" fits at the leafs), \code{dt.Npart} (number of trees), \code{dt.minp} (minimum size
#'  of each partition) and \code{dt.ab} (the tree prior parameter) model parameters.
#'  \item lm [Default]: linear global regression using \code{model$bases} (required) basis functions 
#'  (+ constant) which is a function pointer.
#'  }
#' @export
#' @return a list containing
#' \itemize{
#' \item \code{fit} a list containing all the models generated at each time-step. \code{fit[[1]]} is the emulator
#' at \eqn{t=\Delta t}, the last one is \code{fit[[M-1]]} which is emulator for \eqn{T-\Delta t}.
#' \item \code{val}: the in-sample pathwise rewards
#' \item \code{test}: the out-of-sample pathwise rewards
#' \item \code{p}: the final price (2-vector for in/out-of-sample)
#' \item \code{timeElapsed} total running time in seconds, based on \code{Sys.time}
#' }
#' @details
#'  Works with a probabilistic design that requires storing all paths in memory. Specifying \code{subset}
#'  allows to compute in parallel with the original computation an out-of-sample estimate of the value function
#'  
#'  Calls \code{model$payoff.func}, so the latter must be set prior to calling.
#'  Also needs \code{model$dt}, \code{model$T} for simulation and \code{model$r} for discounting
#'  
#'  Calls \code{model$sim.func} to generate forward paths
#'  
#'  Emulator is trained only on paths where payoffs are strictly positive
#' @examples
#' set.seed(1)
#' model2d <- list(look.ahead=1,K=40,x0=rep(40,2),sigma=rep(0.2,2),r=0.06,
#'  div=0, T=1,dt=0.04,dim=2, sim.func=sim.gbm, payoff.func=put.payoff)
#'  bas22 <- function(x) return(cbind(x[,1],x[,1]^2,x[,2],x[,2]^2,x[,1]*x[,2]))
#'  model2d$bases <- bas22
#'  prob.lm <- osp.prob.design(30000,model2d,method="lm",subset=1:15000)
#'  prob.lm$p
#'  # yields [1] 1.440918 1.482422


#------------
osp.prob.design <- function(N,model,subset=1:N,method="lm")
{
  M <- as.integer(round(model$T/model$dt))
  grids <- list()
  all.models <- list()
  if (is.null(subset))
    subset = 1:N

  # divide into in-sample and out-of-sample
  if (length(subset) < N)
    train <- (1:N)[-subset]
  else
    train <- 1:N

  # Build the designs based on a simulation of X_{1:T}
  if (length(model$x0) == model$dim)
     grids[[1]] <- model$sim.func( matrix(rep(model$x0, N), nrow=N,byrow=T), model, model$dt)
  else if (length(model$x0) == N*model$dim)
     grids[[1]] <- model$sim.func( matrix(rep(model$x0, N), nrow=N,byrow=T), model, model$dt)
  else
    warning("Length of the initial condition which must match N")
  for (i in 2:M)
    grids[[i]] <- model$sim.func( grids[[i-1]], model, model$dt)

  # initialize at T
  contValue <- exp(-model$r*model$dt)*model$payoff.func( grids[[M]], model)
  tau <- rep(model$T, N)
  t.start <- Sys.time()

  # Backward stepping in time
  # Estimate T(t,x)
  for (i in (M-1):1)
  {
    # forward predict
    immPayoff <- model$payoff.func(grids[[i]],model)
    
    # train only on in-the-money
    c.train <- train[which(immPayoff[train] > 0)]
    yVal <- contValue[c.train]-immPayoff[c.train]
    if (length(c.train) == 0) { # all paths are out-of-the-money
      contValue <- exp(-model$r*model$dt)*contValue
      next
    }

    ### CASES DEPENDING ON METHOD
    if (method == "spline" & ncol(grids[[i]]) == 1) { # only works in 1D
      if (is.null(model$nk) )
        stop("Missing parameters to pass to stats::smooth.spline function. Need model$nk")
      
      all.models[[i]] <- stats::smooth.spline( x=grids[[i]][c.train],y=yVal,
                                        nknots = model$nk)
      timingValue <- predict(all.models[[i]],grids[[i]])$y
    }

    if (method == "randomforest") {
      if (is.null(model$rf.ntree) | is.null(model$rf.maxnode))
        stop("Missing parameters to pass to randomForest function. Need rf.maxnode and rf.ntree")
      
      all.models[[i]] <-  randomForest::randomForest(x=grids[[i]][c.train,,drop=F],y=yVal,
                                       ntree=model$rf.ntree,replace=F,maxnode=model$rf.maxnode)
      #timingValue <- predict(all.models[[i]],grids[[i]],predict.all=T)$individual
      #stopProb <- apply( (timingValue < 0), 1, sum)/model$rf.ntree  # not used right now
      #timingValue <- apply(timingValue,1,mean)   # median or mean prediction could be used
      timingValue <- predict(all.models[[i]],grids[[i]])
    }

    if (method == "loess" & ncol(grids[[i]]) <= 2) { # LOESS only works in 1D or 2D
        if (is.null(model$lo.span) )
          stop("Missing parameters to pass to stats::loess function. Need model$lo.span")
      
        if (ncol(grids[[i]]) == 1) {
          all.models[[i]] <- stats::loess(y ~ x, data.frame(x=grids[[i]][c.train], y=yVal),
                                 span=model$lo.span, control = loess.control(surface = "direct"))
          stopProb <- predict(all.models[[i]], data.frame(x=grids[[i]]),se=TRUE)
        }
        if (ncol(grids[[i]]) ==2) {
          all.models[[i]] <- stats::loess(y ~ x1+x2, data.frame(x1=grids[[i]][c.train,1], x2=grids[[i]][c.train,2], y=yVal),
                                   span=model$lo.span, control = loess.control(surface = "direct"))
          stopProb <- predict(all.models[[i]], new=data.frame(x1=grids[[i]][,1],x2=grids[[i]][,2]))
        }
        timingValue <- stopProb$fit
        stopProb <- 1-pnorm( (stopProb$fit)/stopProb$se)
    }
    if (method == "earth") {  # Multivariate Adaptive Regression Splines
      if (is.null(model$earth.deg) | is.null(model$earth.thresh) | is.null(model$earth.nk))
         stop("Missing parameters to pass to earth::earth function. Need earth.deg, earth.nk, earth.thresh")
       all.models[[i]] <- earth::earth(x=grids[[i]][c.train,,drop=F],y=yVal,
                                degree=model$earth.deg,nk=model$earth.nk,thresh=model$earth.thresh)
       timingValue <- predict(all.models[[i]],grids[[i]])
    }
    if (method == "deepnet") {  # Neural Network via the deepnet library
      if (is.null(model$nn.layers) )
        stop("Missing parameters to pass to deepnet::nn.train function. Need nn.layers")
      all.models[[i]] <- deepnet::nn.train(x=grids[[i]][c.train,,drop=F],y=yVal,
                               hidden=model$nn.layers)
      timingValue <- deepnet::nn.predict(all.models[[i]],grids[[i]])
    }
    if (method == "nnet") {  # Neural Network via the nnet library
      all.models[[i]] <- nnet::nnet(x=grids[[i]][c.train,,drop=F],y=yVal,
                              size=model$nn.nodes, linout=TRUE, maxit=1000,trace=FALSE)
      timingValue <- predict(all.models[[i]],grids[[i]], type="raw")
    }
    if (method == "lagp") { # approximate GP
      if (is.null(model$lagp.type))
        lagp.type <- "alcray"
      else
        lagp.type <- model$lagp.type
      if (is.null(model$lagp.end))
        lagp.end <- 50
      else
        lagp.end <- model$lagp.end
      if (i==(M-1))
        startd <- NULL
      else
        # initial guess for lengthscale based on previous step, make sure it's within allowed range
        startd <- pmax(pmin(darg(NULL,grids[[i]][c.train, , drop = F])$max*0.9, mean(all.models[[i+1]]$mle$d)),
                       darg(NULL,grids[[i]][c.train, , drop = F])$min*1.05)
      # do not estimate the nugget
      all.models[[i]] <- aGP(X=grids[[i]][c.train,,drop=F], Z=yVal, XX=grids[[i]], g=NULL,
                             end=lagp.end, method=lagp.type,d=startd,
                             omp.threads=4,verb=0,Xi.ret=FALSE)
      timingValue <- all.models[[i]]$mean 
      all.models[[i]] <- list(x=grids[[i]][c.train,,drop=F], y=yVal,
                              mle=all.models[[i]]$mle,time=all.models[[i]]$time)
      class(all.models[[i]]) <- "agp"
    }
    if (method == "ligp") {
        scale_list <- find_shift_scale_stats(as.matrix(grids[[i]][c.train,,drop=F]), yVal)
        Xsc <- shift_scale(list(X=as.matrix(grids[[i]])), 
                           shift=scale_list$xshift, scale=scale_list$xscale)$X
  
        Ysc <- (yVal-scale_list$yshift)/scale_list$yscale
        lhs_design <- randomLHS(model$ligp.lhs,model$dim)
        Xmt <- scale_ipTemplate(Xsc, model$ligp.n, space_fill_design=lhs_design, method='qnorm')$Xm.t
        
        out <- liGP(XX=Xsc, X=Xsc[c.train,,drop=F], Y=Ysc, Xm=Xmt, N=model$ligp.n, theta=model$ligp.theta,
                    g=list(start=1e-4,min=1e-6, max=model$ligp.gmax))
        timingValue <- out$mean*scale_list$yscale + scale_list$yshift
        
        all.models[[i]] <-  list(scale=scale_list, Xsc=Xsc[c.train,,drop=F], Ysc = Ysc, Xmt=Xmt)
        class(all.models[[i]]) <- "ligp"
        
      
    }
    if (method =="lm") {  # Linear regression with specified basis functions
      matb <- model$bases(grids[[i]][c.train,,drop=F])
      all.models[[i]] <- lm(yVal ~ matb)
      lenn <- length(all.models[[i]]$coefficients)
      timingValue <-  all.models[[i]]$coefficients[1] +
        model$bases(grids[[i]]) %*% all.models[[i]]$coefficients[2:lenn]
    }
    if (method =="rvm") {   # Relevance Vector Machine Regression
      if (is.null(model$rvm.kernel))
          rvmk <- "rbfdot"
      else
          rvmk <- model$rvm.kernel
      all.models[[i]] <- kernlab::rvm(x=grids[[i]][c.train,,drop=F], y=yVal,kernel=rvmk)
      timingValue <- predict(all.models[[i]], new=grids[[i]])
    }
    if (method == "npreg") {
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

        all.models[[i]] <- np::npregbw(xdat = grids[[i]][c.train,,drop=F], ydat = yVal,
                           regtype=regtype, ckertype=kertype,ckerorder=kerorder)
        timingValue <- fitted(npreg(bws=all.models[[i]],exdat=grids[[i]]))
    }
    if (method=="dynatree") {
      if (is.null(model$dt.Npart) | is.null(model$dt.minp) | is.null(model$dt.ab) | is.null(model$dt.type))
        stop("Missing parameters to pass to dynaTree::dynaTree function. Need dt.type, dt.minp, dt.Npart, dt.ab")
      
      all.models[[i]] <- dynaTree::dynaTree(grids[[i]][c.train,,drop=F], yVal, 
                                            model=model$dt.type, N=model$dt.Npart,minp=model$dt.minp,
                                            ab=model$dt.ab ,verb=0)
                #icept="augmented",
      timingValue <- predict(all.models[[i]], grids[[i]],quants=F)$mean  
    }

    # paths on which stop right now
    stopNdx <- which( timingValue <= 0 & immPayoff > 0)
    contValue[stopNdx] <- immPayoff[stopNdx]
    tau[stopNdx] <- i*model$dt

    # else continue and discount
    contValue <- exp(-model$r*model$dt)*contValue
  }
  test <- NULL

  # in/sample and out-of-sample average at x0
  if (length(subset) < N & length(subset) > 0) {
      price <- c(mean(contValue[train]),mean(contValue[subset]))
      print(sprintf("in-sample v_0 %3f; and out-of-sample: %3f", price[1], price[2]))
      test <- contValue[subset]
  }
  else
    price <- mean(contValue)

  # returns a list containing
  # fit are all the models generated at each time-step, stored as a list
  # p is the final price (2-vector for in/out-of-sample)
  # val are the in-sample pathwise rewards
  # test are the out-of-sample pathwise rewards
  # timeElapsed: total running time
  return( list(fit=all.models,p=price, val=contValue[train], test=test,
               timeElapsed=Sys.time()-t.start))
}


####################################
#' RMC based on a batched non-adaptive design with a variety of regression methods
#'
#' @title Generic dynamic emulation of OSP with a non-sequential design
#' @param model a list containing all the model parameters, see Details. 
#' @param input.domain the domain of the emulator. Several options are available. Default in \code{NULL}
#' All the empirical domains rely on pilot paths generated using \code{pilot.nsims}>0 model parameter.
#' \itemize{
#' \item  NULL will use an empirical probabilistic design based on the pilot paths (default);
#' \item if a vector of length 2*model$dim then specifies the bounding rectangle
#' \item a single positive number, then build a bounding rectangle based on the \eqn{\alpha}-quantile of the pilot paths
#' \item a single negative number, then build a bounding rectangle based on the full range of the pilot paths
#' \item a vector specifies the precise design, used as-is (\emph{overrides design size})
#' }
#' @param method regression method to use (defaults to \code{km})
#' \itemize{
#' \item km: Gaussian process with fixed hyperparams  uses \pkg{DiceKriging} via \code{km} (default)
#' \item trainkm: GP w/trained hyperparams: use \pkg{DiceKriging} via \code{km}
#' \item mlegp Local approximate GP from the \pkg{laGP}. Requires
#' \item homgp Homoskedastic GP: use \pkg{hetGP} with  \code{mleHomGP}
#' \item hetgp Heteroskedastic GP: use \pkg{hetGP} with \code{mleHetGP}
#' \item spline: Smoothing Splines, use \code{smooth.spline} (only in 1D). Requires number of 
#' knots via \code{model$nk} 
#' \item loess: Local Regression: use \code{loess} with \code{lo.span} parameter (only in 1D or 2D)
#' \item rvm: Relevance Vector Machine: uses \code{rvm} from \pkg{kernlab}. Can optionally provide
#' \code{rvm.kernel} parameter (default is 'rbfdot')
#' \item npreg: kernel regression using \pkg{npreg} package. Can optionally provide \code{np.kertype} 
#'  (default is "gaussian"); \code{np.regtype} (default is "lc"); \code{np.kerorder} (default is 2)
#' \item lm: linear model from \pkg{stats}
#' }
#' @param inTheMoney.thresh which paths are kept, out-of-the-money is dropped.
#' Defines threshold in terms of \code{model$payoff.func}
#' @param stop.freq *experimental, currently disabled* frequency of stopping decisions (default is \code{model$dt}).
#' Can be used to stop less frequently.
#' @return a list containing:
#' \itemize{
#' \item \code{fit} a list containing all the models generated at each time-step. \code{fit[[1]]} is the emulator
#' at \eqn{t=\Delta t}, the last one is \code{fit[[M-1]]} which is emulator for \eqn{T-\Delta t}.
#' \item \code{timeElapsed}: total running time based on \code{Sys.time}
#' }
#' @details The design can be replicated through \code{batch.nrep} model parameter. Replication allows to use
#' nonparametric techniques which would be too expensive otherwise, in particular LOESS, GP and RVM.
#' All designs are restricted to in-the-money region, see \code{inTheMoney.thresh} parameter (modify at your own risk)
#' Thus, actual design size will be smaller than specified. By default, no forward evaluation is provided, ie the
#' method only builds the emulators. Thus, to obtain an actual estimate of the opton price
#' combine with \code{\link{forward.sim.policy}}.
#' 
#' @author Mike Ludkovski
#' @export
#' 
#' @examples
#' set.seed(1)
#' model2d <- list(K=40,x0=rep(40,2),sigma=rep(0.2,2),r=0.06,div=0,
#'  T=1,dt=0.04,dim=2, sim.func=sim.gbm, payoff.func=put.payoff,pilot.nsims=1000,
#'  batch.nrep=100,kernel.family="matern5_2",N=400)
#  kmSolve <- osp.fixed.design(model2d, input.dom=0.02,method="trainkm")
osp.fixed.design <- function(model,input.domain=NULL, method ="km",inTheMoney.thresh = 0, stop.freq=model$dt)
{
  M <- as.integer(round(model$T/model$dt))
  t.start <- Sys.time()
  
  if ( is.null(model$stop.times))
    model$stop.times = 1:M

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
 

  #----- step back in time
  design.size <- rep(0,M)

  for (i in (M-1):1) {
    if (i %in% model$stop.times == FALSE)
       next   # no stopping right now   
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
      init.grid <- init.grid[ model$payoff.func(grids[[i]], model) > inTheMoney.thresh,,drop=F]

      init.grid <- init.grid[sample(1:min(design.size[i],dim(init.grid)[1]), design.size[i], rep=F),,drop=F]
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
      
      if (is.null(model$min.lengthscale) & model$dim == 1 )
        model$min.lengthscale <- (my.domain[2]-my.domain[1])/100
      if (is.null(model$min.lengthscale) & model$dim > 1 )
        model$min.lengthscale <- (my.domain[,2]-my.domain[,1])/100
      if (is.null(model$max.lengthscale) )
        model$max.lengthscale <- 100*model$min.lengthscale

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
    big.grid <- matrix(rep(t(init.grid), n.reps), ncol = ncol(init.grid), byrow = TRUE)

    fsim <- forward.sim.policy( big.grid, M-i,fits[i:M],model,offset=0)
    #fsim <- policy.payoff( big.grid,(M-i)*model$dt/stop.freq,fits[i:M],model,offset=0,path.dt=stop.freq,interp=t.interp)
    immPayoff <- model$payoff.func( init.grid, model)
    cur.sim <- cur.sim + fsim$nsims

    # pre-averaged mean/variance
    for (jj in 1:design.size[i]) {
      all.X[jj,model$dim+1] <- mean( fsim$payoff[ jj + seq(from=0,len=n.reps,by=design.size[i])]) - immPayoff[ jj]
      all.X[jj,model$dim+2] <- var( fsim$payoff[ jj + seq(from=0,len=n.reps,by=design.size[i])])
    }

    all.X[,1:model$dim] <- init.grid  # use the first dim+1 columns for the batched GP regression.

    # create the km object
    if (n.reps > 10 & method == "km")
      fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[,model$dim+1]),
                                   noise.var=all.X[,model$dim+2]/n.reps, control=list(trace=F),
                                   coef.trend=0,coef.cov=model$km.cov, coef.var=model$km.var, covtype=model$kernel.family)
    else if (method == "km")  # manually estimate the nugget for small batches
      fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[,model$dim+1]),
                                   control=list(trace=F), lower=model$min.lengthscale, upper=model$max.lengthscale,
                                   coef.trend=0, coef.cov=model$km.cov, coef.var=model$km.var,
                                   nugget.estim=TRUE, nugget=sqrt(mean(all.X[,model$dim+2])), covtype=model$kernel.family)
    else if (method =="trainkm")
      fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[,model$dim+1]),
                                   control=list(trace=F), lower=model$min.lengthscale, upper=model$max.lengthscale,
                                   noise.var=all.X[,model$dim+2]/n.reps, covtype=model$kernel.family)

    else if (method == "mlegp")  { # laGP library implementation of regular GP
      fits[[i]]  <- laGP::newGP(X=init.grid, Z=all.X[,model$dim+1],
                                d=model$lagp.d, g=1e-6,dK=TRUE)
    
     laGP::jmleGP(fits[[i]], drange=c(model$min.lengthscale,model$max.lengthscale), grange=c(1e-8, 0.001))
    }
    else if(method =="hetgp") {
      big.payoff <- model$payoff.func(big.grid,model)
      hetData <- hetGP::find_reps(big.grid, fsim$payoff-big.payoff)
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
    if (method == "lagp") { # approximate GP
      if (is.null(model$lagp.type))
        lagp.type <- "alcray"
      else
        lagp.type <- model$lagp.type
      if (is.null(model$lagp.end))
        lagp.end <- 50
      else
        lagp.end <- model$lagp.end
      
      big.payoff <- model$payoff.func(big.grid,model)
      hetData <- hetGP::find_reps(big.grid, fsim$payoff-big.payoff)
      # normalize before running laGP
      scale_list <- find_shift_scale_stats(hetData$X0, hetData$Z0)
      # initial guess for lengthscale; do not estimate the nugget
      startd <- darg(NULL, shift_scale(list(X=init.grid), 
                                  shift=scale_list$xshift, scale=scale_list$xscale)$X)
      startd$start <- 1
     
      Ysc <- (fsim$payoff-big.payoff-scale_list$yshift)/scale_list$yscale
      Xsc <- shift_scale(list(X=big.grid), 
                           shift=scale_list$xshift, scale=scale_list$xscale)$X
      fits[[i]] <- list(Xsc=Xsc, Ysc=Ysc, scale=scale_list, d=startd)
      class(fits[[i]]) <- "agp"
    }
    
    else if (method == "ligp") {  # local inducing point Gaussian Process via the ligp library
      big.payoff <- model$payoff.func(big.grid,model)
      hetData <- hetGP::find_reps(big.grid, fsim$payoff-big.payoff)
      scale_list <- find_shift_scale_stats(hetData$X0, hetData$Z0)
      Xsc <- hetData$X0
      for (j in 1:model$dim)
        Xsc[,j] <- (Xsc[,j]-scale_list$xshift[j])/scale_list$xscale[j]
      
      
      Ysc <- (fsim$payoff-big.payoff-scale_list$yshift)/scale_list$yscale
      lhs_design <- randomLHS(model$ligp.lhs,model$dim)
      
      Xmt <- scale_ipTemplate(Xsc, model$ligp.n, 
                              space_fill_design=lhs_design, method='qnorm')$Xm.t
      
      Xsc <- shift_scale(list(X=big.grid), 
                         shift=scale_list$xshift, scale=scale_list$xscale)$X
      #Xsc <- big.grid
      #for (j in 1:model$dim)
      #  Xsc[,j] <- (Xsc[,j]-scale_list$xshift[j])/scale_list$xscale[j]
      hetData2 <- hetGP::find_reps(Xsc, Ysc)
      
      fits[[i]] <-  list(scale=scale_list, reps=hetData2, Xmt=Xmt)
      class(fits[[i]]) <- "ligprep"
    }
    else if (model$dim == 1 & method=="spline")  # only possible in 1D
      fits[[i]] <- stats::smooth.spline(x=init.grid,y=all.X[,2],nknots=model$nk)
    else if (method == "rvm") {
        if (is.null(model$rvm.kernel))
          rvmk <- "rbfdot"
        else
          rvmk <- model$rvm.kernel
        fits[[i]] <- kernlab::rvm(x=init.grid, y=all.X[,model$dim+1],kernel=rvmk)
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

        fits[[i]] <- np::npregbw(xdat = init.grid, ydat = all.X[,model$dim+1],
                           regtype=regtype, ckertype=kertype,ckerorder=kerorder)
    }

  }  # end of loop over time-steps

  return (list(fit=fits,timeElapsed=Sys.time()-t.start,nsims=cur.sim))
}


####################################
#' Swing option solver based on a batched non-adaptive design with a variety of regression methods
#'
#' @title Generic dynamic emulation of a multiple-stopping problem with a non-sequential design
#' @param model a list defining all the model parameters
#' @param input.domain the domain of the emulator. Several options are available. Default in \code{NULL}
#' All the empirical domains rely on pilot paths generated using \code{pilot.nsims}>0 model parameter.
#' \itemize{
#' \item NULL will use an empirical probabilistic design based on the pilot paths (default);
#' \item if a vector of length 2*model$dim then specifies the bounding rectangle
#' \item a single positive number, then build a bounding rectangle based on the \eqn{\alpha}-quantile of the pilot paths
#' \item a single negative number, then build a bounding rectangle based on the full range of the pilot paths
#' \item a vector specifies the precise design, used as-is (\emph{overrides design size})
#' }
#' @param method regression method to use (defaults to \code{km})
#' \itemize{
#' \item km: [(default] Gaussian process with fixed hyperparams  uses \pkg{DiceKriging} 
#' via \code{km}. Requires \code{km.cov} (vector of lengthscales)
#' and \code{km.var} (scalar process variance)  
#' \item trainkm: GP w/trained hyperparams: use \pkg{DiceKriging} via \code{km}. 
#' Requires to specify kernel family via \code{kernel.family}
#' \item mlegp Local GP from \pkg{laGP} (uses Gaussian squared exponential kernel)
#' \item homgp Homoskedastic GP: use \pkg{hetGP} with  \code{mleHomGP}. 
#' Requires to specify kernel family via \code{kernel.family}
#' \item hetgp Heteroskedastic GP: use \pkg{hetGP} with \code{mleHetGP}
#' Requires to specify kernel family via \code{kernel.family}
#' \item spline: Smoothing Splines, use \code{smooth.spline} from \pkg{base}  
#' with the user-specified \code{nk} number of knots (1D only)
#' \item cvspline: \code{smooth.spline} from \pkg{base} with automatically chosen 
#'  (via cross-validation) degrees of freedom/number of knots. Only works \emph{in 1D}
#' \item loess: Local polynomial regression: use \code{loess} with \code{lo.span} parameter
#' \item rvm: Relevance Vector Machine: use \pkg{kernlab} with \code{rvm}
#' \item lm: linear model from \pkg{stats} using \code{model$bases}
#' }
#' @param inTheMoney.thresh which paths are kept, out-of-the-money is dropped.
#' Defines threshold in terms of \code{model$payoff.func}
#' @return a list containing:
#' \itemize{
#' \item \code{fit} a list containing all the models generated at each time-step. \code{fit[[1]]} is the emulator
#' at \eqn{t=\Delta t}, the last one is \code{fit[[M-1]]} which is emulator for \eqn{T-\Delta t}.
#' \item \code{val}: the in-sample pathwise rewards
#' \item \code{test}: the out-of-sample pathwise rewards
#' \item \code{p}: the final price (2-vector for in/out-of-sample)
#' \item \code{timeElapsed} (based on \code{Sys.time})
#' }
#' @details Solves for a swing with \code{n.swing} exercise rights. The payoff function is 
#' saved in \code{swing.payoff}. Also assumes a refraction period of \code{refract} between consecutive 
#' exercises. The experimental design is based on  \code{\link{osp.fixed.design}}. By default, no forward evaluation is provided, ie the
#' method only builds the emulators. Thus, to obtain an actual estimate of the value
#' combine with \code{\link{swing.policy}}.
#' @export
#' 
#' @examples
#' set.seed(1)
#' swingModel <- list(dim=1, sim.func=sim.gbm, x0=100,
#' swing.payoff=put.payoff, n.swing=3,K=100, 
#' sigma=0.3, r=0.05, div=0,
#' T=1,dt=0.02,refract=0.1,
#' N=800,pilot.nsims=1000,batch.nrep=25)
#' swingModel$nk=16  # number of knots for the smoothing spline
#' spl.swing <- swing.fixed.design(swingModel,input.domain=0.03, method ="spline")
swing.fixed.design <- function(model,input.domain=NULL, method ="km",inTheMoney.thresh = 0)
{
  M <- as.integer(round(model$T/model$dt))
  t.start <- Sys.time()
  refractN <- model$refract/model$dt
  
  fits <- matrix(rep(list(),M*model$n.swing), nrow=M, ncol=model$n.swing)   # matrix of fits at each step
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
   for (kk in 1:model$n.swing)  {
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
      
      init.grid <- init.grid[sample(1:min(design.size[i],dim(init.grid)[1]), design.size[i], rep=F),,drop=F]
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
    
    init.grid <- init.grid[ model$swing.payoff(init.grid, model) > inTheMoney.thresh,,drop=F]
    
    design.size[i] <- dim(init.grid)[1]
    all.X <- matrix( rep(0, (model$dim+2)*design.size[i]), ncol=model$dim+2)
    
    # construct replicated design
    big.grid <- matrix(rep(t(init.grid), n.reps), ncol = ncol(init.grid), byrow = TRUE)
    
    fsim <- swing.policy( big.grid, M-i,fits[i:M,],model,offset=0,n.swing=kk)
    #fsim <- policy.payoff( big.grid,(M-i)*model$dt/stop.freq,fits[i:M],model,offset=0,path.dt=stop.freq,interp=t.interp)
    # if use one now, have kk-1 left (possibly zero is ok)
    immPayoff <- model$swing.payoff( big.grid, model) 
    #if (i < 30)
    #   browser()
    if (i+refractN < M & kk > 1) {
        delayedPayoff <- swing.policy(big.grid, M-i-refractN, fits[(i+refractN):M,],model,offset=0,n.swing=kk-1)
        immPayoff <- immPayoff + exp(-model$r*model$dt*refractN)*delayedPayoff$totPayoff
        cur.sim <- cur.sim + delayedPayoff$nsims
    }
    qValue = fsim$totPayoff - immPayoff
    cur.sim <- cur.sim + fsim$nsims #+ delayedPayoff$nsims
    
    # pre-averaged mean/variance
    for (jj in 1:design.size[i]) {
      all.X[jj,model$dim+1] <- mean( qValue[ jj + seq(from=0,len=n.reps,by=design.size[i])])
      all.X[jj,model$dim+2] <- var( qValue[ jj + seq(from=0,len=n.reps,by=design.size[i])])
    }
    
    all.X[,1:model$dim] <- init.grid  # use the first dim+1 columns for the batched GP regression.
    
    # create the km object
    if (method == "km")
      fits[[i,kk]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[,model$dim+1]),
                                   noise.var=all.X[,model$dim+2]/n.reps, control=list(trace=F),
                                   coef.trend=0,coef.cov=model$km.cov, coef.var=model$km.var, covtype=model$kernel.family)
    else if (method =="trainkm")
      fits[[i,kk]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[,model$dim+1]),
                                   control=list(trace=F), lower=model$min.lengthscale, upper=model$max.lengthscale,
                                   noise.var=all.X[,model$dim+2]/n.reps, covtype=model$kernel.family)
    
    else if (n.reps < 10 & method == "mlegp")  # laGP library implementation
      fits[[i,kk]]  <- laGP::newGP(X=init.grid, Z=all.X[,model$dim+1],
                                d=list(mle=FALSE, start=model$km.cov), g=list(start=1, mle=TRUE))
    else if(method =="hetgp") {
      hetData <- hetGP::find_reps(big.grid, qValue)
      fits[[i,kk]] <- hetGP::mleHetGP(X = list(X0=hetData$X0, Z0=hetData$Z0,mult=hetData$mult), Z= hetData$Z,
                                   lower = model$min.lengthscale, upper = model$max.lengthscale, 
                                   noiseControl=list(g_bounds=c(1e-4,100)), covtype=model$kernel.family)
      #ehtPred <- predict(x=check.x, object=hetModel)
    }
    else if (method =="homgp") {
      hetData <- hetGP::find_reps(big.grid, qValue)
      fits[[i,kk]] <- hetGP::mleHomGP(X = list(X0=hetData$X0, Z0=hetData$Z0,mult=hetData$mult), Z= hetData$Z,
                                   lower = model$min.lengthscale, upper = model$max.lengthscale, 
                                   covtype=model$kernel.family)
    }
    else if (method =="lm") {  # Linear regression with specified basis functions
      matb <- model$bases(init.grid)
      fits[[i,kk]] <- stats::lm(all.X[,model$dim+1] ~ matb)
    }
    else if (model$dim == 1 & method=="spline")  # only possible in 1D
      fits[[i,kk]] <- smooth.spline(x=init.grid,y=all.X[,2],nknots=model$nk)
    else if (model$dim == 1 & method=="cvspline")  # only possible in 1D
      fits[[i,kk]] <- smooth.spline(x=init.grid,y=all.X[,2])  # cross-validated DF
    else if (method == "rvm") {
      if (is.null(model$rvm.kernel))
        rvmk <- "rbfdot"
      else
        rvmk <- model$rvm.kernel
      fits[[i,kk]] <- kernlab::rvm(x=init.grid, y=all.X[,model$dim+1],kernel=rvmk)
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
      
      fits[[i,kk]] <- np::npregbw(xdat = init.grid, ydat = all.X[,model$dim+1],
                         regtype=regtype, ckertype=kertype,ckerorder=kerorder)
    }
  }  # end of loop over swing rights kk  
  }  # end of loop over time-steps
  
  return (list(fit=fits,timeElapsed=Sys.time()-t.start,nsims=cur.sim))
}

#############################
#' RMC using TvR along a global set of paths.
#' All designs are kept in memory
#' @title Tsitsiklis van Roy RMC algorithm with a variety of regression methods
#' @param model a list defining the simulator and reward model, with the two main model hooks being 
#' \code{payoff.func} (plus parameters) and \code{sim.func} (plus parameters). 
#' 
#' Also \code{x0}
#' is a required part of the \code{model}. Can be either a vector of length \code{model$dim} 
#' or a vector of length \code{model$dim}*N
#' @param N the number of forward paths to train on
#' @param subset To reserve out-of-sample paths, specify \code{subset} (eg 1:1000) to use for testing.
#' By default everything is in-sample.
#' 
#' @param method a string specifying regression method to use
#' \itemize{
#'  \item spline: \code{smooth.spline} from \pkg{base} which only works \emph{in 1D}
#'  \item cvspline: \code{smooth.spline} from \pkg{base} with automatically chosen 
#'  (via cross-validation) degrees of freedom/number of knots. Only works \emph{in 1D}
#'  \item randomforest: (from \pkg{randomForest} package) requires \code{rf.maxnode}
#'  and \code{rf.ntree} (number of trees) model parameters
#'  \item loess: only works in \emph{1D or 2D}, requires \code{lo.span} model parameter
#'  \item earth: multivariate regression splines (MARS) using \pkg{earth} package.
#'  requires \code{earth.deg} (interaction degree), \code{earth.nk} (max number of terms to keep),
#'  \code{earth.thresh} params
#'  \item rvm: relevance vector machine from \pkg{kernlab} package. Optional \code{rvm.kernel}
#'  model parameter to decide which kernel family to utilize. Default kernel is rbfdot
#'  \item deepnet: neural network using \pkg{deepnet}. Specify \code{nn.layers} as a vector 
#'  to describe the number of nodes across hidden layers
#'  \item lm [Default]: linear global regression using \code{model$bases} (required) basis functions (+ constant)
#'  }
#' @export
#' @return a list containing
#' \itemize{
#' \item \code{fit} a list containing all the models generated at each time-step. \code{fit[[1]]} is the emulator
#' at \eqn{t=\Delta t}, the last one is \code{fit[[M-1]]} which is emulator for \eqn{T-\Delta t}.
#' \item \code{val}: the in-sample pathwise rewards
#' \item \code{test}: the out-of-sample pathwise rewards
#' \item \code{p}: the final price (2-vector for in/out-of-sample)
#' \item \code{timeElapsed} (based on \code{Sys.time})
#' }
#' @details
#'  Works with a probabilistic design that requires storing all paths in memory. Specifying \code{subset}
#'  allows to compute in parallel with the original computation an out-of-sample estimate of the value function
#'  
#'  Calls \code{model$payoff.func}, so the latter must be set prior to calling.
#'  Also needs \code{model$dt} and \code{model$r} for discounting
#'  
#'  Calls \code{model$sim.func} to generate forward paths
#'  
#'  Emulator is trained on all paths, even those that are out-of-the-money
#'  
#' @examples
#' set.seed(1)
#' require(earth)
#' model2d <- list(K=40,x0=rep(40,2),sigma=rep(0.2,2),r=0.06,div=0,
#'  T=1,dt=0.04,dim=2, sim.func=sim.gbm, payoff.func=put.payoff,pilot.nsims=1000,
#'  earth.deg=2,earth.nk=200,earth.thresh=1E-8)
#' tvrSolve <- osp.tvr(N=41000,model2d, subset=1:1000,method="earth")
#' # "in-sample v_0 1.224009; and out-of-sample: 1.233986"
###############################
osp.tvr <- function(N,model,subset=1:N,method="lm")
{
  M <- model$T/model$dt
  if ( is.null(model$stop.times))
    model$stop.times = 1:M
  grids <- list()
  all.models <- list()
  if (is.null(subset))
     subset = 1:N
  
  # divide into in-sample and out-of-sample
  if (length(subset) < N)
    train <- (1:N)[-subset]
  else
    train <- 1:N
  
  # Build the designs based on a simulation of X_{1:T}
  if (length(model$x0) == model$dim)
    grids[[1]] <- model$sim.func( matrix(rep(model$x0, N), nrow=N,byrow=T), model, model$dt)
  else if (length(model$x0) == N*model$dim)
    grids[[1]] <- model$sim.func( matrix(rep(model$x0, N), nrow=N,byrow=T), model, model$dt)
  else
    warning("Length of the initial condition must match N")
  for (i in 2:M)
    grids[[i]] <- model$sim.func( grids[[i-1]], model, model$dt)
  
  # initialize at T
  contValue <- exp(-model$r*model$dt)*model$payoff.func( grids[[M]], model)
  tau <- rep(model$T, N)
  t.start <- Sys.time()
  
  # Backward stepping in time
  # Estimate T(t,x)
  for (i in (M-1):1)
  {
    if (i %in% model$stop.times == FALSE)
      next   # no stopping right now 
    # forward predict
    immPayoff <- model$payoff.func(grids[[i]],model)
    
    yVal <- contValue-immPayoff
    
    ### CASES DEPENDING ON METHOD
    if (method == "spline" & ncol(grids[[i]]) == 1) { # only works in 1D
      all.models[[i]] <- stats::smooth.spline( x=grids[[i]],y=yVal,
                                        nknots = model$nk)
      timingValue <- predict(all.models[[i]],grids[[i]])$y
    }
    if (method == "cvspline" & ncol(grids[[i]]) == 1) { # only works in 1D
      all.models[[i]] <- stats::smooth.spline( x=grids[[i]],y=yVal)  # cross-validate DF
      timingValue <- predict(all.models[[i]],grids[[i]])$y
    }
    
    if (method == "randomforest") {
      if (is.null(model$rf.ntree) | is.null(model$rf.maxnode))
        stop("Missing parameters to pass to randomForest function. Need rf.maxnode and rf.ntree")
      
      all.models[[i]] <-  randomForest::randomForest(x=grids[[i]],y=yVal,
                                       ntree=model$rf.ntree,replace=F,maxnode=model$rf.maxnode)
      timingValue <- predict(all.models[[i]],grids[[i]],predict.all=T)$individual
      timingValue <- apply(timingValue,1,mean)   # median or mean prediction could be used
    }
    
    if (method == "ligp") {  # local inducing point Gaussian Process via the ligp library
      scale_list <- find_shift_scale_stats(as.matrix(grids[[i]]), yVal)
      Xsc <- shift_scale(list(X=as.matrix(grids[[i]])), 
                             shift=scale_list$xshift, scale=scale_list$xscale)$X

      Ysc <- (yVal-scale_list$yshift)/scale_list$yscale
      lhs_design <- randomLHS(model$ligp.lhs,model$dim)

      Xmt <- scale_ipTemplate(Xsc, model$ligp.n, space_fill_design=lhs_design, method='qnorm')$Xm.t

      out <- liGP(XX=Xsc, X=Xsc, Y=Ysc, Xm=Xmt, N=model$ligp.n, theta=model$ligp.theta,
                  g=list(start=1e-4,min=1e-6, max=model$ligp.gmax))
      timingValue <- out$mean*scale_list$yscale + scale_list$yshift

      all.models[[i]] <-  list(scale=scale_list, Xsc=Xsc, Ysc = Ysc, Xmt=Xmt)
      class(all.models[[i]]) <- "ligp"
    }
    
    if (method == "loess" & ncol(grids[[i]]) <= 2) { # LOESS only works in 1D or 2D
      if (ncol(grids[[i]]) == 1) {
        all.models[[i]] <- stats::loess(y ~ x, data.frame(x=grids[[i]], y=yVal),
                                 span=model$lo.span, control = loess.control(surface = "direct"))
        timingValue <- predict(all.models[[i]], data.frame(x=grids[[i]]),se=TRUE)$fit
      }
      if (ncol(grids[[i]]) ==2) {
        all.models[[i]] <- stats::loess(y ~ x1+x2, data.frame(x1=grids[[i]][,1], x2=grids[[i]][,2], y=yVal),
                                 span=model$lo.span, control = loess.control(surface = "direct"))
        timingValue <- predict(all.models[[i]], new=data.frame(x1=grids[[i]][,1],x2=grids[[i]][,2]))$fit
      }
    }
    if (method == "earth") {  # Multivariate Adaptive Regression Splines
      if (is.null(model$earth.deg) | is.null(model$earth.nk) | is.null(model$earth.thresh))
        stop("Missing parameters to pass to earth::earth function. Need earth.deg, earth.nk, earth.thresh model fields")
      
      all.models[[i]] <- earth::earth(x=grids[[i]],y=yVal,
                               degree=model$earth.deg,nk=model$earth.nk,thresh=model$earth.thresh)
      timingValue <- predict(all.models[[i]],grids[[i]])
    }
    if (method == "deepnet") {  # Neural Network via the deepnet library
      if (is.null(model$nn.layers) )
        stop("Missing parameters to pass to deepnet::nn.train function. Need model$nn.layers")
      
      all.models[[i]] <- deepnet::nn.train(x=grids[[i]][c.train,,drop=F],y=yVal,
                                  hidden=model$nn.layers)
      timingValue <- deepnet::nn.predict(all.models[[i]],grids[[i]])
    }
    if (method == "nnet") {  # Neural Network via the nnet library
      if (is.null(model$nn.nodes) )
        stop("Missing parameters to pass to nnet::nnet function. Need model$nn.nodes")
      
      all.models[[i]] <- nnet::nnet(x=grids[[i]][c.train,,drop=F],y=yVal,
                                  size=model$nn.nodes, linout=TRUE, maxit=1000,trace=FALSE)
      timingValue <- predict(all.models[[i]],grids[[i]], type="raw")
    }
    
    # Default
    if (method =="lm") {  # Linear regression with specified basis functions
      matb <- model$bases(grids[[i]])
      all.models[[i]] <- stats::lm(yVal ~ matb)
      lenn <- length(all.models[[i]]$coefficients)
      timingValue <-  all.models[[i]]$coefficients[1] +
        model$bases(grids[[i]]) %*% all.models[[i]]$coefficients[2:lenn]
    }
    if (method =="rvm") {   # Relevance Vector Machine Regression
      if (is.null(model$rvm.kernel))
        rvmk <- "rbfdot"
      else
        rvmk <- model$rvm.kernel
      rvModel <- kernlab::rvm(x=grids[[i]], y=yVal,kernel=rvmk)
      timingValue <- predict(rvModel, new=grids[[i]])
    }
    
    # paths on which stop right now
    stopNdx <- which( timingValue <= 0 & immPayoff > 0)
    contValue[stopNdx] <- immPayoff[stopNdx]
    tau[stopNdx] <- i*model$dt
    
    # else continue. plug-in predicted timingValue and discount
    contValue <- exp(-model$r*model$dt)*(timingValue + immPayoff)
  }
  
  # in/sample and out-of-sample average at x0
  test <- NULL
  if (length(subset) < N & length(subset) > 0) {
    price <- c(mean(contValue[train]),mean(contValue[subset]))
    print(sprintf("in-sample v_0 %3f; and out-of-sample: %3f", price[1], price[2]))
    test <- contValue[subset]
  }
  else
    price <- mean(contValue)
  
  # returns a list containing
  # fit are all the models generated at each time-step, stored as a list
  # p is the final price (2-vector for in/out-of-sample)
  # val are the in-sample pathwise rewards
  # test are the out-of-sample pathwise rewards
  # timeElapsed: total running time
  return( list(fit=all.models,p=price, val=contValue[train], test=test,
               timeElapsed=Sys.time()-t.start))
}


#############################
#' RMC for impulse control.
#' Training design specified explicitly by the user
#' @title LS-flavor RMC algorithm with a variety of regression methods for stochastic impulse control
#' @param model a list defining the simulator and reward model, with the two main model hooks being 
#' \code{impulse.func} (plus parameters) and \code{sim.func} (plus parameters). 
#' 
#' 
#' @param method a string specifying regression method to use
#' \itemize{
#'  \item spline [Default]: \code{smooth.spline} from \pkg{base} which only works \emph{in 1D}
#'  \item randomforest: (from \pkg{randomForest} package) requires \code{rf.maxnode}
#'  and \code{rf.ntree} (number of trees) model parameters
#'  \item loess: only works in \emph{1D or 2D}, requires \code{lo.span} model parameter
#'  \item deepnet: neural network using \pkg{deepnet}. Specify \code{nn.layers} as a vector 
#'  to describe the number of nodes across hidden layers
#'  \item homgp Homoskedastic GP: use \pkg{hetGP} with  \code{mleHomGP}
#' \item hetgp Heteroskedastic GP: use \pkg{hetGP} with \code{mleHetGP}
#'  \item lm: linear global regression using \code{model$bases} (required) basis functions (+ constant)
#'  }
#' @export
#' @return a list containing
#' \itemize{
#' \item \code{fit} a list containing all the models generated at each time-step. \code{fit[[1]]} is the emulator
#' at \eqn{t=\Delta t}, the last one is \code{fit[[M-1]]} which is emulator for \eqn{T-\Delta t}.
#' \item \code{timeElapsed} (based on \code{Sys.time})
#' }
#' @details
#'  Works with a design specified by the user
#'  
#'  Calls \code{model$impulse.func}, so the latter must be set prior to calling.
#'  Also needs \code{model$dt} and \code{model$r} for discounting. 
#'  
#'  Calls \code{model$sim.func} to generate forward paths. Use in conjunction with
#'  \code{\link{forward.impulse.policy}} 
#'  
#' @author 
#' Mike Ludkovski
#'  
#' @examples
#' set.seed(1)
#' require(DiceKriging)
#' modelBelak <- list(dim=1, sim.func=sim.bm, r=0.5, drift=0, sigma=1, 
#' x0=1,  impulse.fixed.cost = 1,impulse.target = 0,impulse.func = forest.impulse,
#' imp.type = "forest",T=5, dt=0.05,pilot.nsims=0,batch.nrep = 10,nk = 30,N = 601)
#' belSolve <- osp.impulse.control(modelBelak, input.domain = seq(-0.5,2.5,by=0.005),method="spline")
###############################

osp.impulse.control <- function(model,input.domain=NULL, method ="spline",verb=101, mpc=FALSE)
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
  
  
  #----- step back in time
  design.size <- rep(0,M)
  fits[[M]] <- NULL
  
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
      #init.grid <- init.grid[ model$payoff.func(grids[[i]], model) > inTheMoney.thresh,,drop=F]
      
      init.grid <- init.grid[sample(1:min(design.size[i],dim(init.grid)[1]), design.size[i], rep=F),,drop=F]
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
      
      if (is.null(model$min.lengthscale) & model$dim == 1 )
        model$min.lengthscale <- (my.domain[2]-my.domain[1])/100
      if (is.null(model$min.lengthscale) & model$dim > 1 )
        model$min.lengthscale <- (my.domain[,2]-my.domain[,1])/100
      if (is.null(model$max.lengthscale) )
        model$max.lengthscale <- 100*model$min.lengthscale
      
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
    else if (length(input.domain) == M){
      init.grid <- matrix(input.domain[[M]],nrow=length(input.domain[[M]])/model$dim)
      design.size[i] <- nrow(init.grid)
      
    }
    else  {   # fixed pre-specified design
      init.grid <- matrix(input.domain,nrow=length(input.domain)/model$dim)
      design.size[i] <- nrow(init.grid)
    }
    
    #init.grid <- init.grid[ model$payoff.func(init.grid, model) > inTheMoney.thresh,,drop=F]
    
    design.size[i] <- dim(init.grid)[1]
    all.X <- matrix( rep(0, (model$dim+2)*design.size[i]), ncol=model$dim+2)
    
    # construct replicated design
    big.grid <- matrix(rep(t(init.grid), n.reps), ncol = ncol(init.grid), byrow = TRUE)
    
    fsim <- forward.impulse.policy( big.grid, M-i,fits[(i+1):M],model, mpc=mpc)
    fsim <- pmax( fsim$payoff, 0)
    
    #if (i < (M-1))
    #   immPayoff <- model$impulse.func( init.grid, model, fits[[i+1]])
    #else
    #   immPayoff <- rep(0, design.size[i])
    
    # pre-averaged mean/variance
    for (jj in 1:design.size[i]) {
      all.X[jj,model$dim+1] <- mean( fsim[ jj + seq(from=0,len=n.reps,by=design.size[i])]) #- immPayoff[ jj]
      all.X[jj,model$dim+2] <- var( fsim[ jj + seq(from=0,len=n.reps,by=design.size[i])])
    }
    
    all.X[,1:model$dim] <- init.grid  # use the first dim+1 columns for the batched GP regression.
    if (i %% verb == 0)
      browser()
    
    # create the km object
    if (n.reps > 10 & method == "km")
      fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[,model$dim+1]),
                                   noise.var=all.X[,model$dim+2]/n.reps, control=list(trace=F),
                                   coef.trend=0,coef.cov=model$km.cov, coef.var=model$km.var, covtype=model$kernel.family)
    else if (method == "km")  # manually estimate the nugget for small batches
      fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[,model$dim+1]),
                                   control=list(trace=F), lower=model$min.lengthscale, upper=model$max.lengthscale,
                                   coef.trend=0, coef.cov=model$km.cov, coef.var=model$km.var,
                                   nugget.estim=TRUE, nugget=sqrt(mean(all.X[,model$dim+2])), covtype=model$kernel.family)
    else if (method =="trainkm")
      fits[[i]] <- DiceKriging::km(y~0, design=data.frame(x=init.grid), response=data.frame(y=all.X[,model$dim+1]),
                                   control=list(trace=F), lower=model$min.lengthscale, upper=model$max.lengthscale,
                                   noise.var=pmax(1e-6,all.X[,model$dim+2]/n.reps), covtype=model$kernel.family)
    
    else if (method == "mlegp")  { # laGP library implementation of regular GP
      fits[[i]]  <- laGP::newGP(X=init.grid, Z=all.X[,model$dim+1],
                                d=model$lagp.d, g=1e-6,dK=TRUE)
      
      laGP::jmleGP(fits[[i]], drange=c(model$min.lengthscale,model$max.lengthscale), grange=c(1e-8, 0.001))
    }
    else if(method =="hetgp") {
      big.payoff <- model$payoff.func(big.grid,model)
      hetData <- hetGP::find_reps(big.grid, fsim$payoff-big.payoff)
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
    else if (method == "ligp") {  # local inducing point Gaussian Process via the ligp library
      big.payoff <- model$payoff.func(big.grid,model)
      hetData <- hetGP::find_reps(big.grid, fsim$payoff-big.payoff)
      scale_list <- find_shift_scale_stats(hetData$X0, hetData$Z0)
      Xsc <- hetData$X0
      for (j in 1:model$dim)
        Xsc[,j] <- (Xsc[,j]-scale_list$xshift[j])/scale_list$xscale[j]
      
      
      Ysc <- (fsim$payoff-big.payoff-scale_list$yshift)/scale_list$yscale
      lhs_design <- randomLHS(model$ligp.lhs,model$dim)
      
      Xmt <- scale_ipTemplate(Xsc, model$ligp.n, 
                              space_fill_design=lhs_design, method='qnorm')$Xm.t
      
      Xsc <- shift_scale(list(X=big.grid), 
                         shift=scale_list$xshift, scale=scale_list$xscale)$X
      #Xsc <- big.grid
      #for (j in 1:model$dim)
      #  Xsc[,j] <- (Xsc[,j]-scale_list$xshift[j])/scale_list$xscale[j]
      hetData2 <- hetGP::find_reps(Xsc, Ysc)
      
      fits[[i]] <-  list(scale=scale_list, reps=hetData2, Xmt=Xmt)
      class(fits[[i]]) <- "ligprep"
    }
    else if (model$dim == 1 & method=="spline")  # only possible in 1D
      fits[[i]] <- stats::smooth.spline(x=init.grid,y=all.X[,2],nknots=model$nk)
    else if (model$dim == 1 & method=="cvspline")  # only possible in 1D
      fits[[i]] <- stats::smooth.spline(x=init.grid,y=all.X[,2])
    else if (method == "rvm") {
      if (is.null(model$rvm.kernel))
        rvmk <- "rbfdot"
      else
        rvmk <- model$rvm.kernel
      fits[[i]] <- kernlab::rvm(x=init.grid, y=all.X[,model$dim+1],kernel=rvmk)
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
      
      fits[[i]] <- np::npregbw(xdat = init.grid, ydat = all.X[,model$dim+1],
                             regtype=regtype, ckertype=kertype,ckerorder=kerorder)
    }
    
  }  # end of loop over time-steps
  
  return (list(fit=fits,timeElapsed=Sys.time()-t.start))
}
