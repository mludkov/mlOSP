############################################
#' Simulate h(X_tau) using FIT (can be a dynaTree, smooth.spline, MARS or RandomForest or RVM or hetGP)
#' @title Forward simulation based on a sequence of emulators
#' @param model a list containing all model parameters
#' @param offset (internal for debugging purposes)
#' @param x     is a matrix of starting values
#'
#' if input \code{x} is a list, then use the grids specified by x
#' @param M     number of time steps to forward simulate
#' @param fit   a list of fitted emulators that determine the stopping classifiers to be used
#' @param compact flag; if FALSE returns additional information about forward x-values.
#' @param use.qv boolean to indicate whether to plug-in continuation value for allpaths 
#' still alive at the last time-step. Default is set to \code{FALSE}
#' @export
#' @return a list containing:
#' \itemize{
#'  \item \code{payoff} is the resulting payoff NPV from t=0
#'  \item \code{fvalue[i]} is a list of resulting payoffs (on paths still not stopped) NPV from t=i
#'  \item \code{tau} are the times when stopped
#'  \item \code{sims} is a list; \code{sims[[i]]} are the forward x-values of paths at t=i (those not stopped yet)
#' \code{nsims} number of total 1-step simulations performed
#' }
#' @details Should be used in conjuction with the \code{osp.xx} functions that build the emulators. Also called
#' internally from \code{\link{osp.fixed.design}}
forward.sim.policy <- function( x,M,fit,model,offset=1,compact=TRUE,use.qv=FALSE)
{
  if (is.null(model$look.ahead)) 
      model$look.ahead <- 1
  
  nsim <- 0
  if (is(x, "matrix") | is(x, "numeric") ) {
    curX <- model$sim.func( x,model,model$dt)
    nsim <- nrow(x)
  }
  if (is(x, "list" ) )
    curX <- x[[1]]
  payoff <- rep(0, nrow(curX))
  tau <- rep(0, nrow(curX))
  sims <- vector(mode="list",len=model$look.ahead+1)
  save.ndx <- vector(mode="list",len=model$look.ahead+1)
  fvalue <- vector(mode="list",len=model$look.ahead+1)
  
  contNdx <- 1:nrow(curX)
  i <- 1
  payoff[contNdx]  <- exp(-model$dt*model$r)*model$payoff.func( curX[contNdx,,drop=F], model)
  
  # main loop forward
  while (i < (M+(use.qv==TRUE)) & length(contNdx) > 0) {
    
    if (use.qv == TRUE & i== M) 
      in.the.money <- 1:length(contNdx)  # predict everywhere since need q(t,x)
    else
     in.the.money <-  which( payoff[contNdx] > 0) 
    if (length(in.the.money) >0 & is.null( fit[[i+1-offset]]) == FALSE) {
      
      if (is(fit[[i+1-offset]],"earth") )
        rule <- predict(fit[[i+1-offset]],curX[contNdx[in.the.money],,drop=F]) # for use with  MARS
      if (is(fit[[i+1-offset]],"deepnet") )
        rule <- nn.predict(fit[[i+1-offset]],curX[contNdx[in.the.money],,drop=F]) # for use with deepnet
      if (is(fit[[i+1-offset]],"nnet") )
        rule <- predict(fit[[i+1-offset]],curX[contNdx[in.the.money],,drop=F], type="raw") # for use with nnet
      if (is(fit[[i+1-offset]],"smooth.spline") )
        rule <- predict(fit[[i+1-offset]],curX[contNdx[in.the.money],,drop=F])$y # for use with  splines
      if (is(fit[[i+1-offset]],"randomForest") ) {
        obj <- predict(fit[[i+1-offset]],curX[contNdx[in.the.money],,drop=F],predict.all=T)$individual
        rule <- apply(obj,1,median)
      }
      if (is(fit[[i+1-offset]],"dynaTree")  ){
        obj <- predict(fit[[i+1-offset]], curX[contNdx[in.the.money],,drop=F],quants=F)
        rule <- obj$mean
      }
      if (is(fit[[i+1-offset]],"km") ) # DiceKriging
        rule <- predict(fit[[i+1-offset]],data.frame(x=curX[contNdx[in.the.money],,drop=F]),type="UK")$mean
      if (is(fit[[i+1-offset]],"integer") )  # laGP
        rule <- predGP(fit[[i+1-offset]],XX=curX[contNdx[in.the.money],,drop=F], lite=TRUE ,nonug=TRUE)$mean
      if (is(fit[[i+1-offset]],"lm") ) {
        lenn <- length(fit[[i+1-offset]]$coefficients)
        rule <-  fit[[i+1-offset]]$coefficients[1] + model$bases(curX[contNdx[in.the.money],,drop=F]) %*%
          fit[[i+1-offset]]$coefficients[2:lenn]
      }
      if( class(fit[[i+1-offset]])=="homGP" | class(fit[[i+1-offset]]) == "hetGP" | class(fit[[i+1-offset]]) == "homTP")
        rule <- predict(x=curX[contNdx[in.the.money],,drop=F], object=fit[[i+1-offset]])$mean
      if( class(fit[[i+1-offset]])=="rvm")
        rule <-  predict(fit[[i+1-offset]], new=curX[contNdx[in.the.money],,drop=F])
      if (class(fit[[i+1-offset]]) == "npregression")
        rule <- predict(fit[[i+1-offset]], new=curX[contNdx[in.the.money],,drop=F])
      if(class(fit[[i+1-offset]]) == "ligp") {
        newX <- shift_scale(list(X=as.matrix(curX[contNdx[in.the.money],,drop=F])), 
        #                   #list(sigma_vec = matrix(sigma_vec, nrow=1)),
                           shift=fit[[i+1-offset]]$scale$xshift, scale=fit[[i+1-offset]]$scale$xscale)$X
        
        #Xmt <- scale_ipTemplate(Xsc, model$ligp.n, space_fill_design=lhs_design, method='chr')$Xm.t
        out <- liGP(XX=newX, X=fit[[i+1-offset]]$Xsc, Y=fit[[i+1-offset]]$Ysc, 
                     Xm=fit[[i+1-offset]]$Xmt, N=model$ligp.n, theta=model$ligp.theta,
                    g=list(start=1e-4,min=1e-6, max=1e-2), epsQ=1e-4)
        rule <- out$mean*fit[[i+1-offset]]$scale$yscale + fit[[i+1-offset]]$scale$yshift
        
      }
      if(class(fit[[i+1-offset]]) == "ligprep") {
        newX <- shift_scale(list(X=as.matrix(curX[contNdx[in.the.money],,drop=F])), 
                            shift=fit[[i+1-offset]]$scale$xshift, scale=fit[[i+1-offset]]$scale$xscale)$X
        #d <- darg(NULL, fit[[i+1-offset]]$reps$X0)
        #d$start <- 0.1
        g <- garg(list(mle=TRUE), fit[[i+1-offset]]$reps$Z0)
        g$max <- model$ligp.gmax
        #ligp_qnorm <- liGP(XX, Xm.t=Xmt.qnorm$Xm.t, N=N,g=g, theta=d,
        #                   num_thread=num_thread, reps=reps_list)
        out <- liGP(XX=newX, Xm=fit[[i+1-offset]]$Xmt, N=model$ligp.n, theta=model$ligp.theta,
                    g=g, reps=fit[[i+1-offset]]$reps, num_thread=model$ligp.cores)
        rule <- out$mean*fit[[i+1-offset]]$scale$yscale + fit[[i+1-offset]]$scale$yshift
        
      }
      
      if (use.qv == TRUE & i== M) {
        payoff[contNdx] = payoff[contNdx] + rule  # continuation value of paths that didn't stop yet
        break
      }
      
      # stop if the expected gain is negative
      #contNdx <- contNdx[ which (rule > 0 | payoff[contNdx] == 0)]
      if (length(which(rule <0)) > 0)
        contNdx <- contNdx[-in.the.money[which(rule < 0)]  ]
    }
    tau[contNdx] <- (i)*model$dt
    
    if (compact == F) {
      sims[[min(i,model$look.ahead+1)]] <- curX[contNdx,,drop=F]
      save.ndx[[min(i,model$look.ahead+1)]] <- contNdx
    }
    # update the x values by taking a step of length model$dt for those paths not exercised yet
    i <- i+1
    
    if (length(contNdx) > 0) {
      if (is(x, "matrix") | is(x, "numeric") ) {
        curX[contNdx,] <- model$sim.func( curX[contNdx,,drop=F],model,model$dt )
        nsim <- nsim + length(contNdx)
      }
      
      if (is(x, "list") )  # stored list of paths
        curX[contNdx,] <- x[[i]][contNdx,,drop=F]
      # payoff for next timestep -- used for terminal payoff at i=M
      payoff[contNdx]  <- exp(-(i)*model$dt*model$r)*model$payoff.func( curX[contNdx,,drop=F], model)
    }
  }
  for (i in 2:(model$look.ahead))   # payoff for a trajectory starting at x^n_{t+i} which was still alive then
    fvalue[[i]] <- payoff[ save.ndx[[i]] ]*exp((i-1)*model$dt*model$r)
  
  
  return( list(payoff=payoff,sims=sims,fvalue=fvalue,tau=tau+model$dt, nsims=nsim))
  # payoff is the resulting payoff NPV from t=0
  # fvalue[i] is a list of resulting payoffs (on paths still not stopped) NPV from t=i
  # tau are the times when stopped -- Added model$dt since smallest is after 1 step
  # sims is a list; sims[[i]] are the forward x-values of paths at t=i (those not stopped yet)
}

########################################################################################

######
#' Two-dimensional raster+contour+point plot of an mlOSP emulator at a single time step.
#' @title Visualize a 2D emulator + stopping region
#' @details Uses the raster plot from \pkg{ggplot2}. For GP-based objects, also shows the unique design
#' sites via geom_point. See \code{\link{plt.2d.surf.batch}} for a similar plot 
#' for \code{osp.seq.batch.design} emulators.
#'
#' @param x,y locations to use for the \code{predict()} functions. Default is a 200x200 fine grid.
#' Passed to \code{expand.grid}
#' @param fit a fitted emulator. can be any of the types supported by \code{\link{forward.sim.policy}}
#' @param show.var if \code{TRUE} then plot posterior surrogate variance instead of surrogate mean [default = FALSE]
#' This only works for \code{km} and \code{het/homGP/homTP} objects
#' @param only.contour -- just the zero-contour, no raster plot
#' @param contour.col (default is "red") -- color of the zero contour
#' @param bases (only used for lm objects)
#' @return a ggplot handle for the created plot.
#' @export
plt.2d.surf <- function( fit, x=seq(31,43,len=201),y = seq(31,43,len=201),
                         show.var=FALSE, only.contour=FALSE, contour.col="red",
                         bases=NULL, strike=0)
{
  gr <- expand.grid(x=x,y=y)
  
  if (is(fit,"randomForest") ) {
    obj <- predict(fit,cbind(gr$x,gr$y),predict.all=T)$individual
    obj <- apply(obj,1,median)
  }
  if (is(fit,"earth") )
    obj <- predict(fit,cbind(gr$x,gr$y)) # for use with  MARS
  if (is(fit,"deepnet") )
    obj <- nn.predict(fit,cbind(gr$x,gr$y)) # for use with deepnet
  if (is(fit,"nnet") )
    obj <- predict(fit,cbind(gr$x,gr$y), type="raw") # for use with nnet
  
  if (is(fit,"dynaTree")  )
    obj <- predict(fit,cbind(gr$x,gr$y),quants=F)$mean
  if (is(fit,"lm") ) {
    obj <-  fit$coefficients[1] + bases(cbind(gr$x,gr$y)) %*%
      fit$coefficients[2:length(fit$coefficients)]
  }
  
  if (class(fit)=="km" & show.var == FALSE)
    obj <- predict(fit,data.frame(x=cbind(gr$x,gr$y)), type="UK")$mean
  if (is(fit,"km") & show.var == TRUE)
    obj <- predict(fit,data.frame(x=cbind(gr$x,gr$y)), type="UK")$sd
  if( (class(fit)=="homGP" | class(fit) == "hetGP" | class(fit) == "homTP") & show.var == FALSE)
    obj <- predict(x=cbind(gr$x,gr$y), object=fit)$mean
  if( (class(fit)=="homGP" | class(fit) == "hetGP" | class(fit) == "homTP") & show.var == TRUE)
    obj <- sqrt(predict(x=cbind(gr$x,gr$y), object=fit)$sd2)
  if( class(fit)=="rvm")
    obj <-  predict(fit, new=cbind(gr$x,gr$y))
  if (class(fit) == "npregression")
    obj <- predict(fit,exdat=cbind(gr$x,gr$y))
  if (class(fit)== "ligprep") {
    newX <- shift_scale(list(X=as.matrix(cbind(gr$x, gr$y))), 
                        shift=fit$scale$xshift, scale=fit$scale$xscale)$X
    #d <- darg(NULL, ft$reps$X0)
    #d$start <- 0.1
    g <- garg(list(mle=TRUE), fit$reps$Z0)
    g$max <- 1
    out <- liGP(XX=newX, Xm=fit$Xmt, N=100, theta=1,
                g=g, reps=fit$reps, num_thread=4)
    obj <- out$mean*fit$scale$yscale + fit$scale$yshift
  }
  # to cut-out the OTM region for max-Call
  if (strike > 0){
    payoffs <-  pmax(apply(gr,1,max)-strike,0)
    obj[ which(payoffs == 0)] <- NaN
  }
  
  if (only.contour==TRUE) {
    #imx <- as.image(x=cbind(gr$x,gr$y),obj,nr=100,nc=100)
    #contour(imx$x,imx$y,imx$z,levels=0,add=T,drawlab=F,lwd=2,col=contour.col)
    g1 <- ggplot( data.frame(x=gr$x, y=gr$y,z=obj)) + 
      geom_contour(breaks=0,color=contour.col,aes(x,y,z=z),size=1.4) +
      scale_x_continuous(expand=c(0,0),limits=range(x)) + 
      scale_y_continuous(expand=c(0,0),limits=range(y)) +
      labs(x=expression(X[t]^1),y=expression(X[t]^2))
  }
  else {
    g1 <- ggplot( data.frame(x=gr$x, y=gr$y,z=obj)) + 
      scale_x_continuous(expand=c(0,0),limits=range(x)) + 
      scale_y_continuous(expand=c(0,0),limits=range(y)) +
      labs(x=expression(X[t]^1),y=expression(X[t]^2)) +
      geom_raster(aes(x,y,fill=z)) + 
      geom_contour(breaks=0,color=contour.col,aes(x,y,z=z),size=1.4) +
      theme(legend.title=element_blank(),legend.key.width = unit(0.35,"cm"),
            legend.text = element_text(size = 11), legend.spacing.x = unit(0.2, 'cm'),
            legend.key.height = unit(1,"cm"), axis.text=element_text(size=11),
            axis.title=element_text(size=12,face="bold") ) +
      scale_fill_gradientn(colours = tim.colors(64),na.value="white") 
  }
  #  quilt.plot(gr$x, gr$y, pmin(ub,obj),xlim=range(x),ylim=range(y),
  #             xlab=expression(X[t]^1), ylab=expression(X[t]^2),cex.lab=1.2, cex.axis=1.1,...)
  #imx <- as.image(x=cbind(gr$x,gr$y),obj,nr=100,nc=100)
  #contour(imx$x,imx$y,imx$z,levels=0,add=T,drawlab=F,lwd=2,col=contour.col)
  #if (class(fit)=="dynaTree") { #is.null(fit@X) == F & is(fit,"km")==F) {
  #  kk <- kde2d(fit$X[,1],fit$X[,2],lims=c(range(x),range(y)))
  #  contour(kk$x, kk$y, kk$z,nlev=12,lwd=0.5,lty=2,drawlab=F,add=T,col=contour.col)
  #}
  if (is(fit,"km") & (only.contour== FALSE))
    g1 <- g1 + geom_point(data=data.frame(xx=fit@X[,1],yy=fit@X[,2]),color="orange", aes(xx,yy), size=2)
  if( (class(fit)=="hetGP" | class(fit)=="homGP" | class(fit) == "homTP") & (only.contour== FALSE))
    g1 <- g1 + geom_point(data=data.frame(xx=fit$X0[,1],yy=fit$X0[,2]),color="orange", aes(xx,yy), size=2)
    
  return(g1)
 
}



######################################
#' @title Create a Bouchard-Warin equal prob grid
#'
#' @details Recursively sort along each of the d-coordinates
#' At the end do local linear regression at each leaf
#' This is a recursive algorithm!
#' first column is reserved for the y-coordinate (timingValue)
#' It's safest if nrows(grid) is divisible by \code{model$nChildren^model$dim}
#' @param curDim dimension of the grid
#' @param model a list containing all model parameters. In particular must have \code{model$nChildren}
#' @param grid dataset of x-values
#' @param test testPaths to predict along as well
#' @export
treeDivide.BW <- function(grid,curDim, model,test)
{
  lenChild <- round(dim(grid)[1]/model$nChildren)
  fit.in <- array(0, dim(grid)[1])
  fit.test <- array(0, dim(test)[1])
  
  
  if (curDim <= (model$dim+1))
  {
    sorted <- sort(grid[,curDim],index.ret=T)$ix
    for (j in 1:model$nChildren) {
      ndxChild <- sorted[ (1+(j-1)*lenChild):min((j*lenChild),length(sorted))]
      # on the test data the index is shifted back by 1 since there is no 'y' column
      test.ndx <- which( test[,curDim-1] > grid[sorted[ (1+(j-1)*lenChild)],curDim] & test[,curDim-1] <= grid[sorted[j*lenChild],curDim])
      fitc <- treeDivide.BW(grid[ndxChild,], curDim+1, model, test[test.ndx,,drop=F])
      # save the in-sample and the out-of-sample results
      fit.in[ndxChild] <- fitc$in.sample
      fit.test[test.ndx] <- fitc$out.sample
    }
  }
  else  {
    #fit <- predict(earth(grid,degree=2,nk=100,thresh=1e-8,nprune=10))
    fit <- lm(grid)
    fit.in <- predict(fit)
    fit.test <- predict(fit, data.frame(grid=test))
  }
  return (list(in.sample=fit.in,out.sample=fit.test))
}


############################################
#' 1-d version of \code{\link{treeDivide.BW}} that stores all the local fits
#' @inheritParams treeDivide.BW
#' @export
treeDivide.BW.1d <- function(grid,curDim, model, test)
{
  lenChild <- floor(dim(grid)[1]/model$nChildren)
  fitchild <- array(0, dim(grid)[1])
  fittest <- array(0, length(test))
  lm.model <- list()
  bounds <- array(0, model$nChildren)
  
  sorted <- sort(grid[,curDim],index.ret=T)$ix
  for (j in 1:model$nChildren) {
    ndxChild <- sorted[ (1+(j-1)*lenChild):(j*lenChild)]
    test.ndx <- which( test > grid[sorted[ (1+(j-1)*lenChild)],curDim] & test <= grid[sorted[j*lenChild],curDim])
    
    bounds[j] <- max( grid[ndxChild,curDim])
    mod <- lm(grid[ndxChild,])
    lm.model[[j]] <- as.numeric(mod$coef)
    fitchild[ndxChild] <- predict(mod)
    fittest[test.ndx] <- predict(mod, data.frame(grid=test[test.ndx]))
  }
  return (list(in.sample=fitchild,out.sample=fittest,lm.model=lm.model,boundaries=bounds))
}
######################################################################################

###############
#' Expected Loss for Contour Finding
#'
#' @title Compute expected loss using the optimal stopping loss function.
#' @param objMean predicted mean response
#' @param objSd posterior standard deviation of the response
#' @references 
#' Mike Ludkovski, Kriging Metamodels and Experimental Design for Bermudan Option Pricing
#'  Journal of Computational Finance, 22(1), 37-77, 2018
#' @author Mike Ludkovski
#' @export
cf.el <- function(objMean,objSd)
{
  el <- pmax(0, objSd*dnorm( -abs(objMean)/objSd ) - abs(objMean)*pnorm( -abs(objMean)/objSd) )
  return (el)
}

#####################
#' SUR (Stepwise Uncertainty Reduction) acquisition function for Contour Finding from Ludkovski (2018)
#'
#' @title Compute EI for Contour Finding using the ZC-SUR formula
#' @param objMean predicted mean response
#' @param objSd posterior standard deviation of the response
#' @param nugget GP nugget parameter, see under details
#' @details   compute the change in ZC = sd*(1-sqrt{(nugget)})/sqrt{(nugget + sd^2)}
#' @references 
#' Mike Ludkovski, Kriging Metamodels and Experimental Design for Bermudan Option Pricing
#'  Journal of Computational Finance, 22(1), 37-77, 2018
#' @author Mike Ludkovski
#' @export
cf.sur <- function(objMean, objSd, nugget)
{
  a = abs(objMean)/objSd  # normalized distance to zero contour
  M = objMean*pnorm(-a)-objSd*dnorm(-a)
  var_reduce <- objSd*(1-sqrt(nugget)/sqrt(nugget+objSd^2))  # reduction in posterior stdev from a new observation
  newSigma <- (objSd - var_reduce)   # look-ahead variance
  a_new <- abs(objMean)/newSigma   # new distance to zero-contour
  M_new = objMean*pnorm(-a_new)-newSigma*dnorm(-a_new)  # new ZC measure
  # difference between next-step ZC and current ZC
  return( M_new-M)
}

#####################
#' cSUR for Contour Finding based on Ludkovski (2018)
#'
#' @title Compute reduction in contour-distance
#' @param objMean predicted mean response
#' @param objSd posterior standard deviation of the response
#' @param nugget the noise variance to compute the ALC factor
#' @details   compute the change in ZC = sd*(1-sqrt{(nugget)})/sqrt{(nugget + sd^2)}
#' @references 
#' Mike Ludkovski, Kriging Metamodels and Experimental Design for Bermudan Option Pricing
#'  Journal of Computational Finance, 22(1), 37-77, 2018
#' @author Mike Ludkovski
#' @seealso [osp.seq.design]
#' @export
cf.csur <- function(objMean, objSd, nugget)
{
  a = pnorm(-abs(objMean)/objSd)  # normalized distance to zero contour
  var_reduce <- objSd*(1-sqrt(nugget)/sqrt(nugget+objSd^2))  # reduction in posterior stdev from a new observation
  newSigma <- (objSd - var_reduce)   # look-ahead variance
  a_new <- pnorm(-abs(objMean)/newSigma)   # new distance to zero-contour
  # difference between next-step ZC and current ZC
  return( a-a_new)
}


#####################
#' MCU for Contour Finding. DEPRECATED.
#'
#' @title Maximum Contour Uncertainty criterion
#' @param objMean predicted mean response
#' @param objSd posterior standard deviation of the response
#' @details   compute normalized distance to zero-contour |mu|/sd
#' @references 
#' Mike Ludkovski, Kriging Metamodels and Experimental Design for Bermudan Option Pricing
#'  Journal of Computational Finance, 22(1), 37-77, 2018
#' @author Mike Ludkovski
#' @seealso [osp.seq.design]
#' @export
cf.mcu <- function(objMean, objSd)
{
  return(pnorm(-abs(objMean)/objSd))
}

#####################
#' straddle MCU with a specified variance weight
#'
#' @title Straddle Maximum Contour Uncertainty criterion
#' @param objMean predicted mean response
#' @param objSd posterior standard deviation of the response
#' @param gamma weight on the variance
#' @details   compute the UCB criterion with constant weight: gamma*s(x) - |f(x)|
#' @references 
#' Mike Ludkovski, Kriging Metamodels and Experimental Design for Bermudan Option Pricing
#'  Journal of Computational Finance, 22(1), 37-77, 2018
#'  
#'  X.Lyu, M Binois, M. Ludkovski (2020+) Evaluating Gaussian Process Metamodels and Sequential Designs for 
#'  Noisy Level Set Estimation <https://arxiv.org/abs/1807.06712>
#' @author Mike Ludkovski
#' @seealso [osp.seq.design]
#' @export
cf.smcu <- function(objMean, objSd, gamma=1.96)
{
  return( gamma*objSd - abs(objMean) )
}

#####################
#' tMSE for Contour Finding
#'
#' @title targeted Mean Squared Error criterion
#' @param objMean predicted mean response
#' @param objSd posterior standard deviation of the response
#' @param seps epsilon in the tMSE formula. By default taken to be zero.
#' @details   compute predictive density at the contour, smoothed by seps
#' 
#' @references 
#' Mike Ludkovski, Kriging Metamodels and Experimental Design for Bermudan Option Pricing
#'  Journal of Computational Finance, 22(1), 37-77, 2018\cr
#'  
#'  X.Lyu, M Binois, M. Ludkovski (2020+) Evaluating Gaussian Process Metamodels and Sequential Designs for 
#'  Noisy Level Set Estimation <https://arxiv.org/abs/1807.06712>
#'  
#' @author Mike Ludkovski
#' @seealso [osp.seq.design]
#' @export
cf.tMSE <- function(objMean, objSd, seps = 0)
{
  
  w <- 1/sqrt(2 * pi * (objSd^2 + seps)) * exp(-0.5 * (objMean / sqrt(objSd^2 + seps))^2)
  
  return(w * objSd^2)
}

#################################
# more flexible forward.sim
policy.payoff <- function( x,M,fit,model,offset=1,path.dt=model$dt,use.qv=FALSE,interp=0)
{
  nsim <- 0
  if (is(x, "matrix") | is(x, "numeric") ) {
    curX <- model$sim.func( x,model,path.dt)
    nsim <- nrow(x)
  }
  if (is(x, "list" ) )
    curX <- x[[1]]
  payoff <- rep(0, nrow(curX))
  tau <- rep(0, nrow(curX))
  fvalue <- vector(mode="list",len=model$look.ahead+1)
  
  contNdx <- 1:nrow(curX)
  i <- 1
  payoff[contNdx]  <- exp(-(i)*path.dt*model$r)*model$payoff.func( curX[contNdx,,drop=F], model)
  
  # main loop forward
  while (i < (M+(use.qv==TRUE)) & length(contNdx) > 0) {
    
    
    in.the.money <-  which( payoff[contNdx] > 0) # 1:length(contNdx) #
    stop.possible <- (interp > 0) + (i %% (model$dt/path.dt) == 0)
    #if (interp > 0 & ceiling( (i+1-offset)/(model$dt/path.dt) ) == model$T/model$dt)
    #  stop.possible <- 0
    if ( length(fit) <= (2-offset))  # no actual rules setup yet
       stop.possible <- 0
    
    
    if (length(in.the.money) >0 & stop.possible > 0) {
      if (interp == 0)
         fit.ndx <- i/(model$dt/path.dt) + 1-offset
      if (interp == 1 | interp == 2)
         fit.ndx <- min(ceiling( (i)/(model$dt/path.dt) )+1-offset, length(fit)-1)
      
      if (is(fit[[fit.ndx]],"earth") )
        rule <- predict(fit[[fit.ndx]],curX[contNdx[in.the.money],,drop=F]) # for use with  MARS
      if (is(fit[[fit.ndx]],"deepnet") )
        rule <- nn.predict(fit[[fit.ndx]],curX[contNdx[in.the.money],,drop=F]) # for use with  deepnet
      if (is(fit[[fit.ndx]],"smooth.spline") )
        rule <- predict(fit[[fit.ndx]],curX[contNdx[in.the.money],,drop=F])$y # for use with  splines
      if (is(fit[[fit.ndx]],"randomForest") ) {
        obj <- predict(fit[[fit.ndx]],curX[contNdx[in.the.money],,drop=F],predict.all=T)$individual
        rule <- apply(obj,1,median)
      }
      if (is(fit[[fit.ndx]],"dynaTree")  ){
        obj <- predict(fit[[fit.ndx]], curX[contNdx[in.the.money],,drop=F],quants=F)
        rule <- obj$mean
      }
      if (is(fit[[fit.ndx]],"km") ) # DiceKriging
        rule <- predict(fit[[fit.ndx]],data.frame(x=curX[contNdx[in.the.money],,drop=F]),type="UK")$mean
      if (is(fit[[fit.ndx]],"gpi") )  # laGP
        rule <- predGP(fit[[fit.ndx]],XX=curX[contNdx[in.the.money],,drop=F], lite=TRUE)$mean
      if (is(fit[[fit.ndx]],"lm") ) {
        lenn <- length(fit[[fit.ndx]]$coefficients)
        rule <-  fit[[fit.ndx]]$coefficients[1] + model$bases(curX[contNdx[in.the.money],,drop=F]) %*%
          fit[[fit.ndx]]$coefficients[2:lenn]
      }
      if( class(fit[[fit.ndx]])=="homGP" | class(fit[[fit.ndx]]) == "hetGP" | class(fit[[fit.ndx]]) == "homTP")
        rule <- predict(x=curX[contNdx[in.the.money],,drop=F], object=fit[[fit.ndx]])$mean
      if( class(fit[[fit.ndx]])=="rvm")
        rule <-  predict(fit[[fit.ndx]], new=curX[contNdx[in.the.money],,drop=F])
      if (class(fit[[fit.ndx]]) == "npregression")
        rule <- predict(fit[[fit.ndx]], new=curX[contNdx[in.the.money],,drop=F])
      
      if (use.qv == TRUE & i== M) {
        payoff[contNdx] = payoffCont[contNdx] + rule  # continuation value of paths that didn't stop yet)
        break
      }
      
      if (interp == 2 & fit.ndx <= (length(fit) - 1) & fit.ndx > (1+1-offset) ) {
        if (i > (length(fit)-(1-offset)-1)*model$dt/path.dt ) # (model$T-model$dt)/path.dt) 
        {
          #rule <- predict(x=curX[contNdx[in.the.money],,drop=F], object=fit[[fit.ndx]])$mean
          rule3 <- exp(-(i)*path.dt*model$r)*model$payoff.func( curX[contNdx[in.the.money],,drop=F])
          weight <- (i*path.dt - model$T+model$dt)/model$dt
          rule <- (rule3+ rule)*(1-weight) - rule3
          
        }
        else if (i %% (model$dt/path.dt) > 0) {
          rule2 <- predict(x=curX[contNdx[in.the.money],,drop=F], object=fit[[fit.ndx-1]])$mean
          weight <- (i %% (model$dt/path.dt))/(model$dt/path.dt)
          rule3 <- rule
          rule <- rule3*weight + rule2*(1-weight)
        }
      } 
      
      # stop if the expected gain is negative
      #contNdx <- contNdx[ which (rule > 0 | payoff[contNdx] == 0)]
      if (length(which(rule <0)) > 0)
        contNdx <- contNdx[-in.the.money[which(rule < 0)]  ]
    }
    tau[contNdx] <- (i)*path.dt
    
   
    # update the x values by taking a step of length dt
    i <- i+1
    
    if (is(x, "matrix") | is(x, "numeric") ) {
      curX[contNdx,] <- model$sim.func( curX[contNdx,,drop=F],model,path.dt)
      nsim <- nsim + length(contNdx)
    }
    
    if (is(x, "list") )  # stored list of paths
      curX[contNdx,] <- x[[i]][contNdx,,drop=F]
    # payoff for next timestep -- used for terminal payoff at i=M
    payoff[contNdx]  <- exp(-(i)*path.dt*model$r)*model$payoff.func( curX[contNdx,,drop=F], model)
  }
  if (model$look.ahead > 1)
    for (i in 2:(model$look.ahead))   # payoff for a trajectory starting at x^n_{t+i} which was still alive then
      fvalue[[i]] <- payoff[ save.ndx[[i]] ]*exp((i-1)*path.dt*model$r)
  
  
  return( list(payoff=payoff,tau=tau, nsims=nsim))
  # payoff is the resulting payoff NPV from t=0
  # fvalue[i] is a list of resulting payoffs (on paths still not stopped) NPV from t=i
  # tau are the times when stopped
  # sims is a list; sims[[i]] are the forward x-values of paths at t=i (those not stopped yet)
}

############################################
#' Simulate \eqn{\sum_k h(X_{tau_k})} using \code{fit} emulators
#' @title Forward simulation of a swing payoff based on a sequence of emulators
#' @param x   a matrix of starting values (N x \code{model$dim}).
#' If input \code{x} is a list, then use the grids specified by x
#' @param M     number of time steps to forward simulate
#' @param fit   a list of fitted emulators that determine the stopping classifiers to be used
#' @param offset deprecated
#' @param model List containing all model parameters. In particular uses \code{model$dt,model$r} 
#' for discounting and \code{model$swing.payoff} to compute payoffs
#' @param use.qv experimental, do not use
#' @param verbose for debugging purposes
#' @param n.swing number of swing rights (integer, at least 1)
#' @export
#' @return a list containing:
#' \itemize{
#'  \item \code{payoff}: a vector of length `nrow(x)` containing the resulting payoffs NPV from $t=0$
#'  \item \code{tau} matrix of the times when stopped. Columns represent the rights exercised
#'  \item  \code{nsims} number of total 1-step simulations performed
#' }
#' @details Should be used in conjuction with the \code{\link[mlOSP]{swing.fixed.design}} function that builds the emulators. 
swing.policy <- function( x,M,fit,model,offset=1,use.qv=FALSE,n.swing=1,verbose=FALSE)
{
  nsim <- 0
  if (is(x, "matrix") | is(x, "numeric") ) {
    curX <- model$sim.func( x,model,model$dt)
    nsim <- nrow(x)
  }
  if (is(x, "list" ) )
    curX <- x[[1]]
  payoff <- array(0, dim=c(nrow(curX), n.swing) )
  tau <- array(0, dim=c(nrow(curX), n.swing) )
  refract <- rep(0, nrow(curX))
  refractN <- model$refract/model$dt

  i <- 1
  if (n.swing <= 0)
    return( list(totPayoff=rep(0, nrow(curX))) )
  
  rightsLeft <- rep(n.swing, nrow(curX))  # remaining rights for each path
  
  # main loop forward
  while (i < (M+(use.qv==TRUE))) {
    for (kk in 1:n.swing) {  # loop over number of remaining paths
      curNdx <- which( rightsLeft == kk)
      if (length(curNdx) == 0)
        next
      # immediate payoff
      myFit <- fit[[i+1-offset,kk]]; 
      myx <- curX[curNdx,,drop=F]
      
      if (is(myFit,"earth") ) {
        rule <- predict(myFit,myx) # for use with  MARS
        
      }
      if (is(myFit,"deepnet") ) {
        rule <- nn.predict(myFit,myx) # for use with deepnet
        
      }
      if (is(myFit,"nnet") ) {
        rule <- predict(myFit,myx,type="raw") # for use with nnet
        
      }
      if (is(myFit,"smooth.spline") ) {
        rule <- predict(myFit,myx)$y # for use with  splines
        
      }
      if (is(myFit,"randomForest") ) {
        obj <- predict(myFit,myx,predict.all=T)$individual
        rule <- apply(obj,1,median)  ############### randomForest uses median
        
      }
      if (is(myFit,"dynaTree")  ){
        obj <- predict(myFit, myx,quants=F)
        rule <- obj$mean
        
        
      }
      if (is(myFit,"km") ) { # DiceKriging
        rule <- predict(myFit,data.frame(x=myx),type="UK")$mean
        
      }
      if (is(myFit,"gpi") ) { # laGP
        rule <- predGP(myFit,XX=myx, lite=TRUE)$mean
        
      }
      if (is(myFit,"lm") ) {
        lenn <- length(myFit$coefficients)
        rule <-  myFit$coefficients[1] + model$bases(myx) %*% myFit$coefficients[2:lenn]
        
      }
      if( class(myFit)=="homGP" | class(myFit) == "hetGP" | class(myFit) == "homTP") {
        rule <- predict(x=myx, object=myFit)$mean
        
      }
      if( class(myFit)=="rvm") {
        rule <-  predict(myFit, new=myx)
        
      }
      if (class(myFit) == "npregression") {
        rule <- predict(myFit, new=myx)
        
      }
      
      #if (i < M-refractPeriod) {
      #  delayedFit <- fitVplusDelta[[min(M, i+1+refractN-offset),kk]]
      #  refractPayoff <- ospPredict(delayedFit, myx, model)
      #}
      #else 
      #  refractPayoff <- 0
      
      #imm  <- exp(-(i)*model$dt*model$r)*model$swing.payoff(myx , model) + refractPayoff
      imm <- exp(-(i)*model$dt*model$r)*model$swing.payoff(myx , model) 
      
      # exercise if rule < 0 and no refraction left
      if (length(which(rule <0 & refract[curNdx] == 0)) > 0) {
        stopKK <- which(rule <0 & refract[curNdx] == 0)
        payoff[ curNdx[stopKK], n.swing-kk+1] <- imm[stopKK]
        tau[ curNdx[stopKK], n.swing-kk+1] <- i*model$dt
        rightsLeft[ curNdx [stopKK]] <- kk - 1
        refract[ curNdx[ stopKK]] <- model$refract
      }
    }  #loop over kk
  
    # update the x values by taking a step of length dt
    i <- i+1
    if (verbose==TRUE)
       browser()
    
    if (is(x, "matrix") | is(x, "numeric") ) {
      curX <- model$sim.func( curX,model,model$dt)
      nsim <- nsim + nrow(curX)
    }
    refract <- pmax(0, refract - model$dt)
    
    if (is(x, "list") )  # stored list of paths
      curX <- x[[i]]
    
  }
  curNdx <- which( refract == 0 & rightsLeft > 0)
  totPayoff = apply(payoff,1,sum)
  totPayoff[curNdx] <- totPayoff[curNdx] + exp(-(M)*model$dt*model$r)*model$swing.payoff( curX[curNdx,,drop=F], model)

  return( list(payoff=payoff,sims=nsim,tau=tau,totPayoff = totPayoff))
  # payoff is the resulting payoff NPV from t=0
  # fvalue[i] is a list of resulting payoffs (on paths still not stopped) NPV from t=i
  # tau are the times when stopped
  # sims is a list; sims[[i]] are the forward x-values of paths at t=i (those not stopped yet)
}


###########
ospPredict <- function(myFit,myx,model)
{
  if (is(myFit,"earth") ) {
    prediction <- predict(fit,myx) # for use with  MARS
  }
  if (is(myFit,"deepnet") ) {
    prediction <- nn.predict(fit,myx) # for use with deepnet
  }
  if (is(myFit,"nnet") ) {
    prediction <- predict(fit,myx,type="raw") # for use with nnet
  }
  if (is(myFit,"smooth.spline") ) {
    prediction <- predict(myFit,myx)$y # for use with  splines
  }
  if (is(myFit,"randomForest") ) {
    prediction <- predict( myFit, myx)
  }
  if (is(myFit,"dynaTree")  ){
    prediction <- predict( myFit, myx, quants=F)$mean
    
  }
  if (is(myFit,"km") ) { # DiceKriging
    prediction <- predict(myFit,data.frame(x=myx),type="UK")$mean
  }
  if (is(myFit,"gpi") ) { # laGP
    prediction <- predGP(myFit,XX=myx, lite=TRUE)$mean
  }
  if (is(myFit,"lm") ) {
    lenn <- length(myFit$coefficients)
    prediction <-  myFit$coefficients[1] + model$bases(myx) %*% myFit$coefficients[2:lenn]
  }
  if( class(myFit)=="homGP" | class(myFit) == "hetGP" | class(myFit) == "homTP") {
    prediction <- predict(x=myx, object=myFit)$mean
  }
  if( class(myFit)=="rvm") {
    prediction <-  predict(myFit, new=myx)
  }
  if (class(myFit) == "npregression") {
    prediction <- predict(myFit, new=myx)
  }
  
  return( prediction)
  
  
}  

## Helper for liGP
shift_scale <- function(shift_scale_list = NULL, scale_list = NULL, shift = NULL, scale = NULL){
  
  if (!is.null(shift) & !is.null(scale)){
    if (!is.null(shift_scale_list)){
      dim <- ncol(shift_scale_list[[1]])
    } else {
      dim <- ncol(scale_list[[1]])
    }
    
    if (length(shift) == 1) shift <- rep(shift, dim)
    if (length(scale) == 1) scale <- rep(scale, dim)
    
    for(k in 1:dim) {
      if (!is.null(shift_scale_list))
        for(l in 1:length(shift_scale_list))
          shift_scale_list[[l]][,k] <- (shift_scale_list[[l]][,k]-shift[k])/scale[k]
        if (!is.null(scale_list))
          for(l in 1:length(scale_list))
            scale_list[[l]][,k] <- scale_list[[l]][,k]/scale[k]
    }
  }
  
  return(c(shift_scale_list, scale_list))
}


## find_shift_scale_stats:
##
## Normalizes the input dimensions, then fits a global GP
## to a subset of the data (up to gp_size points). The 
## square root of the lengthscales are used to scale the
## normalized inputs. Returns the shift and scale parameters 
## for both the inputs (X) and response (Y).
find_shift_scale_stats <- function(X, Y, gp_size = 1000){
  Xsc <- X
  Ysc <- (Y - mean(Y))/sd(Y)
  
  xmins <- apply(X, 2, min)
  xmaxs <- apply(X, 2, max)
  
  for(k in 1:ncol(X)) 
    Xsc[,k] <- (X[,k]-xmins[k])/(xmaxs[k]-xmins[k])
  
  ## Fit global model to subset of data to get lengthscales
  nsub <- min(nrow(X), gp_size)
  d2 <- darg(NULL, Xsc)
  subs <- sample(1:nrow(X), nsub, replace=FALSE)
  
  gpsi <- mleHomGP(Xsc[subs,], Ysc[subs], lower = rep(5e-3, ncol(X)), 
                   upper=rep(2, ncol(X))) #, known=list(g=1e-4))
  init.th <- gpsi$theta
  
  ## Set shift and scale parameters
  sd <- sqrt(init.th)*(xmaxs-xmins)
  #browser()
  mu <- xmins
  scale_list <- list(xshift = mu, xscale = sd, yshift = mean(Y) + sd(Y)*gpsi$beta0, yscale = sd(Y))
  
  return(scale_list)
}


############### Impulse for Forest rotation
# Impulse function for forest
forest.impulse <- function(cur_x, model, fit)
{
  payoff <- ( cur_x - model$impulse.fixed.cost - model$impulse.target)
  if (is.null(fit))
     next_step_value <- 0
  else
    next_step_value <- predict(fit, model$impulse.target)$y
  
  return( payoff + next_step_value*exp(-model$dt*model$r))
  
}

forward.impulse.policy <- function( x,M,fit,model)
{
  curX <- model$sim.func( x,model,model$dt)
  payoff <- rep(0, nrow(curX))
  tau <- rep(0, nrow(curX))
  paths <- array(0, dim=c(nrow(curX),ncol(curX),M))
  paths[,,1] <- curX
  imps <- array(NaN, dim=c(nrow(curX),M))
  bnd <- rep(0,M)

  i <- 1
  # main loop forward
  while (i < M) {
    
    if (is(fit[[i]],"smooth.spline") )
        cont <- exp(-model$dt*model$r)*predict(fit[[i]],curX)$y # for use with  splines
   
      if (is(fit[[i]],"lm") ) {
        lenn <- length(fit[[i]]$coefficients)
        cont <-  fit[[i]]$coefficients[1] + model$bases(curX) %*%
          fit[[i]]$coefficients[2:lenn]
      }
      impulse <- forest.impulse(curX,model,fit[[i]])
      impNdx <- which(cont < impulse & curX > model$impulse.fixed.cost)
      if (length(impNdx) > 0) {
         imps[impNdx,i] <- curX[impNdx]
         bnd[i] <- min(curX[impNdx])
         payoff[impNdx] <- payoff[impNdx] + exp(-model$r*model$dt*i)*pmax(curX[impNdx] - model$impulse.fixed.cost- model$impulse.target, 0)
         curX[impNdx] <- model$impulse.target  # reset impulsed x
         tau[impNdx] <- tau[impNdx] + 1
         
      }
    
    i <- i+1
    
    curX <- model$sim.func( curX,model,model$dt )
    paths[,,i] <- curX

  }
  # terminal time
  payoff <- payoff + exp(-model$r*model$dt*M)*pmax(curX - model$impulse.fixed.cost- model$impulse.target, 0)
  
  return( list(payoff=payoff,tau=tau,paths=paths,impulses=imps,bnd=bnd))
  # payoff is the resulting payoff NPV from t=0
}

#' Initial design for the 2D Bermudan Put example
#'
#' A dataset containing 20 initial locations to be used with \code{osp.seq.design}. 
#' The Put has strike K=40 and i.i.d assets, so the experimental design is roughly symmetric.
#' The design is space-filling a triangular in-the-money region
#' 
#'
#' @format An array with 20 rows and 2 columns:
#' \describe{
#'   \item{row}{unique designs}
#'   \item{columns}{values of corresponding x_1 and x_2 coordinates}
#' }
"int_2d"

#' Initial design for the 3D Bermudan max-Call example
#'
#' A dataset containing 300 initial locations to be used with \code{osp.seq.design}. 
#' The max Call has strike K=100 and i.i.d assets, so the experimental design is roughly symmetric.
#' The design is space-filling a hyper-rectangular in-the-money region
#' 
#'
#' @format An array with 300 rows and 3 columns:
#' \describe{
#'   \item{row}{unique designs}
#'   \item{columns}{values of corresponding x_1, x_2 and x_3 coordinates}
#' }
"int300_3d"