############################################
#' Simulate h(X_tau) using FIT (can be a dynaTree, smooth.spline, MARS or RandomForest or RVM or hetGP)
#' @title Forward simulation based on a sequence of emulators
#' @param x     is a matrix of starting values
#'
#' if input \code{x} is a list, then use the grids specified by x
#' @param M     is number of time steps to forward simulate
#' @param fit   is a list of fitted emulators that determine the stopping classifiers to be used
#' @param use.qv boolean to indicate whether to plug-in continuation value for all remaining paths at the last time-step
#' default is set to \code{FALSE}
#' @export
#' @return a list containing:
#' \itemize{
#'  \item \code{payoff} is the resulting payoff NPV from t=0
#'  \item \code{fvalue[i]} is a list of resulting payoffs (on paths still not stopped) NPV from t=i
#'  \item \code{tau} are the times when stopped
#'  \item \code{sims} is a list; \code{sims[[i]]} are the forward x-values of paths at t=i (those not stopped yet)
#' \code{nsims} number of total 1-step simulations performed
#' }
#' @details Should be used in conjuction with the \code{osp.} functions that build the emulators. Also called
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
    if (length(in.the.money) >0 & is.null( fit[[[i+1-offset]) == FALSE) {
      
      if (is(fit[[i+1-offset]],"earth") )
        rule <- predict(fit[[i+1-offset]],curX[contNdx[in.the.money],,drop=F]) # for use with  MARS
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
      if (is(fit[[i+1-offset]],"gpi") )  # laGP
        rule <- predGP(fit[[i+1-offset]],XX=curX[contNdx[in.the.money],,drop=F], lite=TRUE)$mean
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
      
      if (use.qv == TRUE & i== M) {
        payoff[contNdx] = payoff[contNdx] + rule  # continuation value of paths that didn't stop yet
        #browser()
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
    # update the x values by taking a step of length dt
    i <- i+1
    
    if (is(x, "matrix") | is(x, "numeric") ) {
      curX[contNdx,] <- model$sim.func( curX[contNdx,,drop=F],model,model$dt)
      nsim <- nsim + length(contNdx)
    }
    
    if (is(x, "list") )  # stored list of paths
      curX[contNdx,] <- x[[i]][contNdx,,drop=F]
    # payoff for next timestep -- used for terminal payoff at i=M
    payoff[contNdx]  <- exp(-(i)*model$dt*model$r)*model$payoff.func( curX[contNdx,,drop=F], model)
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
#' Two-dimensional raster+contour+point plot
#' @title Visualize 2D emulator + stopping region
#' @details Uses the raster plot from \pkg{ggplot2}
#'
#' @param x,y locations to use for the \code{predict()} functions. Default is a 200x200 fine grid.
#' Passed to \code{expand.grid}
#' @param fit can be any of the types supported by \code{\link{forward.sim.policy}}
#' @param show.var -- if \code{TRUE} then plot posterior surrogate variance instead of surrogate mean [default = FALSE]
#' This only works for \code{km} and \code{het/homGP/homTP} objects
#' @param only.contour -- just add the zero-contour, no raster plot (uses \code{contour}, no ggplot)
#' @param ub clip the surface with an upper bound  to see the zero-contour better
#' @param contour.col (default is "red") -- color of the zero contour
#' @export
plt.2d.surf <- function( fit, x=seq(31,43,len=201),y = seq(31,43,len=201),ub=1e6,
                         show.var=FALSE, only.contour=FALSE, contour.col="red")
{
  gr <- expand.grid(x=x,y=y)
  
  if (is(fit,"randomForest") ) {
    obj <- predict(fit,cbind(gr$x,gr$y),predict.all=T)$individual
    obj <- apply(obj,1,median)
  }
  if (is(fit,"dynaTree")  )
    obj <- predict(fit,cbind(gr$x,gr$y),quants=F)$mean
  if (is(fit,"lm") ) {
    obj <-  fit$coefficients[1] + model$bases(cbind(gr$x,gr$y)) %*%
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
  
  if (only.contour==TRUE) {
    imx <- as.image(x=cbind(gr$x,gr$y),obj,nr=100,nc=100)
    contour(imx$x,imx$y,imx$z,levels=0,add=T,drawlab=F,lwd=2,col=contour.col)
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
      scale_fill_gradientn(colours = tim.colors(64)) 
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

################################
#' Make a 3-panel plot of the 1-d fit from mc.put.adaptree
#' @param show.ci can be 'none', 'km', 'dt-vmean', 'dt-df'
#' @param new.fig : whether to create a fresh panel
#' @param heuristic: calls the corresponding qEI function
#' @export
plt.1d.fit <- function(obj,x.gr = seq(28,40,len=500),new.fig=F,show.ci='none',heuristic='zc',al.weights=NULL,obj2=NULL,...)
{
  al.sort <- sort(x.gr, index.ret=T)
  sort.ndx <- al.sort$ix
  al.sort <- al.sort$x
  if (is(obj, "dynaTree") ) {
    show.fit <- predict(obj, al.sort,quants=F)
    show.fit <- alc(show.fit,al.sort)
  }
  if (is(obj, "km"))
    show.fit <- predict(obj,data.frame(x=al.sort), type="UK")
  
  
  
  if (heuristic == 'sgn')
    show.weights <- qEI.sgn(show.fit)
  if (heuristic == 'rbm')
    show.weights <- qEI.rbm(show.fit,alph=1.5)
  if (heuristic == 'zc')
    show.weights <- qEI.zc(show.fit)
  if (heuristic == 'squared')
    show.weights <- qEI.sq(show.fit)
  if (heuristic == 'given')
    show.weights <- al.weights[sort.ndx]
  
  
  if (new.fig == T)
    dev.new()
  
  if (is.null(show.fit$X) == F & show.ci == 'none')
    par(mfrow=c(3,1),mar=c(3.8,4.7,2.1,0.5), oma=c(1.1,1.1,0.5,1), bty="n")
  else
    par(mfrow=c(2,1),mar=c(3.8,4.7,2.1,0.5), oma=c(1.1,1.1,0.5,1), bty="n")
  plot(al.sort,show.fit$mean,type="l",lwd=2,ylab="Timing Value T(t,x)",xlab=expression(S[t]),xlim=range(x.gr),ylim=c(-0.4,1),
       cex.axis=1.8, cex.lab=1.8,...)
  if (show.ci == 'dt-vmean')
    sd <- pmax(sqrt(show.fit$vmean),sqrt(show.fit$var/show.fit$df))
  if (show.ci == 'dt-df')
    sd <- sqrt(show.fit$var/show.fit$df)
  if (show.ci == 'km')
    sd <- show.fit$sd
  
  if (show.ci != 'none') {
    polygon( c(al.sort, rev(al.sort)), c(show.fit$mean-1.96*sd, rev(show.fit$mean+1.96*sd)),
             col=rgb(t(c((col2rgb('hotpink3')/256+3)/4)),alph=0.7), border=NA)
    
    if (is(obj2, "km")) {
      show.bench <- predict(obj2,data.frame(x=al.sort), type="UK")
      lines( al.sort, show.bench$mean,lty=1, lwd=1.5)
    }
    lines(al.sort,show.fit$mean, col="hotpink3",lwd=3)
    
  }
  abline(h=0,lty=2)
  #sgn.weights <- qEI.sgn(show.fit,anneal=1)
  #lines(al.sort,sgn.weights,col="blue",lty=2)
  if (is.null(show.fit$X) == F & show.ci == 'none') {
    plot(al.sort,show.weights, xlim=range(x.gr),type="l",ylab=expression(w[n](t,x)),xlab="", cex.axis=1.8,cex.lab=1.8)
    hist(show.fit$X,40,freq=F, yaxt="n", xlim=range(x.gr),xlab=expression(X[t]),main="",cex.axis=1.8,cex.lab=1.8)
  }
  else
    plot(al.sort,show.weights, xlim=range(x.gr),type="l",ylab=expression(L[n](x)),xlab="", ylim=c(0,0.02),cex.axis=1.8,cex.lab=1.8)
  
}


######################################
#' Create a Bouchard-Warin equal prob grid
#'
#' Recursively sort along each of the d-coordinates
#' At the end do local linear regression at each leaf
#' This is a recursive algorithm!
#' first column is reserved for the y-coordinate (timingValue)
#' It's safest if nrows(grid) is divisible by model$nChildren^dim
#' @param curDim dimension of the grid
#' @param grid -- dataset of x-values
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
######################################################################################

###############
#' Expected Loss for Contour Finding
#'
#' @title Compute expected loss using the optimal stopping loss function.
#' @param objMean: predicted mean response
#' @param objSd: posterior standard deviation of the response
#' @export
cf.el <- function(objMean,objSd)
{
  el <- pmax(0, objSd*dnorm( -abs(objMean)/objSd ) - abs(objMean)*pnorm( -abs(objMean)/objSd) )
  return (el)
}

#####################
#' SUR for Contour Finding
#'
#' @title Compute EI for Contour Finding using the ZC-SUR formula
#' @param objMean: predicted mean response
#' @param objSd: posterior standard deviation of the response
#' @details   compute the change in ZC = sd*(1-sqrt{(nugget)})/sqrt{(nugget + sd^2)}
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
#' cSUR for Contour Finding
#'
#' @title Compute reduction in contour-distance
#' @param objMean: predicted mean response
#' @param objSd: posterior standard deviation of the response
#' @param nugget the noise variance to compute the ALC factor
#' @details   compute the change in ZC = sd*(1-sqrt{(nugget)})/sqrt{(nugget + sd^2)}
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
#' @param objMean: predicted mean response
#' @param objSd: posterior standard deviation of the response
#' @details   compute normalized distance to zero-contour |mu|/sd
#' @export
cf.mcu <- function(objMean, objSd)
{
  return(pnorm(-abs(objMean)/objSd))
}

#####################
#' straddle MCU with a specified variance weight
#'
#' @title Straddle Maximum Contour Uncertainty criterion
#' @param objMean: predicted mean response
#' @param objSd: posterior standard deviation of the response
#' @param gamma: weight on the variance
#' @details   compute the UCB criterion with constant weight: gamma*s(x) - |f(x)|
#' @export
cf.smcu <- function(objMean, objSd, gamma=1.96)
{
  return( gamma*objSd - abs(objMean) )
}

#####################
#' tMSE for Contour Finding
#'
#' @title targeted Mean Squared Error criterion
#' @param objMean: predicted mean response
#' @param objSd: posterior standard deviation of the response
#' @param seps: epsilon in the tMSE formula. By default taken to be zero.
#' @details   compute predictive density at the contour, smoothed by seps
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

##############
swing.policy <- function( x,M,fit,model,offset=1,use.qv=FALSE,n.swing=1)
{
  nsim <- 0
  if (is(x, "matrix") | is(x, "numeric") ) {
    curX <- model$sim.func( x,model,model$dt)
    nsim <- nrow(x)
  }
  if (is(x, "list" ) )
    curX <- x[[1]]
  payoff <- rep(0, nrow(curX))
  tau <- rep(0, nrow(curX))

  contNdx <- 1:nrow(curX)
  i <- 1
  if (n.swing <= 0)
    return( list(payoff=payoff))
  ns <- rep(n.swing, nrow(curX))  # remaining rights for each path
  
  # main loop forward
  while (i < (M+(use.qv==TRUE)) & length(contNdx) > 0) {
    for (kk in 1:n.swing) {  # loop over number of remaining paths
      curNdx <- which( ns == kk)
      if (length(curNdx) == 0)
         continue
      # immediate payoff
      myFit <- fit[[i+1-offset,kk]]; myx <- curX[curNdx,,drop=F]
      imm  <- exp(-(i)*model$dt*model$r)*model$swing.payoff(myx , model$K)
      if (is(myFit,"earth") )
        rule <- predict(myFit,myx) # for use with  MARS
      if (is(myFit,"smooth.spline") )
        rule <- predict(myFit,myx)$y # for use with  splines
      if (is(myFit,"randomForest") ) {
        obj <- predict(myFit,myx,predict.all=T)$individual
        rule <- apply(obj,1,median)
      }
      if (is(myFit,"dynaTree")  ){
        obj <- predict(myFit, myx,quants=F)
        rule <- obj$mean
      }
      if (is(myFit,"km") ) # DiceKriging
        rule <- predict(myFit,data.frame(x=myx),type="UK")$mean
      if (is(myFit,"gpi") )  # laGP
        rule <- predGP(myFit,XX=myx, lite=TRUE)$mean
      if (is(myFit,"lm") ) {
        lenn <- length(myFit$coefficients)
        rule <-  myFit$coefficients[1] + model$bases(myx) %*% myFit$coefficients[2:lenn]
      }
      if( class(myFit)=="homGP" | class(myFit) == "hetGP" | class(myFit) == "homTP")
        rule <- predict(x=myx, object=myFit)$mean
      if( class(myFit)=="rvm")
        rule <-  predict(myFit, new=myx)
      if (class(myFit) == "npregression")
        rule <- predict(myFit, new=myx)
      
      if (use.qv == TRUE & i== M) {
        payoff[contNdx] = payoffCont[contNdx] + rule  # continuation value of paths that didn't stop yet)
        break
      }
      
      # stop if the expected gain is negative
      #contNdx <- contNdx[ which (rule > 0 | payoff[contNdx] == 0)]
      if (length(which(rule <0)) > 0) {
        stopKK <- which(rule <0)
        payoff[ curNdx[stopKK]] <- payoff[ curNdx[stopKK]] + imm[ stopKK]
        nS[ curNdx [stopKK]] <- nS[ curNdx[stopKK]] - 1
      }
    }
    tau[contNdx] <- (i)*model$dt
    
    if (compact == F) {
      sims[[min(i,model$look.ahead+1)]] <- curX[contNdx,,drop=F]
      save.ndx[[min(i,model$look.ahead+1)]] <- contNdx
    }
    # update the x values by taking a step of length dt
    i <- i+1
    
    if (is(x, "matrix") | is(x, "numeric") ) {
      curX[contNdx,] <- model$sim.func( curX[contNdx,,drop=F],model,model$dt)
      nsim <- nsim + length(contNdx)
    }
    
    if (is(x, "list") )  # stored list of paths
      curX[contNdx,] <- x[[i]][contNdx,,drop=F]
    # payoff for next timestep -- used for terminal payoff at i=M
    payoff[contNdx]  <- exp(-(i)*model$dt*model$r)*model$payoff.func( curX[contNdx,,drop=F], model)
  }
  for (i in 2:(model$look.ahead))   # payoff for a trajectory starting at x^n_{t+i} which was still alive then
    fvalue[[i]] <- payoff[ save.ndx[[i]] ]*exp((i-1)*model$dt*model$r)
  
  
  return( list(payoff=payoff,sims=sims,fvalue=fvalue,tau=tau, nsims=nsim))
  # payoff is the resulting payoff NPV from t=0
  # fvalue[i] is a list of resulting payoffs (on paths still not stopped) NPV from t=i
  # tau are the times when stopped
  # sims is a list; sims[[i]] are the forward x-values of paths at t=i (those not stopped yet)
}
