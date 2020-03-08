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


######################
#' Plots one-dimensional cuts of a 2-d fit from mc.put.adaptree
#' @export
plt.2d.proj <- function( fit, x=seq(30,43,len=201),y=40,method="dt") {
  gr2 <- expand.grid(x=x,y=y)
  if (is(fit,"randomForest") ) {
       obj <- predict(fit,cbind(gr2$x,gr2$y),predict.all=T)$individual
       obj <- apply(obj,1,median)
  }
  if (is(fit,"dynaTree")  )
      obj <- predict(fit,cbind(gr2$x,gr2$y),quants=F)$mean

  plot(gr2$x,obj,type="l",lwd=1.5,xlab=expression(X[t]),ylab="Timing Value")
  abline(h=0,lty=2,lwd=0.5)
  # switch the order around between X and Y
  if (is(fit,"randomForest") ) {
       obj <- predict(fit,cbind(gr2$y,gr2$x),predict.all=T)$individual
       obj <- apply(obj,1,median)
  }
  if (is(fit,"dynaTree")  )
      obj <- predict(fit,cbind(gr2$y,gr2$x),quants=F)$mean
  lines(gr2$x,obj,col="red",lwd=1.5)
}

######
#' two-dimensional image+contour plot
#' @title Visualize 2D emulator + stopping region
#' @details Uses the quilt.plot command from \pkg{fields}
#'
#' @param x,y locations to use for the \code{predict()} functions. Default is a 200x200 fine grid.
#' Passed to \code{expand.grid}
#' @param fit can be any of the types supported by \code{\link{forward.sim.policy}}
#' @param show.var -- if \code{TRUE} then plot posterior surrogate variance instead of surrogate mean [default = FALSE]
#' This only works for \code{km} and \code{het/homGP} objects
#' @param only.contour -- just add the zero-contour, no quilt.plot
#' @param ub to create an upper bound to see the zero-contour better
#' @param .. -- pass additional options to \code{quilt.plot}
#' @param contour.col (default is "red") -- color of the zero contour
#' @export
plt.2d.surf <- function( fit, x=seq(31,43,len=201),y = seq(31,43,len=201),ub=1e6,
    show.var=FALSE, only.contour=FALSE, contour.col="red",...)
{
   gr <- expand.grid(x=x,y=y)

   if (is(fit,"randomForest") ) {
       obj <- predict(fit,cbind(gr$x,gr$y),predict.all=T)$individual
       obj <- apply(obj,1,median)
   }
   if (is(fit,"dynaTree")  )
      obj <- predict(fit,cbind(gr$x,gr$y),quants=F)$mean

   if (class(fit)=="km" & show.var == FALSE)
      obj <- predict(fit,data.frame(x=cbind(gr$x,gr$y)), type="UK")$mean
   if (is(fit,"km") & show.var == TRUE)
      obj <- predict(fit,data.frame(x=cbind(gr$x,gr$y)), type="UK")$sd
   if( (class(fit)=="homGP" | class(fit) == "hetGP") & show.var == FALSE)
      obj <- predict(x=cbind(gr$x,gr$y), object=fit)$mean
   if( (class(fit)=="homGP" | class(fit) == "hetGP") & show.var == TRUE)
     obj <- sqrt(predict(x=cbind(gr$x,gr$y), object=fit)$sd2)
   if( class(fit)=="rvm")
     obj <-  predict(fit, new=cbind(gr$x,gr$y))


  if (only.contour== FALSE)
    quilt.plot(gr$x, gr$y, pmin(ub,obj),xlim=range(x),ylim=range(y),
      xlab=expression(X[t]^1), ylab=expression(X[t]^2),cex.lab=1.2, cex.axis=1.1,...)
  imx <- as.image(x=cbind(gr$x,gr$y),obj,nr=100,nc=100)
  contour(imx$x,imx$y,imx$z,levels=0,add=T,drawlab=F,lwd=2,col=contour.col)
  if (class(fit)=="dynaTree") { #is.null(fit@X) == F & is(fit,"km")==F) {
    kk <- kde2d(fit$X[,1],fit$X[,2],lims=c(range(x),range(y)))
    contour(kk$x, kk$y, kk$z,nlev=12,lwd=0.5,lty=2,drawlab=F,add=T,col=contour.col)
  }
  if (is(fit,"km") & (only.contour== FALSE))
      points( fit@X, col="orange", pch=19)
  if( (class(fit)=="hetGP" | class(fit)=="homGP") & (only.contour== FALSE))
      points(fit$X0, col="orange", pch=19)
}

##################################
#' Show the final fits for each timestep
#' @param stop.tree is a list of dynatrees
#' @param x.gr is the grid on which to show the fit
#' @export
##################################
plt.fits <- function(stop.tree, len2, x.gr = seq(32,40,len=801))
{
   for (i in 1:len) {
       pred <- predict(stop.tree[[i]],as.matrix(x.gr) )$mean
       plot( x.gr, pred,
         col=i, ylim=c(-0.25,0.5),type="l",ylab="Timing value", xlab=paste("Time step", i))
       abline(h=0)
       readline()
   }
}



##########################
#' Visualize the fit of Mike's trees
#'
plt.fit.tree <- function(fit, model,grids=NULL, timingValue=NULL)
{
    if (!is.null(timingValue))
      plot(grids,timingValue,xlim=c(0.8*model$K,model$K+1),ylim=c(-1,2),xlab='x', ylab='y',cex=0.45)
    else
      plot(0,0,xlim=c(0.65*model$K,model$K+1),ylim=c(-2,2),xlab='Stock Price', ylab='Timing Value')

    for (jj in 1:length(fit[[1]]))  {
          lines(c(fit[['lower']][jj], fit[['upper']][jj]),c(fit[['mean']][jj],
          fit[['mean']][jj]),col="red",lwd=2.5);
          if (fit[['mean']][jj] < 0)
             lines(c(fit[['lower']][jj], fit[['upper']][jj]),c(-2,-2),col="blue",lwd=3)

    }
    abline(h=model$target.contour,col="black")
}

##########################
#' Visually compare the BW and tree 1-d fits using same grid
#' @param t.ndx is the time index
#' @param ls.fit is from mc.put.ls (only in 1-d!!)
#' @param tree.fit from mc.put.adaptree
#' @param ... is for the parameters to paths to the plots, such as axis labels
#'
comp.plot <- function(t.ndx,x.grid,ls.fit,tree.fit,model,lims=c(31,39),...)
{
    to.plot <- which( x.grid[[t.ndx]] > lims[1] & x.grid[[t.ndx]] < lims[2])
    plot(x.grid[[t.ndx]][to.plot], ls.fit[[t.ndx]][to.plot]-option.payoff(x.grid[[t.ndx]][to.plot,,drop=F],model$K),
       xlim=lims,ylim=c(-0.25,0.5),xlab="", ylab='Timing value',col="blue",...)
    x.gr <- seq(lims[1],lims[2],len=801)
    lines( x.gr, predict(tree.fit[[t.ndx]], x.gr)$mean,col="red",lwd=2)
    abline(h=0)
    legend("topleft",c('DT','BW11'),col=c('red', 'blue'), lty=1,lwd=2)
    text(x=30,y=0.1,labels='Stopping Region')
    hist(tree.fit[[t.ndx]]$X,seq(20,41,by=0.5),xlab="X_t", main="", col=rgb(1,0,0,0.3),freq=F,xlim=lims)
    hist(x.grid[[t.ndx]],seq(17,91,by=0.5),col=rgb(0,0,1,0.3),freq=F,add=T)
}

