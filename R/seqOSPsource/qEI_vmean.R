#################
#' EI using Ranjan, Bingham and Michailidis
qEI.vmean <- function(object, XX=NULL, q=0, alpha=2, verb=0)
  {
    ## sanity check
    if(object$model == "class") stop("only used by regression models")

  	## check XX argument
  	if(!is.null(XX)) object <- predict(object, XX, quants=FALSE, verb=verb)
  	else if(is.null(object$mean)) stop("must call predict first, or specify XX")

    if (is.null(object$alc))
       object <- alc(object,XX)

  	## for timing purposes
    p1 <- proc.time()[3]

    ## check q, alpha
    if(length(q) != 1) stop("q must be a scalar")
    if(length(alpha != 1) && alpha <= 0) stop("alpha must be a positive scalar")

    ## utility functions
    sd <- sqrt(object$vmean)
    sd <- sqrt(object$alc*mean(object$vmean)/mean(object$alc) )
    mean <- object$mean
    epsilon = alpha * sd;
    u1 <- (q - mean - epsilon)/sd;
    u2 <- (q - mean + epsilon)/sd;

    ## calculate quantile expected improvement
    qei <- (epsilon^2 - (mean - q)^2 - sd^2)*(pnorm(u2) - pnorm(u1));
    qei <- qei + (sd^2)*(u2*dnorm(u2) - u1*dnorm(u1));
    qei <- qei + 2.0*(mean - q)*sd*(dnorm(u2) - dnorm(u1));

    ## check for nans, etc
    qei[!is.finite(qei)] <- 0
    qei <- pmax(qei,0)

    ## put in object
    object$qEI <- qei

    ## update time
    object$time <- object$time + proc.time()[3] - p1

    ## done
    return(object)
  }

##############
#' Wrapper on \code{\link{qEI.vMean}} to guarantee positive weights
#' @export
  qEI.rbm <- function(obj,alph=2)
  {
      ei.alph <- alph
      obj <- qEI.vmean(obj, alph=ei.alph)
      while (max(obj$qEI) < 1e-9)  {
              ei.alph <- ei.alph + 1
              obj$vmean <- obj$var/(obj$df+1)
              obj <- qEI.vmean(obj,alph=ei.alph)
      }
      return (obj$qEI)
  }

###############
#' Expected Improvement for Contour Finding
#'
#' Compute EI using the probability of getting the sign wrong
#' @param method Possible values are 'vmean'.
#' @keywords internal
#' @examples
#' qEI.sgn()
#' @export
  qEI.sgn <- function(obj,anneal.pow=1,method="vmean")
  {
      qei <- obj$var/(obj$df+1)
      if (method=="vmean" & diff(range(obj$vmean)) > 1e-9)
         #qei <- obj$vmean
         qei <- obj$alc*mean(obj$vmean)/mean(obj$alc)

      qei <- pnorm( -abs(obj$mean)/sqrt(qei) )
      qei <- qei^anneal.pow
      return (qei)
  }

###############
#' Expected Improvement for Contour Finding
#'
#'  Compute EI using the optimal stopping loss function.
#' @param obj  must be a DynaTree or a list with an sd field
#' @export
  qEI.zc <- function(obj,method="posterior")
  {
      if (is(obj, "dynaTree"))
        sd <- sqrt(obj$var/(obj$df+1))
      if (is(obj, "list"))
         sd <- obj$sd

      qei <- pmax(0, sd*dnorm( -abs(obj$mean)/sd ) - abs(obj$mean)*pnorm( -abs(obj$mean)/sd) )
      return (qei)
  }

  ###############
  #' Expected Loss for Contour Finding
  #'
  #'  Compute expected loss using the optimal stopping loss function.
  #' @param obj  must be a DynaTree or a list with an sd field
  #' @export
  cf.el <- function(objMean,objSd)
  {
    el <- pmax(0, objSd*dnorm( -abs(objMean)/objSd ) - abs(objMean)*pnorm( -abs(objMean)/objSd) )
    return (el)
  }


##########
#' Expected Improvement for Contour Finding
#'
#' Compute EI using the UCB on ALC where the weight is gamma
#' @param obj  must be a DynaTree
#' @param gamma is the UCB weight balancing alc and the mean
#' @export
  qEI.alc <- function(obj,gamma=100)
  {
      if (is(obj, "dynaTree"))
        alc <- obj$alc

      qei <- -abs(obj$mean) + gamma*alc
      return (qei)
  }


###################
#' Compute EI using the squared mean
#' @export
 qEI.sq <- function(obj,anneal=1)
 {
    return( exp(- al.fit[[i]]$mean^2/al.fit$var*anneal) )
 }


#####################
#' Expected Improvement for Contour Finding
#'
#' Compute EI using the SUR formula
#' @param obj  must have mean and sd fields
#' @param nugget is the noise variance to compute the ALC factor
#'    alc = sd*(1-sqrt{(nugget)})/sqrt{(nugget + sd)}
#' @export
   qEI.sur <- function(tmp,nugget)
   {
         #plotinfo[[stp]] <<- cbind(X = as.vector(xs), fit = tmp_o$mean, se = tmp_o$sd, lower = tmp_o$lower95, upper = tmp_o$upper95)
         a = abs(tmp$mean)/tmp$sd
         M = tmp$mean*pnorm(-a)-tmp$sd*dnorm(-a)
         alc <- tmp$sd*(1-sqrt(nugget)/sqrt(nugget+tmp$sd))
         newSigma <- (tmp$sd - alc)
         a_new <- abs(tmp$mean)/newSigma
         M_new = tmp$mean*pnorm(-a_new)-newSigma*dnorm(-a_new)
         return( M_new-M)
   }

   #####################
   #' SUR for Contour Finding
   #'
   #' @title Compute EI for Contour Finding using the ZC-SUR formula
   #' @param objMean: predicted mean response
   #' @param objSd: posterior standard deviation of the response
   #' @param nugget the noise variance to compute the ALC factor
   #' @details   compute the change in ZC = sd*(1-sqrt{(nugget)})/sqrt{(nugget + sd)
   #' @export
   cf.sur <- function(objMean, objSd, nugget)
   {
     #plotinfo[[stp]] <<- cbind(X = as.vector(xs), fit = tmp_o$mean, se = tmp_o$sd, lower = tmp_o$lower95, upper = tmp_o$upper95)
     a = abs(objMean)/objSd  # normalized distance to zero contour
     M = objMean*pnorm(-a)-objSd*dnorm(-a)
     var_reduce <- objSd*(1-sqrt(nugget)/sqrt(nugget+objSd))  # reduction in posterior stdev from a new observation
     newSigma <- (objSd - var_reduce)   # look-ahead variance
     a_new <- abs(objMean)/newSigma   # new distance to zero-contour
     M_new = objMean*pnorm(-a_new)-newSigma*dnorm(-a_new)  # new ZC measure
     # difference between next-step ZC and current ZC
     return( M_new-M)
   }


######################
#' Expected Improvement for Contour Finding
#'
#'  Compute EI using the UCB formula on the OST loss function
#' @param tmp  must have mean/sd fields
#' @param gamma is the UCB weight balancing alc and the mean
#'
#' @export
   qEI.ucb <- function(tmp, gamma)
   {
         #plotinfo[[stp]] <<- cbind(X = as.vector(xs), fit = tmp_o$mean, se = tmp_o$sd, lower = tmp_o$lower95, upper = tmp_o$upper95)
         a = abs(tmp$mean)/tmp$sd
         M = tmp$mean*pnorm(-a)-tmp$sd*dnorm(-a)
         return( -M + gamma*tmp$sd)
   }
