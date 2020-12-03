#####################
#' ABSUR overhead
#'
#' @title Calculates c_over in ABSUR
#' @param theta0,theta1,theta2 parameters in linear regression
#' @param n current design size
#' @details   compute the estimated overhead in ABSUR
#'
CalcOverhead <- function(theta0, theta1, theta2, n) {
  overhead = theta0 + theta1*n + theta2*n^2;
}

#####################
#' ABSUR for Adaptive Batching
#'
#' @title Calculates weights for location and batch size in ABSUR
#' @param objMean predicted mean response
#' @param objSd posterior standard deviation of the response
#' @param nugget the noise variance to compute the ALC factor
#' @param r_cand candidate batch size vector
#' @param overhead estimated computation overhead in GP
#' @param t0 overhead for individual simulation
#' @references 
#' M. Ludkovski, X. Lyu (2020+) Adaptive Batching for Gaussian Process Surrogates with Application 
#' in Noisy Level Set Estimation, <http://arxiv.org/abs/2003.08579>
#' @author Xiong Lyu
#' @seealso \code{\link[mlOSP]{osp.seq.batch.design}}
#' @export
cf.absur <- function(objMean, objSd, nugget, r_cand, overhead, t0) {
  # expand mean and sd vectors to matrix of size len(x_cand) * len(r_cand)
  r_len <- length(r_cand)
  x_len <- length(objMean)
  objMean_matrix <- matrix(objMean, nrow=x_len, ncol=r_len)
  objSd_matrix <- matrix(objSd, nrow=x_len, ncol=r_len)
  r_matrix <- matrix(r_cand, nrow = x_len, ncol=r_len, byrow=TRUE)

  # EI at cross combination of candidate input and batch size
  nugget_matrix <- nugget / sqrt(r_matrix)
  a = pnorm(-abs(objMean_matrix)/objSd_matrix)  # normalized distance to zero contour
  new_objSd2 <- nugget_matrix * objSd_matrix / sqrt(nugget_matrix ^ 2 + objSd_matrix ^ 2) # new posterior variance
  a_new <- pnorm(-abs(objMean_matrix)/new_objSd2)   # new distance to zero-contour

  # difference between next-step ZC and current ZC weighted by the overhead
  return( (a-a_new) / (r_matrix * t0 + overhead) )
}

#####################
#' Ratchet Batching heuristic
#'
#' @title Calculates weights for batch size in RB
#' @param objSd_at_x_optim posterior standard deviation of the response at the selected new input
#' @param r_cand candidate batch size vector
#' @param nugget the noise variance to compute the ALC factor
#' @param last_r the last batch size
#' @param gamma threshold compared with sd
#' @seealso \code{\link[mlOSP]{osp.seq.batch.design}}
#' @references 
#' M. Ludkovski, X. Lyu (2020+) Adaptive Batching for Gaussian Process Surrogates with Application 
#' in Noisy Level Set Estimation, http://arxiv.org/abs/2003.08579
#' @export
batch.rb <- function(objSd_at_x_optim, r_cand, last_r, nugget, gamma) {
  eta = 0.8
  rindex = which(r_cand == last_r)[1]

  # if r reaches the upper bound, just keep it at the same level
  if (rindex == length(r_cand)) {
    roptim = r_cand[rindex]
  } else {

    # compare whether to stay at the current r or move to the next level
    rcand_two_levels = c(r_cand[rindex], r_cand[rindex + 1])
    nugget2_est_by_level = nugget ^ 2 / rcand_two_levels
    post_var = objSd_at_x_optim^2 * nugget2_est_by_level / (nugget2_est_by_level + objSd_at_x_optim^2)
    post_sd = sqrt(post_var)

    # if no solution, update value of gamma
    if (post_sd[1] < gamma) {
      k = ceiling(log(post_sd[1] / gamma) / log(eta))
      gamma = gamma * (eta ^ k)
    }

    # choose whether to stay depending on comparison of sd with gamma
    r_index = rindex - 1 + sum(post_sd > gamma)
    roptim = r_cand[r_index]
  }
  return(list(roptim = roptim, gamma = gamma))
}

#####################
#' Multi-Level Batching Heuristic
#'
#' @description Calculates weights for batch size in MLB when called from \link{osp.seq.batch.design}
#' @param objSd_at_x_optim posterior standard deviation of the response at the selected new input
#' @param r_cand candidate batch size vector
#' @param nugget the noise variance to compute the ALC factor
#' @param gamma threshold compared with sd
#' @return  list containing:
#'  \itemize{ 
#'  \item \code{roptim}: new replication count
#'  \item \code{gamma}: new gamma variable
#'  }
#' @references 
#' M. Ludkovski, X. Lyu (2020+) Adaptive Batching for Gaussian Process Surrogates with Application 
#' in Noisy Level Set Estimation, http://arxiv.org/abs/2003.08579
#' 
#' @seealso \code{\link[mlOSP]{osp.seq.batch.design}}
#' @export
batch.mlb <- function(objSd_at_x_optim, r_cand,  nugget, gamma) {
  eta = 0.5;

  # calculate the new posterior variance at the new site location
  nugget2_est_by_level = nugget ^ 2 / r_cand
  post_var = objSd_at_x_optim^2 * nugget2_est_by_level / (nugget2_est_by_level + objSd_at_x_optim^2)
  post_sd = sqrt(post_var)

  # if no solution, update value of gamma
  if (post_sd[1] < gamma) {
    k = ceiling(log(post_sd[1] / gamma) / log(eta))
    gamma = gamma * (eta ^ k)
  }

  r_index = sum(post_sd > gamma)
  roptim = r_cand[r_index]
  return(list(roptim = roptim, gamma = gamma))
}

#####################
#' ADSA for Adaptive Batching
#'
#' @title Calculates reallocated batch size or new input location for Adaptive Design with Sequential Allocation
#' @param fit GP/TP fit
#' @param r_seq batch size vector for existing inputs
#' @param xtest testing points to compare reallocation and adding a new inputs
#' @param xt_dens density of xtest
#' @param x_new new input location selected by the EI criteria
#' @param r0 total number of new simulations
#' @param nugget the noise variance to compute the ALC factor
#' @param method \code{km} or \code{trainkm} or \code{hetgp} or \code{homtp}
#' @return a list containing
#' \itemize{
#' \item \code{xoptim}: new design input (NULL if re-allocation is chosen)
#' \item \code{roptim}: added replications (a scalar r0 if new input chosen, 
#' a vector containing the re-allocation amounts otherwise
#' }
#' @references 
#' M. Ludkovski, X. Lyu (2020+) Adaptive Batching for Gaussian Process Surrogates with Application 
#' in Noisy Level Set Estimation, <http://arxiv.org/abs/2003.08579>
#' @author Xiong Lyu
#' @seealso \code{\link[mlOSP]{osp.seq.batch.design}}
#' @export
batch.adsa <- function(fit, r_seq, xtest, xt_dens, x_new, r0, nugget, method) {
  x_optim = NULL
  ddsa.res = batch.ddsa(fit, r_seq, xtest, xt_dens, r0, method)
  r_new = ddsa.res$r_new
  K = ddsa.res$K
  L = ddsa.res$L
  # determine between allocating in existing samples or moving to a new location

  # i1 - reallocation
  Delta_R = 1/r_seq - 1/r_new
  Delta_R = diag(Delta_R)

  v = solve(t(L)) %*% (solve(L) %*% K)
  I1 = diag(t(v) %*% Delta_R %*% v)

  # i2 - adding a new input
  if (method == "km" | method == "trainkm") {
    x_all = rbind(xtest, fit@X, x_new)
    C = covMatrix(fit@covariance, x_all)$C
    k = C[1:nrow(xtest), (nrow(xtest) + nrow(fit@X) + 1)]
    k_new = C[(nrow(xtest) + 1):(nrow(xtest) + nrow(fit@X)), (nrow(xtest) + nrow(fit@X) + 1)]
  } else {
    x_all = rbind(xtest, fit$X0, x_new)

    if (method == "homtp") {
      # Lower triangle matrix of covariance matrix
      C <- hetGP::cov_gen(x_all, theta=fit$theta, type="Gaussian") * fit$sigma2
    } else {
      C <- hetGP::cov_gen(x_all, theta=fit$theta, type="Gaussian") * fit$nu_hat
    }
    k = C[1:nrow(xtest), (nrow(xtest) + nrow(fit$X0) + 1)]
    k_new = C[(nrow(xtest) + 1):(nrow(xtest) + nrow(fit$X0)), (nrow(xtest) + nrow(fit$X0) + 1)]
  }

  if (method == "km" | method == "trainkm") {
    ss = fit@covariance@sd2
  } else {
    if (method == 'homtp') {
      ss = fit$sigma2
    } else {
      ss = fit$nu_hat
    }
  }

  a = solve(t(L)) %*% (solve(L) %*% k_new)
  var_new = ss - t(k_new) %*% a
  cov = k - t(K) %*% a
  I2 = vec(cov) ^ 2 / (nugget^2 / r0 + var_new[1, 1])

  if (sum(I1 * xt_dens) > sum(I2 * xt_dens)) {
    # allocate on existing samples
    r_optim = r_new
  } else {
    # choose a new location
    r_optim = r0
    x_optim = x_new
  }
  return(list(x_optim = x_optim, r_optim = r_optim))
}

#####################
#' DDSA for Adaptive Batching
#'
#' @title Calculates reallocated batch size for DDSA
#' @param fit gp/tp fit
#' @param r_seq batch size vector for existing inputs
#' @param xtest testing points to compare reallocation and adding a new inputs
#' @param xt_dens density of xtest
#' @param r0 total number of new simulations to add
#' @param method "km" or "trainkm" or "hetgp" or "homtp"
#' @references 
#' M. Ludkovski, X. Lyu (2020+) Adaptive Batching for Gaussian Process Surrogates with Application 
#' in Noisy Level Set Estimation, <http://arxiv.org/abs/2003.08579>
#' @author Xiong Lyu
#' @seealso \code{\link[mlOSP]{osp.seq.batch.design}}
#' @export
batch.ddsa <- function(fit, r_seq, xtest, xt_dens, r0, method) {
  if (method == "km" | method == "trainkm") {
    pred.test <- predict(fit, data.frame(x=xtest), type="UK")

    # Covariance of xtest and x
    K <- pred.test$c
    # Lower triangle matrix of covariance matrix
    L <- t(fit@T)
  } else {   # hetGP/homTP
    pred.test <- predict(x=xtest, object=fit, xprime = fit$X0)

    # Covariance of xtest and x
    K <- t(pred.test$cov)

    if (method == "homtp") {
      # Lower triangle matrix of covariance matrix
      C <- hetGP::cov_gen(fit$X0, theta=fit$theta, type="Gaussian") * fit$sigma2 + diag(rep(fit$g, dim(fit$X0)[1]))
    } else {
      C <- fit$nu_hat * (hetGP::cov_gen(fit$X0, theta=fit$theta, type="Gaussian") + fit$Lambda * diag(1/fit$mult))
    }
    L <- t(chol(C))
  }


  # allocate between existing samples
  r_new = r.reallocate(L, K, xt_dens, r_seq, r0)
  return(list(r_new = r_new, K = K, L = L))
}

#####################
#' New batch size calculator
#'
#' @title Calculates reallocated batch size
#' @param L lower triangle of cholesky decomposition of covariance matrix
#' @param K covariance matrix
#' @param xt_dens density of xtest
#' @param r_seq batch size vector for existing inputs
#' @param r0 total number of new simulations
r.reallocate <- function(L, K, xt_dens, r_seq, r0) {
  U <- solve(t(L)) %*% (solve(L) %*% K) %*% xt_dens
  r_seq_new <- pegging.alg(r0 + sum(r_seq), U, r_seq)
  return(r_seq_new)
}

#####################
#' Pegging algorithm for ADSA/DDSA
#'
#' @title Calculates reallocated batch size for ADSA
#' @param r total number of simulations
#' @param U weighted matrix for pegging algorithm
#' @param r_seq batch size vector for existing inputs
#' @references 
#' M. Ludkovski, X. Lyu (2020+) Adaptive Batching for Gaussian Process Surrogates with Application 
#' in Noisy Level Set Estimation, http://arxiv.org/abs/2003.08579
#' @author Xiong Lyu
pegging.alg <- function(r, U, r_seq) {
  is_end = FALSE
  indexes = seq(1, length(r_seq))
  r_new = r_seq
  r_total = r

  while(!is_end) {
    r_new[indexes] = r_total * U[indexes] / sum(U[indexes])

    if(sum(r_new[indexes] >= r_seq[indexes]) == length(indexes)) {
      is_end = TRUE
    } else {
      idx = which(r_new <= r_seq)
      unchanged = intersect(idx, indexes)
      r_new[unchanged] = r_seq[unchanged]
      indexes = which(r_new > r_seq)
      r_total = r - sum(r_new[idx])
    }
  }
  r_new_int = round(r_new)
  if (sum(r_new_int == r_seq) == length(r_seq)) {
    idx_max = which(r_new == max(r_new[r_new != r_new_int]))
    r_new_int[idx_max] = r_new_int[idx_max] + 1
  }
  return(r_new_int)
}


######
#' two-dimensional image of contour + site + batch plot for two fits
#' @title Visualize and compare 2D emulator + stopping region for two fits
#'
#' @param fit1,fit2 can be any of the types supported by \code{\link{forward.sim.policy}}
#' @param r1,r2 batch vectors for the two fits
#' @param x,y locations to use for the \code{predict()} functions. Default is a 200x200 fine grid.
#' Passed to \code{expand.grid}
#' @param batch1,batch2 batch heristics for two fits; Passed to \code{ggplot()}
#' This only works for \code{km} and \code{het/homGP} objects
#' @export
plt.2d.surf.batch <- function( fit1, fit2, r1, r2, batch1, batch2, x=seq(31,43,len=201), y = seq(31,43,len=201))
{
  gr <- expand.grid(x=x,y=y)
  r <- c(r1, r2)
  batch.samples <- c(rep(batch1, length(r1)), rep(batch2, length(r2)))
  batch.fitted <- c(rep(batch1, length(x) ^ 2), rep(batch2, length(x) ^ 2))

  # calculate posterior mean and standard deviation for grid
  if (class(fit1)=="km") {
    x <- rbind(fit1@X, fit2@X)
    m <- c(predict(fit1,data.frame(x=cbind(gr$x,gr$y)), type="UK")$mean,
           predict(fit2,data.frame(x=cbind(gr$x,gr$y)), type="UK")$mean)
    sd<- c(predict(fit1,data.frame(x=cbind(gr$x,gr$y)), type="UK")$sd,
           predict(fit2,data.frame(x=cbind(gr$x,gr$y)), type="UK")$sd)
  }
  if( (class(fit1)=="homGP" | class(fit1) == "hetGP")) {
    x <- rbind(fit1$X0, fit2$X0)
    m <- c(predict(x=cbind(gr$x,gr$y), object=fit1)$mean,
           predict(x=cbind(gr$x,gr$y), object=fit2)$mean)
    sd <- sqrt(c(predict(x=cbind(gr$x,gr$y), object=fit1)$sd2,
                 predict(x=cbind(gr$x,gr$y), object=fit2)$sd2))
  }

  # credible interval
  lb <- m - 1.96 * sd
  ub <- m + 1.96 * sd

  fitted.data.2d <- data.frame(x1 = gr$x, x2 = gr$y, m = m, lbound = lb, ubound = ub, batch = batch.fitted)
  samples <- data.frame(x1 = x[,1], x2 = x[,2], r = r, batch = batch.samples)

  p <- ggplot(data = samples, aes(x = x1, y = x2)) +
    geom_raster(data = fitted.data.2d, aes(x = x1, y = x2, fill = m)) +
    scale_fill_gradient2(low = "black", high = "yellow", mid = "red") +
    geom_point(aes(color = r), size = 2) +
    scale_color_gradient(low = "cyan", high = "purple") +
    facet_grid(col = vars(batch)) + labs(x=expression(X[t]^1), y= expression(X[t]^2)) +
    scale_x_continuous(name = expression(X[t]^1),  expand = c(0, 0)) +
    scale_y_continuous(name = expression(X[t]^2),  expand = c(0, 0))

  p + stat_contour(data = fitted.data.2d, aes(x = x1, y = x2, z = m), breaks = 0, color = "black", inherit.aes=FALSE, size = 1.4) +
    facet_grid(col = vars(batch)) +
    stat_contour(data = fitted.data.2d, aes(x = x1, y = x2, z = lbound), color = "black", linetype = "longdash", breaks = 0, size = 0.8) +
    stat_contour(data = fitted.data.2d, aes(x = x1, y = x2, z = ubound), color = "black",linetype = "longdash", breaks = 0, size = 0.8) +
    plot_style()
}

##########################
#' GGplot style
#' @param base_size is the font size
#' @param base_family is the font style
#' @param ... is for the parameters to specialize the ggplot style
#' @noMd
plot_style <- function(base_size = 14, base_family = "Helvetica",...) {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_blank(),
          axis.line = element_line(colour = "black",size = 0.5),
          strip.background = element_rect(fill = "white",colour = "white",size = 1),
          legend.direction   = "vertical",
          legend.position = "right",
          legend.box = "vertical",
          legend.text.align = 0,
          legend.text = element_text(size=8, margin = margin(r = 5, unit = "pt")),
          legend.box.just = "left",
          legend.key = element_blank(),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.margin = margin(c(2, 2, 2, 0)),
          legend.key.width = unit(0.35,"cm"),
          legend.key.height = unit(0.7,"cm"), 
          axis.text.x = element_text( colour = 'black', size = 10, hjust = 0.5, vjust = 0.5),
          axis.title.x = element_text(size = 10, hjust = 0.5, vjust = 0.2),
          axis.text.y = element_text(colour = 'black', size = 10),
          axis.title.y = element_text(size = 10, angle = 90, hjust = 0.5, vjust = 0.2),
          strip.text.y = element_text(size = 10, hjust = 0.5,  vjust = 0.5, face = 'bold'),
          strip.text.x = element_text(size = 10, hjust = 0.5,  vjust = 0.5, face = 'bold'),
          ...)
}

######
#' two-dimensional image+contour plot with replication counts for an \code{osp.seq.batch.design} fit
#' @title Visualize 2D emulator + stopping region + batch amounts 
#'
#' @param x,y locations to use for the \code{predict()} functions. Default is a 200x200 fine grid.
#' Passed to \code{expand.grid}
#' @param fit An emulator object. Can be any of the types supported by \code{\link[mlOSP]{forward.sim.policy}}
#' @param show.var if \code{TRUE} then plot posterior surrogate variance instead of 
#' surrogate mean (Default = FALSE)
#' This only works for \code{km} and \code{het/homGP} objects
#' @param only.contour -- just add the zero-contour, no geom_raster
#' @param contour.col (default is "red") -- color of the zero contour
#' @param batch_size matrix of replication counts for each input generated by \code{\link[mlOSP]{osp.seq.batch.design}}
#' 
#' @author Xiong Lyu
#' @export
plt.2d.surf.with.batch <- function( fit, batch_size, x=seq(25,50,len=201),y = seq(25,50,len=201),
                         show.var=FALSE, only.contour=FALSE, contour.col="red")
{
  gr <- expand.grid(x=x,y=y)

  # calculate posterior mean and standard deviation for grid
  if (class(fit)=="km") {
    m <- predict(fit,data.frame(x=cbind(gr$x,gr$y)), type="UK")$mean
    sd<- predict(fit,data.frame(x=cbind(gr$x,gr$y)), type="UK")$sd
    samples <- data.frame(x1 = fit@X[,1], x2 = fit@X[,2], r = batch_size)
  }
  if( (class(fit)=="homGP" | class(fit) == "hetGP")) {
    m <- predict(x=cbind(gr$x,gr$y), object=fit)$mean
    sd <- sqrt(predict(x=cbind(gr$x,gr$y), object=fit)$sd2)
    samples <- data.frame(x1 = fit$X0[,1], x2 = fit$X0[,2], r = batch_size)
  }

  fitted.data.2d <- data.frame(x1 = gr$x, x2 = gr$y, m = m, 
                               lbound = m - 1.96 * sd , ubound = m + 1.96 * sd)
  
  cols <- RColorBrewer::brewer.pal(n = 9, name = "PuBuGn")

  p <- ggplot(fitted.data.2d) +
    geom_raster(aes(x = x1, y = x2, fill = m)) +
    #scale_fill_gradient2(low = "black", high = "yellow", mid = "red") +
    scale_fill_gradientn(colours = fields::tim.colors(64)) +
    geom_point(data=samples, aes(x = x1, y = x2, color = r), size = 2) +
    # scale_color_gradient(low = "cyan", high = "blue", limits = c(0, 200)) +
    scale_colour_gradientn(colours = cols, values = scales::rescale(c(seq(20, 140, 8), 188)), guide = "colorbar", limits=c(10, 188)) +
    scale_x_continuous(name = expression(X[t]^1),  expand = c(0, 0),limits=range(x)) +
    scale_y_continuous(name = expression(X[t]^2), expand = c(0, 0),limits=range(y))

  p + stat_contour(data = fitted.data.2d, aes(x = x1, y = x2, z = m), breaks = 0, color = contour.col, inherit.aes=FALSE, size = 1.4) +
    stat_contour(data = fitted.data.2d, aes(x = x1, y = x2, z = lbound), color = "black", linetype = "longdash", breaks = 0, size = 0.8) +
    stat_contour(data = fitted.data.2d, aes(x = x1, y = x2, z = ubound), color = "black",linetype = "longdash", breaks = 0, size = 0.8) +
    plot_style()
  
}
