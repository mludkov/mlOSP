#####################
#' SUR for Contour Finding
#'
#' @title Calculates c_over in ABSUR
#' @param theta0, theta1 and theta2: parameters in linear regression
#' @param n: current design size
#' @details   compute the estimated overhead in ABSUR 
#' @export
CalcOverhead <- function(theta0, theta1, theta2, n) {
  overhead = theta0 + theta1*n + theta2*n^2;
}

#####################
#' ABSUR for Adaptive Batching
#'
#' @title Calculates weights for location and batch size in ABSUR
#' @param objMean: predicted mean response
#' @param objSd: posterior standard deviation of the response
#' @param nugget: the noise variance to compute the ALC factor
#' @param r_cand: candidate batch size vector
#' @param overhead: estimated computation overhead in GP
#' @param t0: overhead for individual simulation
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
#' RB for Adaptive Batching
#'
#' @title Calculates weights for batch size in RB
#' @param objSd_at_x_optim: posterior standard deviation of the response at the selected new input
#' @param r_cand: candidate batch size vector
#' @param nugget: the noise variance to compute the ALC factor
#' @param last_r: the last batch size
#' @param gamma: threshold compared with sd
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
#' MLB for Adaptive Batching
#'
#' @title Calculates weights for batch size in MLB
#' @param objSd_at_x_optim: posterior standard deviation of the response at the selected new input
#' @param r_cand: candidate batch size vector
#' @param nugget: the noise variance to compute the ALC factor
#' @param gamma: threshold compared with sd
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
#' @title Calculates reallocated batch size or new input location for ADSA
#' @param fit: gp/tp fit
#' @param r_seq: batch size vector for existing inputs
#' @param xtest: testing points to compare reallocation and adding a new inputs
#' @param xt_dens: density of xtest
#' @param x_new: new input location selected by one EI criteria
#' @param r0: total number of new simulations
#' @param nugget: the noise variance to compute the ALC factor
#' @param method: "km" or "trainkm" or "hetgp" or "homtp"
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
#' @param fit: gp/tp fit
#' @param r_seq: batch size vector for existing inputs
#' @param xtest: testing points to compare reallocation and adding a new inputs
#' @param xt_dens: density of xtest
#' @param r0: total number of new simulations
#' @param method: "km" or "trainkm" or "hetgp" or "homtp"
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
#' @param L: lower triangle of cholesky decomposition of covariance matrix
#' @param K: covariance matrix
#' @param xt_dens: density of xtest
#' @param r_seq: batch size vector for existing inputs
#' @param r0: total number of new simulations
#' @export
r.reallocate <- function(L, K, xt_dens, r_seq, r0) {
  U <- solve(t(L)) %*% (solve(L) %*% K) %*% xt_dens
  r_seq_new <- pegging.alg(r0 + sum(r_seq), U, r_seq)
  return(r_seq_new)
}

#####################
#' Pegging algorithm
#'
#' @title Calculates reallocated batch size
#' @param r: total number of simulations
#' @param U: weighted matrix for pegging algorithm
#' @param r_seq: batch size vector for existing inputs
#' @export
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
    r_new_int(idx_max) = r_new_int(idx_max) + 1
  }
  return(r_new_int)
}
