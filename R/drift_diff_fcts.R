
require(Rcpp)
require(RcppArmadillo)
require(mvtnorm)
require(gtools)


#' Half Cauchy pdf
#' 
#' @param x points where to evaluate the pdf
#' @param sigma scale parameter
#' @return pdf for the half Cauchy distribution with location 0 and scale $\sigma$
dhalfcauhy <- function(x, sigma, log = T){
  
  out <- log(2) + log(sigma) - log(pi) - log(sigma^2 + x^2)
  if (!log){
    out <- exp(out)
  }
  
  return(out)
}


#' Construct the covariance matrix P of the smoothness inducing prior for the
#' spline coefficients
#' 
#' @param K number of spline knots
#' @return covariance of the smoothness inducing prior (penalizing first
#' differences in beta)
P_smooth1 <- function(K){
  D <- diag(rep(1,K))
  D <- diff(D)
  P <- crossprod(D)
  
  return(P)
}


#' Computes the g() function necessary for locally informed moves.
#' Details are given in Section S.4 in the supplementary materials
#' 
#' @param x input of the function
#' @param log should the result be computed on a log scale?
#' @return g(x)
g_HB <- function(x, log = T){
  # In the case we do not want to use the locally informed version, 
  # set out = 0
  out <- x/2
  
  if (log){
    return(out)
  }
  else{
    return(exp(out))
  }
}


#' Construct the basis functions for the splines evaluated on a grid
#' 
#' @param xgrid grid where we want to evaluate the spline functions, length T
#' @param knots vector of knots for the splines
#' @return (T x K+1) - matrix representing the value of each basis function 
B_basis <- function(xgrid, knots){
  n <- length(xgrid) # number of grid points where to evaluate the spline
  K <- length(knots) # number of knots
  delta <- knots[2] - knots[1]
  
  B <- array(0, dim = c(n, K + 1))
  for (j in 1:(K-1)){
    act_idx <- (1:n)[(xgrid >= knots[j]) & (xgrid <= knots[j+1])]
    act_x <- xgrid[act_idx]
    resc_x <- (act_x - knots[j]) / (knots[j+1] - knots[j])
    
    B[act_idx,j] <- (1/2) * (1 - resc_x)^2
    B[act_idx,j+1] <- -(resc_x^2) + resc_x + 1/2
    B[act_idx,j+2] <- (resc_x^2)/2
  }
  
  return(B)
}


#' Computes the Hamming Ball centered at x with radius r
#' 
#' @param x center of the Hamming Ball
#' @param S number of states 
#' @param r radius of the Hamming Ball
#' @return Hamming Ball
H_ball <- function(x, S, r){
  K <- length(x) # number of chains
  card_m <- (S - 1)^(0:r) * choose(K, 0:r) 
  # sum(card_m) is the cardinality of each HB with radius r
  
  # First, put the current vector in its HB
  HB_temp <- matrix(x, K, 1)
  
  for (m in 1:r){ # loop the possible HB radius
    HB_diff <- matrix(rep(x, card_m[m+1]), K, card_m[m+1])
    index <- utils::combn(K, m) # what elements of the vector to change?
    
    for (j in 1:ncol(index)){
      vec <- NULL
      for (i in 1:nrow(index)){ 
        vec <- c(vec, (1:S)[-x[index[i,j]]])
      }
      prop_col <- t(gtools::permutations(n = length(unique(vec)), r = m, 
                                         v = vec, repeats.allowed = TRUE))
      keep_col <- prop_col[,which(colSums(prop_col != x[index[,j]]) == m)]
      HB_diff[index[,j], ((j-1)*(card_m[m+1]/ncol(index))+1):(j*card_m[m+1]/ncol(index))] <- keep_col
    }
    HB_temp <- cbind(HB_temp,HB_diff) # save these in the Hamming Ball
  }
  
  return (HB_temp)
}


#' Samples a configuration uniformly within the HB
#' 
#' @param x center of the Hamming Ball
#' @param S number of states 
#' @param r radius of the Hamming Ball
#' @return sampled state
H_ball_unif <- function(x, S, r){
  HB_temp <- H_ball(x, S, r)
  HB_samp <- HB_temp[,sample(1:ncol(HB_temp), 1)]
  
  return (HB_samp)
}


#' Metropolis Hastings step to update the variance parameter and the smoothness
#' parameter for the random effects
#'
#' @param sigma2_ua_old variance parameter at the previous iteration
#' @param sigma2_us_old smoothness parameter at the previous iteration
#' @param beta_u_old random effects at the current iteration
#' @param hypers hyperparameters for the prior variances
#' @param n_ind number of participants (random effects)
#' @return sampled state
sample_smooth_var <- function(sigma2_ua_old, sigma2_us_old, 
                              beta_u_old, hypers, n_ind){
  
  J <- ncol(beta_u_old)
  
  # Update \sigma_{u,a}^{2}
  sigma2_ua_prop <- exp(rnorm(1, mean = log(sigma2_ua_old), sd = 0.2))
  
  lu <- log(runif(1))
  log_prop <- 0.5 * n_ind * log(det(diag(J)/sigma2_ua_prop + hypers$P/sigma2_us_old)) -
    0.5 * sum(diag(beta_u_old %*% tcrossprod(diag(J)/sigma2_ua_prop, beta_u_old))) -
    log(1 + sigma2_ua_prop^2)
  log_old <- 0.5 * n_ind * log(det(diag(J)/sigma2_ua_old + hypers$P/sigma2_us_old)) -
    0.5 * sum(diag(beta_u_old %*% tcrossprod(diag(J)/sigma2_ua_old, beta_u_old))) -
    log(1 + sigma2_ua_old^2)
  alpha <- min(c(0, log_prop + log(sigma2_ua_prop) -
                   log_old - log(sigma2_ua_old)))
  if (lu < alpha){
    sigma2_ua_old <- sigma2_ua_prop
  }
  
  # Update \sigma_{u,s}^{2}
  sigma2_us_prop <- exp(rnorm(1, mean = log(sigma2_us_old), sd = 0.1))
  while (sigma2_us_prop > 0.3){
    sigma2_us_prop <- exp(rnorm(1, mean = log(sigma2_us_old), sd = 0.1))
  }
  
  lu <- log(runif(1))
  log_prop <- 0.5 * n_ind * log(det(diag(J)/sigma2_ua_old + hypers$P/sigma2_us_prop)) -
    0.5 * sum(diag(beta_u_old %*% tcrossprod(hypers$P/sigma2_us_prop, beta_u_old))) -
    # (1 + a_sigma) * log(sigma2_us_prop) - b_sigma/sigma2_us_prop
    log(1 + sigma2_us_prop^2)
  log_old <- 0.5 * n_ind * log(det(diag(J)/sigma2_ua_old + hypers$P/sigma2_us_old)) -
    0.5 * sum(diag(beta_u_old %*% tcrossprod(hypers$P/sigma2_us_old, beta_u_old))) -
    # (1 + a_sigma) * log(sigma2_us_old) - b_sigma/sigma2_us_old
    log(1 + sigma2_us_old^2)
  alpha <- min(c(0, log_prop + log(sigma2_us_prop) -
                   log_old - log(sigma2_us_old)))
  if (lu < alpha){
    sigma2_us_old <- sigma2_us_prop
  }
  
  return(list('sigma2_us_old' = sigma2_us_old, 
              'sigma2_ua_old' = sigma2_ua_old))
  
}



#' Independent DDM fitted for every block separately
#'
#' @param tau vector of response times
#' @param ind vector of indicators denoting the participants
#' @param time vector of training blocks
#' @param trial vector of indicators denoting the trial number within each block
#' @param cens vector of censoring indicators
#' @param D matrix of covariates X = \{s, d\}
#' @param Niter number of iterations
#' @param burnin how many samples to discard
#' @param thin thinning of the chain
#' @return list with MCMC samples
indep_DDM_deltai <- function(tau, ind, time, trial, cens, D, 
                             Niter = 5000, burnin = 2000, 
                             thin = 5){
  
  samp_size <- (Niter - burnin)/thin # sample size
  p <- ncol(D) # number of covariates
  d_j <- rep(0, p) # Number of levels for each covariate
  for (j in 1:p){
    d_j[j] <- length(unique(D[!is.na(D[,j]),j]))
  }
  n <- length(tau) # total number of observations
  n_ind <- length(unique(ind)) # number of individuals
  L <- max(trial) # max number of trials
  
  # Set MCMC objects
  post_mean_delta <- array(NA, dim = c(samp_size, d_j[1]))
  post_mean_mu <- array(NA, dim = c(samp_size, T_max, d_j))
  post_mean_b <- array(NA, dim = c(samp_size, T_max, d_j))
  post_ind_delta <- array(NA, dim = c(samp_size, n_ind, d_j[1]))
  post_ind_mu <- array(NA, dim = c(samp_size, T_max, n_ind, d_j))
  post_ind_b <- array(NA, dim = c(samp_size, T_max, n_ind, d_j))
  pred_time <- array(NA, dim = c(T_max, d_j[1], samp_size))
  pred_ans <- array(NA, dim = c(T_max, d_j[1], samp_size))
  pred_time_ind <- array(NA, dim = c(n_ind, T_max, d_j[1], samp_size))
  pred_ans_ind <- array(NA, dim = c(n_ind, T_max, d_j[1], samp_size))
  sigma_u_b_chain <- array(NA, dim = samp_size)
  sigma_u_mu_chain <- array(NA, dim = samp_size)
  u_S_mu_chain <- array(NA, dim = c(samp_size, T_max, n_ind))
  u_F_mu_chain <- array(NA, dim = c(samp_size, T_max, n_ind))
  u_S_b_chain <- array(NA, dim = c(samp_size, T_max, n_ind))
  u_F_b_chain <- array(NA, dim = c(samp_size, T_max, n_ind))
  
  # Set initial values
  delta_old <- array(NA, dim = c(d_j[1], n_ind))
  mu_old <- array(0, dim = c(T_max, d_j))
  b_old <- array(0, dim = c(T_max, d_j))
  u_S_mu_old <- array(0, dim = c(T_max, n_ind))
  u_S_b_old <- array(0, dim = c(T_max, n_ind))
  u_F_mu_old <- array(0, dim = c(T_max, n_ind))
  u_F_b_old <- array(0, dim = c(T_max, n_ind))
  delta_dat <- array(NA, dim = n)
  for (s_temp in 1:d_j[1]){
    for (i_temp in 1:n_ind){
      idx_temp <- which((D[,1] == s_temp) & (ind == i_temp))
      delta_old[s_temp,i_temp] <- min(tau[which((D[,1] == s_temp) & (ind == i_temp))])/2
      delta_dat[idx_temp] <- delta_old[s_temp,i_temp]
    }
  }
  
  # MH proposal parameters
  sd_MH_delta <- array(0.1, dim = c(d_j[1], n_ind))
  sd_MH_mu <- array(0.1, dim = c(T_max, d_j))
  sd_MH_b <- array(0.1, dim = c(T_max, d_j))
  sd_MH_u_mu <-  array(0.6, dim = c(T_max, n_ind))
  sd_MH_u_b <- array(0.6, dim = c(T_max, n_ind))
  acc_delta <- array(0, dim = c(d_j[1], n_ind))
  acc_mu <- array(0, dim = c(T_max, d_j))
  acc_b <- array(0, dim = c(T_max, d_j))
  acc_u_mu <- array(0, dim = c(T_max, n_ind))
  acc_u_b <- array(0, dim = c(T_max, n_ind))
  n_batch <- 0
  
  # Auxiliary variables
  f_mu_dat <- array(0, dim = c(n, d_j[1]))
  u_mu_dat <- array(0, dim = c(n, d_j[1]))
  f_b_dat <- array(0, dim = c(n, d_j[1]))
  u_b_dat <- array(0, dim = c(n, d_j[1]))
  b_dat <- exp(f_b_dat + u_b_dat)
  mu_dat <- exp(f_mu_dat + u_mu_dat)
  sigma_u_b <- 0.5
  sigma_u_mu <- 0.5
  
  
  # Gibbs Sampler
  it <- 1
  pb <- txtProgressBar(style = 3)
  for (iter in 1:Niter){
    
    # (1) Update of the delta parameter: \delta_{x}: MH with log normal proposal
    for (s_temp in 1:d_j[1]){
      for (i_temp in 1:n_ind){
        idx_temp <- which((D[,1] == s_temp) & (ind == i_temp))
        tau_temp <- tau[idx_temp]
        cens_temp <- cens[idx_temp]
        D_temp <- D[idx_temp,]
        
        # log-normal proposal distribution centered on the current value
        delta_prop <- exp(rnorm(1, log(delta_old[s_temp,i_temp]), sd_MH_delta[s_temp,i_temp]))
        
        loglik_prop <- log_likelihood_cpp(tau_temp, mu_dat[idx_temp,], 
                                          b_dat[idx_temp,], 
                                          rep(delta_prop, length(idx_temp)), 
                                          cens_temp, D_temp, TRUE)
        loglik_old <- log_likelihood_cpp(tau_temp, mu_dat[idx_temp,], 
                                         b_dat[idx_temp,], 
                                         rep(delta_old[s_temp,i_temp], length(idx_temp)), 
                                         cens_temp, D_temp, TRUE)
        
        alpha_acc <- min(0, loglik_prop + log(delta_prop) -
                           loglik_old - log(delta_old[s_temp,i_temp]))
        l_u <- log(runif(1))
        
        if (l_u < alpha_acc){
          delta_old[s_temp,i_temp] <- delta_prop
          delta_dat[idx_temp] <- delta_old[s_temp,i_temp]
          acc_delta[s_temp,i_temp] <- acc_delta[s_temp,i_temp] + 1
        }
        
      }
    }
    
    
    for (k in T_min:T_max){ 
      for (s_temp in 1:d_j[1]){
        idx_i <- which( (D[,1] == s_temp) & (time == k) )
        tau_temp <- tau[idx_i]
        cens_temp <- cens[idx_i]
        D_temp <- D[idx_i,]
        for (d_temp in 1:d_j[2]){
          # (2) Update of b parameter: b{x,t}: MH with normal proposal
          b_prop <- rnorm(1, b_old[k,s_temp,d_temp], sd_MH_b[k,s_temp,d_temp])
          b_prop_dat <- b_dat[idx_i,]
          b_prop_dat[,d_temp] <- exp(b_prop + u_b_dat[idx_i,d_temp])
          
          logpost_prop <- log_likelihood_cpp(tau_temp, mu_dat[idx_i,], 
                                             b_prop_dat, delta_dat[idx_i], 
                                             cens_temp, D_temp, TRUE)
          logpost_old <- log_likelihood_cpp(tau_temp, mu_dat[idx_i,], 
                                            b_dat[idx_i,], delta_dat[idx_i], 
                                            cens_temp, D_temp, TRUE)
          
          alpha_acc <- min(0, logpost_prop - logpost_old)
          l_u <- log(runif(1))
          
          if (l_u < alpha_acc){
            b_old[k,s_temp,d_temp] <- b_prop
            b_dat[idx_i,] <- b_prop_dat
            f_b_dat[idx_i,d_temp] <- b_prop
            acc_b[k,s_temp,d_temp] <- acc_b[k,s_temp,d_temp] + 1
          }
          
          
          # (3) Update of \mu parameter: \mu{x,t}: MH with normal proposal
          mu_prop <- rnorm(1, mu_old[k,s_temp,d_temp], sd_MH_mu[k,s_temp,d_temp])
          mu_prop_dat <- mu_dat[idx_i,]
          mu_prop_dat[,d_temp] <- exp(mu_prop + u_mu_dat[idx_i,d_temp])
          
          logpost_prop <- log_likelihood_cpp(tau_temp, mu_prop_dat, 
                                             b_dat[idx_i,], delta_dat[idx_i], 
                                             cens_temp, D_temp, TRUE)
          logpost_old <- log_likelihood_cpp(tau_temp, mu_dat[idx_i,], 
                                            b_dat[idx_i,], delta_dat[idx_i], 
                                            cens_temp, D_temp, TRUE)
          
          alpha_acc <- min(0, logpost_prop - logpost_old)
          l_u <- log(runif(1))
          
          if (l_u < alpha_acc){
            mu_old[k,s_temp,d_temp] <- mu_prop
            mu_dat[idx_i,] <- mu_prop_dat
            f_mu_dat[idx_i,d_temp] <- mu_prop
            acc_mu[k,s_temp,d_temp] <- acc_mu[k,s_temp,d_temp] + 1
          }
          
        }
      }
    }
    
    
    for (k in T_min:T_max){ # loop over locations
      for (ind_temp in 1:n_ind){
        idx_i <- which( (ind == ind_temp) & (time == k) )
        tau_temp <- tau[idx_i]
        cens_temp <- cens[idx_i]
        D_temp <- D[idx_i,]
        
        if (length(idx_i) > 1){
          
          # (4) Update of u_b_F, u_b_S parameters:
          u_F_b_prop <- rnorm(1, u_F_b_old[k,ind_temp], sd_MH_u_b[k,ind_temp])
          u_S_b_prop <- rnorm(1, u_S_b_old[k,ind_temp], sd_MH_u_b[k,ind_temp])
          
          b_prop_dat <- f_b_dat[idx_i,] + u_F_b_prop
          b_prop_dat[cbind(1:nrow(D_temp), D_temp[,1])] <- b_prop_dat[cbind(1:nrow(D_temp), D_temp[,1])] + u_S_b_prop - u_F_b_prop
          b_prop_dat <- exp(b_prop_dat)
          
          logpost_prop <- log_likelihood_cpp(tau_temp, mu_dat[idx_i,],
                                             b_prop_dat, delta_dat[idx_i],
                                             cens_temp, D_temp, TRUE) +
            dnorm(u_F_b_prop, 0, sigma_u_b, log = TRUE) + 
            dnorm(u_S_b_prop, 0, sigma_u_b, log = TRUE)
          logpost_old <- log_likelihood_cpp(tau_temp, mu_dat[idx_i,],
                                            b_dat[idx_i,], delta_dat[idx_i],
                                            cens_temp, D_temp, TRUE) +
            dnorm(u_F_b_old[k,ind_temp], 0, sigma_u_b, log = TRUE) + 
            dnorm(u_S_b_old[k,ind_temp], 0, sigma_u_b, log = TRUE)
          
          alpha_acc <- min(0, logpost_prop - logpost_old)
          l_u <- log(runif(1))
          
          if (l_u < alpha_acc){
            u_F_b_old[k,ind_temp] <- u_F_b_prop
            u_S_b_old[k,ind_temp] <- u_S_b_prop
            b_dat[idx_i,] <- b_prop_dat
            u_b_dat[idx_i,] <- u_F_b_prop
            u_b_dat[cbind(idx_i, D_temp[,1])] <- u_b_dat[cbind(idx_i, D_temp[,1])] - u_F_b_prop + u_S_b_prop
          }
          
          
          # (5) Update of u_mu_F, u_mu_S parameters:
          u_F_mu_prop <- rnorm(1, u_F_mu_old[k,ind_temp], sd_MH_u_mu[k,ind_temp])
          u_S_mu_prop <- rnorm(1, u_S_mu_old[k,ind_temp], sd_MH_u_mu[k,ind_temp])
          
          mu_prop_dat <- f_mu_dat[idx_i,] + u_F_mu_prop
          mu_prop_dat[cbind(1:nrow(D_temp), D_temp[,1])] <- mu_prop_dat[cbind(1:nrow(D_temp), D_temp[,1])] + u_S_mu_prop - u_F_mu_prop
          mu_prop_dat <- exp(mu_prop_dat)
          
          logpost_prop <- log_likelihood_cpp(tau_temp, mu_prop_dat,
                                             b_dat[idx_i,], delta_dat[idx_i],
                                             cens_temp, D_temp, TRUE) +
            dnorm(u_F_mu_prop, 0, sigma_u_mu, log = TRUE) + 
            dnorm(u_S_mu_prop, 0, sigma_u_mu, log = TRUE)
          logpost_old <- log_likelihood_cpp(tau_temp, mu_dat[idx_i,],
                                            b_dat[idx_i,], delta_dat[idx_i],
                                            cens_temp, D_temp, TRUE) +
            dnorm(u_F_mu_old[k,ind_temp], 0, sigma_u_mu, log = TRUE) + 
            dnorm(u_S_mu_old[k,ind_temp], 0, sigma_u_mu, log = TRUE)
          
          alpha_acc <- min(0, logpost_prop - logpost_old)
          l_u <- log(runif(1))
          
          if (l_u < alpha_acc){
            u_F_mu_old[k,ind_temp] <- u_F_mu_prop
            u_S_mu_old[k,ind_temp] <- u_S_mu_prop
            mu_dat[idx_i,] <- mu_prop_dat
            u_mu_dat[idx_i,] <- u_F_mu_prop
            u_mu_dat[cbind(idx_i, D_temp[,1])] <- u_mu_dat[cbind(idx_i, D_temp[,1])] - u_F_mu_prop + u_S_mu_prop
          }
        }
      }
    }
    
    
    # (7) Update the variance of the random effects:
    RSS1 <- sum(u_S_b_old^2) + sum(u_F_b_old^2)
    RSS2 <- sum(u_S_mu_old^2) + sum(u_F_mu_old^2)
    sigma_u_b <- sqrt(1/rgamma(1, 10 + n_ind*T_max, 10 + 0.5 * RSS1))
    sigma_u_mu <- sqrt(1/rgamma(1, 10 + n_ind*T_max, 10 + 0.5 * RSS2))
    
    
    # After burnin, save parameters in the chain
    if ( (iter > burnin) && (iter %% thin == 0) ){
      for (k in 1:T_max){
        for (s_temp in 1:d_j[1]){
          mu_temp <- exp(mu_old[k,s_temp,] + 0.5 * sigma_u_mu^2)
          b_temp <- exp(b_old[k,s_temp,] + 0.5 * sigma_u_b^2)
          delta_temp <- mean(delta_old[s_temp,])
          pred_temp <- delta_temp + rinvgaussian(d_j[2], b_temp/mu_temp, 
                                                 b_temp^2)
          pred_ans[k,s_temp,it] <- which.min(pred_temp)
          pred_time[k,s_temp,it] <- min(pred_temp)
          
          for (ind_temp in 1:n_ind){
            mu_temp <- mu_old[k,s_temp,]
            mu_temp[s_temp] <- mu_temp[s_temp] + u_S_mu_old[k,ind_temp]
            mu_temp[-s_temp] <- mu_temp[-s_temp] + u_F_mu_old[k,ind_temp]
            mu_temp <- exp(mu_temp)
            
            b_temp <- b_old[k,s_temp,]
            b_temp[s_temp] <- b_temp[s_temp] + u_S_b_old[k,ind_temp]
            b_temp[-s_temp] <- b_temp[-s_temp] + u_F_b_old[k,ind_temp]
            b_temp <- exp(b_temp)
            
            delta_temp <- delta_old[s_temp,ind_temp]
            
            post_ind_mu[it,k,ind_temp,s_temp,] <- mu_temp
            post_ind_b[it,k,ind_temp,s_temp,] <- b_temp
            post_ind_delta[it,ind_temp,s_temp] <- delta_temp
            
            pred_temp <- delta_temp + rinvgaussian(d_j[2], b_temp/mu_temp,
                                                   b_temp^2)
            pred_ans_ind[ind_temp,k,s_temp,it] <- which.min(pred_temp)
            pred_time_ind[ind_temp,k,s_temp,it] <- min(pred_temp)
          }
        }
      }
      post_mean_delta[it,] <- rowMeans(delta_old)
      post_mean_mu[it,,,] <- exp(mu_old + 0.5 * sigma_u_mu^2)
      post_mean_b[it,,,] <- exp(b_old + 0.5 * sigma_u_b^2)
      sigma_u_b_chain[it] <- sigma_u_b
      sigma_u_mu_chain[it] <- sigma_u_mu
      loglik_chain[it] <- log_likelihood_cpp(tau, mu_dat, b_dat, delta_dat, cens, D, TRUE)
      u_S_mu_chain[it,,] <- u_S_mu_old
      u_F_mu_chain[it,,] <- u_F_mu_old
      u_S_b_chain[it,,] <- u_S_b_old
      u_F_b_chain[it,,] <- u_F_b_old
      
      it <- it + 1
    }
    setTxtProgressBar(pb, iter/Niter)
  }
  
  return(list('post_mean_delta' = post_mean_delta, 
              'post_mean_mu' = post_mean_mu,
              'post_mean_b' = post_mean_b,
              'post_ind_delta' = post_ind_delta,
              'post_ind_mu' = post_ind_mu,
              'post_ind_b' = post_ind_b,
              'pred_ans' = pred_ans, 
              'pred_time' = pred_time,
              'pred_ans_ind' = pred_ans_ind,
              'pred_time_ind' = pred_time_ind,
              'loglik' = loglik_chain, 
              'sigma_u_b' = sigma_u_b_chain, 
              'sigma_u_mu' = sigma_u_mu_chain, 
              'u_S_mu' = u_S_mu_chain, 
              'u_F_mu' = u_F_mu_chain, 
              'u_S_b' = u_S_b_chain, 
              'u_F_b' = u_F_b_chain
  ))
}



#' Main function for the Gibbs sampler for the drift-diffusion model
#'
#' @param tau vector of response times
#' @param ind vector of indicators denoting the participants
#' @param time vector of training blocks
#' @param trial vector of indicators denoting the trial number within each block
#' @param cens vector of censoring indicators
#' @param D matrix of covariates X = \{s, d\}
#' @param knots knots for the spline basis
#' @param xgrid grid where to evaluate the functional parameters
#' @param hypers hyperparameters of the MCMC
#' @param Niter total number of iterations
#' @param burnin burnin of the chain
#' @param thin thinning factor
#' @return list with MCMC samples
loc_clust_fHMM_deltai <- function(tau, ind, time, trial, cens, D, knots, 
                                  xgrid, hypers, Niter = 5000, burnin = 2000, 
                                  thin = 5){
  
  B <- B_basis(data$block, knots)
  Bgrid <- B_basis(xgrid, knots)
  
  samp_size <- (Niter - burnin)/thin # sample size
  p <- ncol(D) # number of covariates
  d_j <- rep(0, p) # Number of levels for each covariate
  for (j in 1:p){
    d_j[j] <- length(unique(D[!is.na(D[,j]),j]))
  }
  J <- ncol(B) # number of locations
  n <- nrow(B) # total number of observations
  n_ind <- length(unique(ind)) # number of individuals
  L <- max(trial) # max number of trials
  
  # Rescale the time steps in \{1, ..., T_max\}
  T_max <- max(time)
  T_min <- min(time)
  time <- time - T_min + 1
  T_max <- max(time)
  T_min <- min(time)
  idx_xy <- t(apply(expand.grid(y = 1:d_j[2], x = 1:d_j[1]), 1, rev))
  colnames(idx_xy) <- NULL
  
  
  # Set HB structures
  r_HB <- 1
  Z_max <- min(prod(d_j), 6)
  dim_HB <- sum((Z_max - 1)^(0:r_HB) * choose(d_j[2], 0:r_HB))
  
  
  # Set MCMC objects
  z <- array(NA, dim = c(prod(d_j), J, samp_size))
  post_mean_delta <- array(NA, dim = c(samp_size, d_j[1]))
  post_mean_mu <- array(NA, dim = c(nrow(Bgrid), d_j, samp_size))
  post_mean_b <- array(NA, dim = c(nrow(Bgrid), d_j, samp_size))
  post_ind_delta <- array(NA, dim = c(d_j[1], n_ind, samp_size))
  post_ind_mu <- array(NA, dim = c(nrow(Bgrid), n_ind, d_j, samp_size))
  post_ind_b <- array(NA, dim = c(nrow(Bgrid), n_ind, d_j, samp_size))
  sigma2_mu_us <- array(NA, dim = samp_size)
  sigma2_mu_ua <- array(NA, dim = samp_size)
  sigma2_b_us <- array(NA, dim = samp_size)
  sigma2_b_ua <- array(NA, dim = samp_size)
  sigma2_1_mu <- array(NA, dim = samp_size)
  sigma2_1_b <- array(NA, dim = samp_size)
  pred_time <- array(NA, dim = c(T_max, d_j[1], samp_size))
  pred_ans <- array(NA, dim = c(T_max, d_j[1], samp_size))
  pred_time_ind <- array(NA, dim = c(n_ind, T_max, d_j[1], samp_size))
  pred_ans_ind <- array(NA, dim = c(n_ind, T_max, d_j[1], samp_size))
  loglik_chain <- array(NA, dim = samp_size)
  
  # Set initial values
  delta_old <- array(NA, dim = c(d_j[1], n_ind))
  beta_mu_old <- array(NA, dim = c(J, d_j))
  beta_b_old <- array(NA, dim = c(J, d_j))
  delta_dat <- array(NA, dim = n)
  for (s_temp in 1:d_j[1]){
    for (i_temp in 1:n_ind){
      idx_temp <- which((D[,1] == s_temp) & (ind == i_temp))
      delta_old[s_temp,i_temp] <- min(tau[which((D[,1] == s_temp) & (ind == i_temp))])/2
      delta_dat[idx_temp] <- delta_old[s_temp,i_temp]
    }
    for(j in 1:d_j[1]){
      beta_mu_old[,s_temp,j] <- rep(0.5*log(mean(tau[D[,1] == s_temp])) - log(sd(tau[D[,1] == s_temp])), J)
      beta_b_old[,s_temp,j] <- rep(1.5*log(mean(tau[D[,1] == s_temp])) - log(sd(tau[D[,1] == s_temp])), J)
    }
  }
  low_bound_mu <- min(beta_mu_old) - 1.5
  upp_bound_mu <- max(beta_mu_old) + 1
  low_bound_b <- min(beta_b_old) - 1.5
  upp_bound_b <- max(beta_b_old) + 1
  
  beta_mu_star_prop <- array(NA, dim = c(J, Z_max))
  beta_mu_u_old <- array(0, dim = c(n_ind, J))
  beta_b_star_prop <- array(NA, dim = c(J, Z_max))
  beta_b_u_old <- array(0, dim = c(n_ind, J))
  sigma2_1_mu_old <- 0.05
  sigma2_1_b_old <- 0.05
  sigma2_b_us_old <- 0.05
  sigma2_mu_us_old <- 0.05
  sigma2_b_ua_old <- 0.05
  sigma2_mu_ua_old <- 0.05
  
  
  # Message passing structures
  beta_mess <- array(NA, dim = c(J, dim_HB))
  z_old <- list()
  for (j in 1:d_j[1]){
    z_old[[j]] <- array(NA, dim = c(d_j[2], J))
    for (jj in 1:d_j[2]){
      if (j == jj){
        z_old[[j]][jj,] <- j
      }
      else{
        z_old[[j]][jj,] <- sample((d_j[1] + 1):Z_max, 1)
      }
    }
  }
  z_temp <- do.call(rbind, z_old)
  
  
  beta_mu_star_old <- array(NA, dim = c(J, Z_max))
  beta_b_star_old <- array(NA, dim = c(J, Z_max))
  for (i in 1:Z_max){
    beta_mu_star_old[,i] <- beta_mu_old[1,idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],1],
                                        idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],2]] 
    beta_b_star_old[,i] <- beta_b_old[1,idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],1],
                                      idx_xy[which(apply(z_temp == i, 1, prod) == 1)[1],2]] 
  }
  
  rand_mat <- array(rnorm(prod(d_j)), dim = d_j)
  idx_succ <- which(rand_mat == diag(rand_mat))
  v_old <- array(NA, dim = c(d_j[2], J))
  z_prop <- list()
  
  # Transition dynamics objects
  alpha_S_old <- alpha_F_old <- 1
  Q_S_old <- array(NA, dim = c(Z_max, Z_max))
  Q_F_old <- array(NA, dim = c(Z_max, Z_max))
  pi_S_0 <- rdirichlet(rep(alpha_S_old / Z_max, Z_max) + 
                         table_int(z_temp[idx_succ,1], Z_max))
  pi_F_0 <- rdirichlet(rep(alpha_F_old / Z_max, Z_max) + 
                         table_int(z_temp[-idx_succ,1], Z_max))
  tr_count_S <- count_assign(z_temp[idx_succ,], Z_max)
  tr_count_F <- count_assign(z_temp[-idx_succ,], Z_max)
  for (h in 1:Z_max){
    Q_S_old[h,] <- rdirichlet(rep(alpha_S_old / Z_max, Z_max) + tr_count_S[h,])
    Q_F_old[h,] <- rdirichlet(rep(alpha_F_old / Z_max, Z_max) + tr_count_F[h,])
  }
  
  
  # Auxiliary variables
  prob_mat <- array(NA, dim = c(dim_HB, dim_HB))
  prob_vec <- array(NA, dim = dim_HB)
  
  # MH proposal parameters
  sd_MH_delta <- array(0.3, dim = c(d_j[1], n_ind))
  sd_MH_beta_mu <- array(0.4, dim = c(J, Z_max))
  sd_MH_beta_b <- array(0.25, dim = c(J, Z_max))
  sd_beta_mu_u <-  array(0.4, dim = c(n_ind, J))
  sd_beta_b_u <-  array(0.4, dim = c(n_ind, J))
  acc_delta <- array(0, dim = c(d_j[1], n_ind))
  acc_beta_mu <- array(0, dim = c(J, Z_max))
  acc_beta_b <- array(0, dim = c(J, Z_max))
  acc_beta_mu_u <- array(0, dim = c(n_ind, J))
  acc_beta_b_u <- array(0, dim = c(n_ind, J))
  n_batch <- 0
  
  
  # Auxiliary variables
  B_beta_mu_dat <- array(0, dim = c(n, d_j[1]))
  B_beta_mu_u_dat <- array(0, dim = c(n, d_j[1]))
  B_beta_b_dat <- array(0, dim = c(n, d_j[1]))
  B_beta_b_u_dat <- array(0, dim = c(n, d_j[1]))
  mu_dat <- exp(B_beta_mu_dat + B_beta_mu_u_dat)
  b_dat <- exp(B_beta_b_dat + B_beta_b_u_dat)
  
  # tau0 <- 1000
  # m0 <- floor(0.75 * burnin)
  # T_ann <- pmax(tau0^(1 - 1:Niter/m0), 1)
  # T_ann <- rep(1, Niter) # for now, do not consider simulated annealing
  
  
  # Gibbs Sampler
  it <- 1
  pb <- txtProgressBar(style = 3)
  for (iter in 1:Niter){
    
    # (0) Adaptively tune the MH variance for the proposals of \delta_{s, i}, 
    # beta_u_mu, beta_u_b
    if (iter %% 20 == 0){
      n_batch <- n_batch + 1
      delta_n <- min(0.01, n_batch^(-0.5))
      
      for (i in 1:n_ind){
        for (x_temp in 1:d_j[1]){
          if (acc_delta[x_temp,i]/iter > 0.44){
            sd_MH_delta[x_temp,i] <- exp(log(sd_MH_delta[x_temp,i]) + delta_n)
          }
          else{
            sd_MH_delta[x_temp,i] <- exp(log(sd_MH_delta[x_temp,i]) - delta_n)
          }
        }
        for (k in 1:J){
          if (acc_beta_mu_u[i,k]/iter > 0.44){
            sd_beta_mu_u[i,k] <- exp(log(sd_beta_mu_u[i,k]) + delta_n)
          }
          else{
            sd_beta_mu_u[i,k] <- exp(log(sd_beta_mu_u[i,k]) - delta_n)
          }
          if (acc_beta_b_u[i,k]/iter > 0.44){
            sd_beta_b_u[i,k] <- exp(log(sd_beta_b_u[i,k]) + delta_n)
          }
          else{
            sd_beta_b_u[i,k] <- exp(log(sd_beta_b_u[i,k]) - delta_n)
          }
        }
      }
    }
    
    
    # (1) Update of the delta parameter: \delta_{s,i}: MH with log normal 
    #     proposal
    for (s_temp in 1:d_j[1]){
      for (i_temp in 1:n_ind){
        idx_temp <- which((D[,1] == s_temp) & (ind == i_temp))
        tau_temp <- tau[idx_temp]
        cens_temp <- cens[idx_temp]
        D_temp <- D[idx_temp,]
        
        # log-normal proposal distribution centered on the current value
        delta_prop <- exp(rnorm(1, log(delta_old[s_temp,i_temp]), sd_MH_delta[s_temp,i_temp]))
        while (delta_prop > min(tau[idx_temp])){
          delta_prop <- exp(rnorm(1, log(delta_old[s_temp,i_temp]), sd_MH_delta[s_temp,i_temp]))
        }
        loglik_prop <- log_likelihood_cpp(tau_temp, mu_dat[idx_temp,], 
                                          b_dat[idx_temp,], 
                                          rep(delta_prop, length(idx_temp)), 
                                          cens_temp, D_temp, TRUE)
        loglik_old <- log_likelihood_cpp(tau_temp, mu_dat[idx_temp,], 
                                         b_dat[idx_temp,], 
                                         rep(delta_old[s_temp,i_temp], length(idx_temp)), 
                                         cens_temp, D_temp, TRUE)
        
        alpha_acc <- min(0, loglik_prop + log(delta_prop) -
                           loglik_old - log(delta_old[s_temp,i_temp]))
        l_u <- log(runif(1))
        
        if (l_u < alpha_acc){
          delta_old[s_temp,i_temp] <- delta_prop
          delta_dat[idx_temp] <- delta_old[s_temp,i_temp]
          acc_delta[s_temp,i_temp] <- acc_delta[s_temp,i_temp] + 1
        }
        
      }
    }
    
    
    # (2) Joint update of b, mu parameters: \b_{x,y}^{(i)}(t), \mu_{x,y}^{(i)}(t)
    for (k in 1:J){ # loop over locations
      if (k == 1){ # only data at t = 1 influence the first coefficient
        idx_time <- T_min
      }
      else if (k == J){ # only data at t = T influence the last coefficient
        idx_time <- T_max
      }
      else { # data at t = {k-1,k} influence the kth coefficient
        idx_time <- (k-1):k
      }
      
      for (h in 1:Z_max){ # loop over possible latent values
        # tuples (combinations of covariates) that are clustered together via 
        # the latent h
        idx_cov <- matrix(idx_xy[which(z_temp[,k] == h),], length(which(z_temp[,k] == h)), p)
        
        X_1k <- unique(idx_cov[,1]) # all possible values of x in this cluster
        X_2k <- unique(idx_cov[,2]) # all possible values of y in this cluster
        
        if (length(X_1k) > 0){ # h \in \mathcal{Z}_{j,k}: posterior update
          # Pick data with covariate levels of x clustered in group h and 
          # at the correct locations
          idx_i <- which( (D[,1] %in% X_1k) & (time %in% idx_time) )
          tau_temp <- tau[idx_i]
          cens_temp <- cens[idx_i]
          D_temp <- D[which((D[,1] %in% X_1k) & (time %in% idx_time)),]
          
          
          # Normal proposal distribution centered on the current value
          if (k == 1){
            beta_b_star_prop[,h] <- beta_b_star_old[,h]
            beta_b_star_prop[k,h] <- rnorm(1, beta_b_star_old[k,h], sd_MH_beta_b[k,h])
            beta_mu_star_prop[,h] <- beta_mu_star_old[,h]
            beta_mu_star_prop[k,h] <- rnorm(1, beta_mu_star_old[k,h], sd_MH_beta_mu[k,h])
          }
          # Normal proposal from the prior
          else if (k == J){
            beta_pre <- beta_b_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            beta_pre1 <- beta_mu_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            
            beta_b_star_prop[,h] <- beta_b_star_old[,h]
            # beta_b_star_prop[k,h] <- rnorm(1, mean(beta_pre), sqrt(sigma2_1_b_old/length(beta_pre)))
            beta_b_star_prop[k,h] <- rnorm(1, beta_pre, sqrt(sigma2_1_b_old))
            beta_mu_star_prop[,h] <- beta_mu_star_old[,h]
            # beta_mu_star_prop[k,h] <- rnorm(1, mean(beta_pre1), sqrt(sigma2_1_mu_old/length(beta_pre1)))
            beta_mu_star_prop[k,h] <- rnorm(1, beta_pre1, sqrt(sigma2_1_mu_old))
          }
          # Normal proposal from the prior
          else {
            beta_pre <- beta_b_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            beta_pre1 <- beta_mu_star_old[k-1, unique(z_temp[which(z_temp[,k] == h),k-1])]
            
            beta_b_star_prop[,h] <- beta_b_star_old[,h]
            # beta_b_star_prop[k,h] <- rnorm(1, mean(beta_pre), sqrt(sigma2_1_b_old/length(beta_pre)))
            beta_b_star_prop[k,h] <- rnorm(1, beta_pre, sqrt(sigma2_1_b_old))
            beta_mu_star_prop[,h] <- beta_mu_star_old[,h]
            # beta_mu_star_prop[k,h] <- rnorm(1, mean(beta_pre1), sqrt(sigma2_1_mu_old/length(beta_pre1)))
            beta_mu_star_prop[k,h] <- rnorm(1, beta_pre1, sqrt(sigma2_1_mu_old))
          }
          B_beta_b_prop_dat <- B_beta_b_dat[idx_i,]
          B_beta_mu_prop_dat <- B_beta_mu_dat[idx_i,]
          
          # Modify the proposed values in the corresponding positions
          for (hh in 1:nrow(idx_cov)){
            B_beta_b_prop_dat[which(D_temp[,1] == idx_cov[hh,1]),idx_cov[hh,2]] <- 
              B[idx_i[which(D_temp[,1] == idx_cov[hh,1])],] %*% beta_b_star_prop[,h]
            B_beta_mu_prop_dat[which(D_temp[,1] == idx_cov[hh,1]),idx_cov[hh,2]] <- 
              B[idx_i[which(D_temp[,1] == idx_cov[hh,1])],] %*% beta_mu_star_prop[,h]
          }
          
          # This is the proposed value for \b_{x,y}^{(i)}(t), \mu_{x,y}^{(i)}(t)
          b_prop_dat <- exp(B_beta_b_prop_dat + B_beta_b_u_dat[idx_i,])
          mu_prop_dat <- exp(B_beta_mu_prop_dat + B_beta_mu_u_dat[idx_i,])
          
          
          if (k == 1){
            beta_post <- beta_b_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            beta_post1 <- beta_mu_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            
            logpost_prop <- log_likelihood_cpp(tau_temp, mu_prop_dat, 
                                               b_prop_dat, delta_dat[idx_i], 
                                               cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_b_old * sum((beta_b_star_prop[k,h] - beta_post)^2) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_prop[k,h] - beta_post1)^2)
            
            logpost_old <- log_likelihood_cpp(tau_temp, mu_dat[idx_i,], 
                                              b_dat[idx_i,], delta_dat[idx_i], 
                                              cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_b_old * sum((beta_b_star_old[k,h] - beta_post)^2) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_old[k,h] - beta_post1)^2)
          }
          else if (k == J){
            logpost_prop <- log_likelihood_cpp(tau_temp, mu_prop_dat, 
                                               b_prop_dat, delta_dat[idx_i], 
                                               cens_temp, D_temp, TRUE)
            logpost_old <- log_likelihood_cpp(tau_temp, mu_dat[idx_i,], 
                                              b_dat[idx_i,], delta_dat[idx_i], 
                                              cens_temp, D_temp, TRUE)
          }
          else {
            beta_post <- beta_b_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            beta_post1 <- beta_mu_star_old[k+1, unique(z_temp[which(z_temp[,k] == h),k+1])]
            
            logpost_prop <- log_likelihood_cpp(tau_temp, mu_prop_dat, 
                                               b_prop_dat, delta_dat[idx_i], 
                                               cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_b_old * sum((beta_b_star_prop[k,h] - beta_post)^2) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_prop[k,h] - beta_post1)^2)
            
            logpost_old <- log_likelihood_cpp(tau_temp, mu_dat[idx_i,], 
                                              b_dat[idx_i,], delta_dat[idx_i], 
                                              cens_temp, D_temp, TRUE) -
              0.5/sigma2_1_b_old * sum((beta_b_star_old[k,h] - beta_post)^2) -
              0.5/sigma2_1_mu_old * sum((beta_mu_star_old[k,h] - beta_post1)^2)
          }
          
          # # Simulated annealing
          # alpha_acc <- min(0, 1/T_ann[iter] * (logpost_prop - logpost_old))
          alpha_acc <- min(0, logpost_prop - logpost_old)
          l_u <- log(runif(1))
          
          if (l_u < alpha_acc){
            beta_b_star_old[k,h] <- beta_b_star_prop[k,h]
            B_beta_b_dat[idx_i,] <- B_beta_b_prop_dat
            b_dat[idx_i,] <- b_prop_dat
            acc_beta_b[k,h] <- acc_beta_b[k,h] + 1
            
            beta_mu_star_old[k,h] <- beta_mu_star_prop[k,h]
            B_beta_mu_dat[idx_i,] <- B_beta_mu_prop_dat
            mu_dat[idx_i,] <- mu_prop_dat
            acc_beta_mu[k,h] <- acc_beta_mu[k,h] + 1
          }
        }
        else { # h \notin \mathcal{Z}_{1,k}: prior sampling
          beta_b_star_old[k,h] <- runif(1, low_bound_b, upp_bound_b)
          beta_mu_star_old[k,h] <- runif(1, low_bound_mu, upp_bound_mu)
        }
      }
    }
    
    
    # (3) Update the cluster assignments
    for (x_temp in 1:d_j[1]){ # loop over possible latent values
      beta_mess <- array(-Inf, dim = c(J, dim_HB))
      beta_mess[J,] <- 1/dim_HB
      
      v_old[,J] <- H_ball_unif(z_old[[x_temp]][,J], S = Z_max, r = r_HB)
      z_prop[[J]] <- H_ball(v_old[,J], S = Z_max, r = r_HB)
      for (k in (J - 1):1){
        idx_i <- which( (time == k) & (D[,1] == x_temp) )
        tau_temp <- tau[idx_i]
        cens_temp <- cens[idx_i]
        D_temp <- D[idx_i,]
        
        # (i) Sample the auxiliary variables
        v_temp <- H_ball(z_old[[x_temp]][,k], S = Z_max, r = r_HB)
        
        probs <- rep(-Inf, dim_HB)
        for (h in 1:dim_HB){
          B_beta_mu_prop_dat <- B_beta_mu_dat[idx_i,]
          B_beta_b_prop_dat <- B_beta_b_dat[idx_i,]
          
          B_beta_mu_prop_dat <- 0.5 * (beta_mu_star_old[k,v_temp[,h]] +
                                         beta_mu_star_old[k+1,z_old[[x_temp]][,k+1]])
          B_beta_b_prop_dat <- 0.5 * (beta_b_star_old[k,v_temp[,h]] +
                                        beta_b_star_old[k+1,z_old[[x_temp]][,k+1]])
          
          mu_dat_prop <- exp(t(B_beta_mu_prop_dat + t(B_beta_mu_u_dat[idx_i,])))
          b_dat_prop <- exp(t(B_beta_b_prop_dat + t(B_beta_b_u_dat[idx_i,])))
          
          probs[h] <- g_HB(log_likelihood_cpp(tau_temp, mu_dat_prop, b_dat_prop,
                                              delta_dat[idx_i], cens_temp, D_temp, TRUE))
        }
        probs <- as.numeric(normalise_log(probs))
        
        v_old[,k] <- v_temp[,sample(1:dim_HB, 1, prob = probs)]
        z_prop[[k]] <- H_ball(v_old[,k], S = Z_max, r = r_HB)
        
        
        # (ii) Pass messages backwards only in the restricted state space given
        #      by the slice
        z_kp1_temp <- which(beta_mess[k+1,] > 0)
        prob_mat <- array(-Inf, dim = c(dim_HB, dim_HB))
        
        for (h1 in z_kp1_temp){
          for (h2 in 1:dim_HB){
            
            B_beta_mu_prop_dat <- B_beta_mu_dat[idx_i,]
            B_beta_b_prop_dat <- B_beta_b_dat[idx_i,]
            B_beta_mu_prop_dat_1 <- 0.5 * (beta_mu_star_old[k,z_prop[[k]][,h2]] +
                                             beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])
            B_beta_b_prop_dat_1 <- 0.5 * (beta_b_star_old[k,z_prop[[k]][,h2]] +
                                            beta_b_star_old[k+1,z_prop[[k+1]][,h1]])
            B_beta_mu_prop_dat_2 <- 0.5 * (beta_mu_star_old[k,v_old[,k]] +
                                             beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])
            B_beta_b_prop_dat_2 <- 0.5 * (beta_b_star_old[k,v_old[,k]] +
                                            beta_b_star_old[k+1,z_prop[[k+1]][,h1]])
            
            mu_dat_prop_1 <- exp(t(B_beta_mu_prop_dat_1 + t(B_beta_mu_u_dat[idx_i,])))
            b_dat_prop_1 <- exp(t(B_beta_b_prop_dat_1 + t(B_beta_b_u_dat[idx_i,])))
            mu_dat_prop_2 <- exp(t(B_beta_mu_prop_dat_2 + t(B_beta_mu_u_dat[idx_i,])))
            b_dat_prop_2 <- exp(t(B_beta_b_prop_dat_2 + t(B_beta_b_u_dat[idx_i,])))
            
            prob_mat[h2,h1] <- log(beta_mess[k+1,h1]) -
              
              0.5 / sigma2_1_mu_old * sum((beta_mu_star_old[k,z_prop[[k]][,h2]] -
                                             beta_mu_star_old[k+1,z_prop[[k+1]][,h1]])^2) -
              0.5 / sigma2_1_b_old * sum((beta_b_star_old[k,z_prop[[k]][,h2]] -
                                            beta_b_star_old[k+1,z_prop[[k+1]][,h1]])^2) +
              
              log_likelihood_cpp(tau_temp, mu_dat_prop_1, b_dat_prop_1,
                                 delta_dat[idx_i], cens_temp, D_temp, TRUE) +
              g_HB(log_likelihood_cpp(tau_temp, mu_dat_prop_2, b_dat_prop_2,
                                      delta_dat[idx_i], cens_temp, D_temp, TRUE)) +
              sum(log(Q_F_old[cbind(z_prop[[k]][-x_temp,h2],z_prop[[k+1]][-x_temp,h1])])) +
              log(Q_S_old[z_prop[[k]][x_temp,h2],z_prop[[k+1]][x_temp,h1]])
          }
        }
        if ( sum(is.infinite(sum_rows_log(prob_mat))) == dim_HB){
          beta_mess[k,] <- 1/dim_HB
        }
        else{
          beta_mess[k,] <- as.numeric(sum_rows_log(prob_mat))
          beta_mess[k,] <- as.numeric(normalise_log(beta_mess[k,]))
        }
      }
      
      
      # (iii) Sample states forward (only on allowed states)
      idx_fail <- (1:d_j[2])[-x_temp]
      # Sample z_1
      prob_vec <- log(beta_mess[1,]) + log(pi_S_0[z_prop[[1]][x_temp,]]) +
        colSums(matrix(log(pi_F_0[z_prop[[1]][-x_temp,]]), d_j[2] - 1, dim_HB))
      prob_vec <- as.numeric(normalise_log(prob_vec))
      
      idx_samp <- sample(1:dim_HB, 1, FALSE, prob_vec)
      z_old[[x_temp]][,1] <- z_prop[[1]][,idx_samp]
      
      # Sample z_k
      for (k in 2:J){
        idx_km1 <- which( (time == k - 1) & (D[,1] == x_temp) )
        tau_temp <- tau[idx_km1]
        cens_temp <- cens[idx_km1]
        D_temp <- D[idx_km1,]
        
        prob_vec <- log(beta_mess[k,]) + 
          log(Q_S_old[cbind(z_old[[x_temp]][x_temp,k-1], z_prop[[k]][x_temp,])])
        for (kkk in idx_fail){
          prob_vec <- prob_vec + log(Q_F_old[cbind(z_old[[x_temp]][kkk,k-1], z_prop[[k]][kkk,])])
        }
        
        for (z_k_temp in which(is.finite(prob_vec))){
          B_beta_mu_prop_dat <- B_beta_mu_dat[idx_km1,]
          B_beta_b_prop_dat <- B_beta_b_dat[idx_km1,]
          
          B_beta_mu_prop_dat_1 <- 0.5 * (beta_mu_star_old[k-1,z_old[[x_temp]][,k-1]] +
                                           beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])
          B_beta_b_prop_dat_1 <- 0.5 * (beta_b_star_old[k-1,z_old[[x_temp]][,k-1]] +
                                          beta_b_star_old[k,z_prop[[k]][,z_k_temp]])
          B_beta_mu_prop_dat_2 <- 0.5 * (beta_mu_star_old[k-1,v_old[,k-1]] +
                                           beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])
          B_beta_b_prop_dat_2 <- 0.5 * (beta_b_star_old[k-1,v_old[,k-1]] +
                                          beta_b_star_old[k,z_prop[[k]][,z_k_temp]])
          
          mu_dat_prop_1 <- exp(t(B_beta_mu_prop_dat_1 + t(B_beta_mu_u_dat[idx_km1,])))
          b_dat_prop_1 <- exp(t(B_beta_b_prop_dat_1 + t(B_beta_b_u_dat[idx_km1,])))
          mu_dat_prop_2 <- exp(t(B_beta_mu_prop_dat_2 + t(B_beta_mu_u_dat[idx_km1,])))
          b_dat_prop_2 <- exp(t(B_beta_b_prop_dat_2 + t(B_beta_b_u_dat[idx_km1,])))
          
          prob_vec[z_k_temp] <- prob_vec[z_k_temp] -
            
            0.5 / sigma2_1_mu_old * sum((beta_mu_star_old[k-1,z_old[[x_temp]][,k-1]] -
                                           beta_mu_star_old[k,z_prop[[k]][,z_k_temp]])^2) -
            0.5 / sigma2_1_b_old * sum((beta_b_star_old[k-1,z_old[[x_temp]][,k-1]] -
                                          beta_b_star_old[k,z_prop[[k]][,z_k_temp]])^2) +
            
            log_likelihood_cpp(tau_temp, mu_dat_prop_1, b_dat_prop_1,
                               delta_dat[idx_km1], cens_temp, D_temp, TRUE) +
            g_HB(log_likelihood_cpp(tau_temp, mu_dat_prop_2, b_dat_prop_2,
                                    delta_dat[idx_km1], cens_temp, D_temp, TRUE))
        }
        prob_vec <- as.numeric(normalise_log(prob_vec))
        
        idx_samp <- sample(1:dim_HB, 1, FALSE, prob_vec)
        z_old[[x_temp]][,k] <- z_prop[[k]][,idx_samp]
      }
      
      
      # (4) Assign the cluster specific curves f_{\mu}
      for (y_temp in 1:d_j[2]){
        beta_mu_old[,x_temp,y_temp] <- beta_mu_star_old[cbind(1:J, z_old[[x_temp]][y_temp,])]
        beta_b_old[,x_temp,y_temp] <- beta_b_star_old[cbind(1:J, z_old[[x_temp]][y_temp,])]
      }
      B_beta_mu_dat[(D[,1] == x_temp),] <- B[D[,1] == x_temp,] %*% beta_mu_old[,x_temp,]
      B_beta_b_dat[(D[,1] == x_temp),] <- B[D[,1] == x_temp,] %*% beta_b_old[,x_temp,]
    }
    z_temp <- do.call(rbind, z_old)
    mu_dat <- exp(B_beta_mu_dat + B_beta_mu_u_dat)
    b_dat <- exp(B_beta_b_dat + B_beta_b_u_dat)
    
    
    
    # (5) Update the transition probabilities
    pi_S_0 <- rdirichlet(rep(alpha_S_old / Z_max, Z_max) + table_int(z_temp[idx_succ,1], Z_max))
    pi_F_0 <- rdirichlet(rep(alpha_F_old / Z_max, Z_max) + table_int(z_temp[-idx_succ,1], Z_max))
    tr_count_S <- count_assign(z_temp[idx_succ,], Z_max)
    tr_count_F <- count_assign(z_temp[-idx_succ,], Z_max)
    for (h in 1:Z_max){
      Q_S_old[h,] <- rdirichlet(rep(alpha_S_old / Z_max, Z_max) + 
                                  tr_count_S[h,])
      Q_F_old[h,] <- rdirichlet(rep(alpha_F_old / Z_max, Z_max) + 
                                  tr_count_F[h,])
    }
    
    # alpha_S_prop <- exp(rnorm(1, log(alpha_S_old), 0.1))
    # while ((alpha_S_prop < 0.01) | (alpha_S_prop > 10)){
    #   alpha_S_prop <- exp(rnorm(1, log(alpha_S_old), 0.1))
    # }
    # alpha_acc <- dgamma(alpha_S_prop, 1, 1, log = T) + 
    #   Z_max * lgamma(alpha_S_prop) - 
    #   Z_max^2 * lgamma(alpha_S_prop/Z_max) + 
    #   (alpha_S_prop/Z_max - 1) * sum(log(Q_S_old)) - (
    #     dgamma(alpha_S_old, 1, 1, log = T) + 
    #       Z_max * lgamma(alpha_S_old) - 
    #       Z_max^2 * lgamma(alpha_S_old/Z_max) + 
    #       (alpha_S_old/Z_max - 1) * sum(log(Q_S_old)))
    # 
    # l_u <- log(runif(1))
    # if (!is.na(alpha_acc)){
    #   if (l_u < alpha_acc){
    #     alpha_S_old <- alpha_S_prop
    #   }
    # }
    # 
    # alpha_F_prop <- exp(rnorm(1, log(alpha_F_old), 0.1))
    # while ((alpha_F_prop < 0.01) | (alpha_F_prop > 10)){
    #   alpha_F_prop <- exp(rnorm(1, log(alpha_F_old), 0.1))
    # }
    # alpha_acc <- dgamma(alpha_F_prop, 1, 1, log = T) + 
    #   Z_max * lgamma(alpha_F_prop) - 
    #   Z_max^2 * lgamma(alpha_F_prop/Z_max) + 
    #   (alpha_F_prop/Z_max - 1) * sum(log(Q_F_old)) - (
    #     dgamma(alpha_F_old, 1, 1, log = T) + 
    #       Z_max * lgamma(alpha_F_old) - 
    #       Z_max^2 * lgamma(alpha_F_old/Z_max) + 
    #       (alpha_F_old/Z_max - 1) * sum(log(Q_F_old)))
    # 
    # l_u <- log(runif(1))
    # if (!is.na(alpha_acc)){
    #   if (l_u < alpha_acc){
    #     alpha_F_old <- alpha_F_prop
    #   }
    # }
    
    # (6) Correction term for random effects.
    corr_mu <- colMeans(beta_mu_u_old)
    corr_b <- colMeans(beta_b_u_old)
    beta_mu_u_old <- t(t(beta_mu_u_old) - corr_mu)
    beta_b_u_old <- t(t(beta_b_u_old) - corr_b)
    
    for (k in 1:J){
      if (k == 1){
        idx_time <- which(time == T_min)
      }
      else if (k == J){
        idx_time <- which(time == T_max)
      }
      else {
        idx_time <- which(time %in% (k-1):k)
      }
      beta_mu_old[k,,] <- beta_mu_old[k,,] + corr_mu[k]
      beta_b_old[k,,] <- beta_b_old[k,,] + corr_b[k]
      beta_mu_star_old[k,] <- beta_mu_star_old[k,] + corr_mu[k]
      beta_b_star_old[k,] <- beta_b_star_old[k,] + corr_b[k]
      B_beta_mu_dat[idx_time,] <- B_beta_mu_dat[idx_time,] + 
        as.numeric(B[idx_time,] %*% corr_mu)
      B_beta_mu_u_dat[idx_time,] <- B_beta_mu_u_dat[idx_time,] - 
        as.numeric(B[idx_time,] %*% corr_mu)
      B_beta_b_dat[idx_time,] <- B_beta_b_dat[idx_time,] + 
        as.numeric(B[idx_time,] %*% corr_b)
      B_beta_b_u_dat[idx_time,] <- B_beta_b_u_dat[idx_time,] -
        as.numeric(B[idx_time,] %*% corr_b)
    }
    
    
    
    # (7) Update the cluster specific smoothness parameters
    RSS_mu <- RSS_b <- 0
    for (h1 in 1:d_j[1]){
      for (h2 in 1:d_j[2]){
        RSS_mu <- RSS_mu + as.numeric(crossprod(beta_mu_old[,h1,h2], hypers$P) %*% 
                                        beta_mu_old[,h1,h2])
        RSS_b <- RSS_b + as.numeric(crossprod(beta_b_old[,h1,h2], hypers$P) %*% 
                                      beta_b_old[,h1,h2])
      }
    }
    
    sigma2_1_mu_temp <- exp(rnorm(1, log(sigma2_1_mu_old), 0.1))
    while(sigma2_1_mu_temp > 0.3){
      sigma2_1_mu_temp <- exp(rnorm(1, log(sigma2_1_mu_old), 0.1))
    }
    l_alpha <- min(c(0, log(sigma2_1_mu_temp) - 
                       0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_mu_temp) - 
                       0.5/sigma2_1_mu_temp * RSS_mu + 
                       dhalfcauhy(sigma2_1_mu_temp, 1, T) -
                       (log(sigma2_1_mu_old) - 
                          0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_mu_old) -
                          0.5/sigma2_1_mu_old * RSS_mu +
                          dhalfcauhy(sigma2_1_mu_old, 1, T))))
    l_u <- log(runif(1))
    if (l_u < l_alpha){
      sigma2_1_mu_old <- sigma2_1_mu_temp
    }
    
    sigma2_1_b_temp <- exp(rnorm(1, log(sigma2_1_b_old), 0.1))
    while(sigma2_1_b_temp > 0.3){
      sigma2_1_b_temp <- exp(rnorm(1, log(sigma2_1_b_old), 0.1))
    }
    l_alpha <- min(c(0, log(sigma2_1_b_temp) - 
                       0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_b_temp) - 
                       0.5/sigma2_1_b_temp * RSS_b + 
                       dhalfcauhy(sigma2_1_b_temp, 1, T) -
                       (log(sigma2_1_b_old) - 
                          0.5 * (prod(d_j)*(J-1)) * log(sigma2_1_b_old) -
                          0.5/sigma2_1_b_old * RSS_b +
                          dhalfcauhy(sigma2_1_b_old, 1, T))))
    l_u <- log(runif(1))
    if (l_u < l_alpha){
      sigma2_1_b_old <- sigma2_1_b_temp
    }
    
    
    # (8a) Update random effects for b
    reff <- sample_reff_b(tau, D, cens, beta_b_u_old, delta_dat, B_beta_b_dat,
                          b_dat, mu_dat, B, hypers$P, ind, time,
                          sigma2_b_us_old, sigma2_b_ua_old, sd_beta_b_u, 
                          acc_beta_b_u)
    beta_b_u_old <- reff$beta_u_old
    B_beta_b_u_dat <- reff$B_beta_u_dat
    b_dat <- exp(B_beta_b_dat + B_beta_b_u_dat)
    acc_beta_b_u <- reff$acc_beta_u
    
    
    # (8b) Update random effects variances: MH with log normal proposal
    ls_var <- sample_smooth_var(sigma2_b_ua_old, sigma2_b_us_old,
                                beta_b_u_old, hypers, n_ind)
    sigma2_b_ua_old <- ls_var$sigma2_ua_old
    sigma2_b_us_old <- ls_var$sigma2_us_old
    
    
    # (9a) Update random effects for mu
    reff <- sample_reff_mu(tau, D, cens, beta_mu_u_old, delta_dat, b_dat,
                           B_beta_mu_dat, mu_dat, B, hypers$P, ind, time,
                           sigma2_mu_us_old, sigma2_mu_ua_old, sd_beta_mu_u,
                           acc_beta_mu_u)
    beta_mu_u_old <- reff$beta_u_old
    B_beta_mu_u_dat <- reff$B_beta_u_dat
    mu_dat <- exp(B_beta_mu_dat + B_beta_mu_u_dat)
    acc_beta_mu_u <- reff$acc_beta_u
    
    
    # (9b) Update random effects variances: MH with log normal proposal
    ls_var <- sample_smooth_var(sigma2_mu_ua_old, sigma2_mu_us_old,
                                beta_mu_u_old, hypers, n_ind)
    sigma2_mu_ua_old <- ls_var$sigma2_ua_old
    sigma2_mu_us_old <- ls_var$sigma2_us_old
    
    
    # After burnin, save parameters in the chain
    if ( (iter > burnin) && (iter %% thin == 0) ){
      # This is the correction for the random effects integration: we need to 
      # compute the variance of the random effects as detailed in the Supplementary
      # Materials
      cov_reff_mu <- solve(diag(J) / sigma2_mu_ua_old + hypers$P / sigma2_mu_us_old)
      corr_term_gr_mu <- rep(0, nrow(Bgrid))
      corr_term_mu <- rep(0, T_max)
      cov_reff_b <- solve(diag(J) / sigma2_b_ua_old + hypers$P / sigma2_b_us_old)
      corr_term_gr_b <- rep(0, nrow(Bgrid))
      corr_term_b <- rep(0, T_max)
      for (k in 1:J){
        corr_term_gr_mu <- corr_term_gr_mu + rowSums(Bgrid[,k] * t(t(Bgrid) * cov_reff_mu[,k]))
        corr_term_mu <- corr_term_mu + rowSums(B_basis(1:T_max, knots)[,k] * 
                                                 t(t(B_basis(1:T_max, knots)) * cov_reff_mu[,k]))
        corr_term_gr_b <- corr_term_gr_b + rowSums(Bgrid[,k] * t(t(Bgrid) * cov_reff_b[,k]))
        corr_term_b <- corr_term_b + rowSums(B_basis(1:T_max, knots)[,k] * 
                                               t(t(B_basis(1:T_max, knots)) * cov_reff_b[,k]))
      }
      
      
      for (h1 in 1:d_j[1]){
        for (h2 in 1:d_j[2]){
          post_mean_mu[,h1,h2,it] <- exp(Bgrid %*% beta_mu_old[,h1,h2] + 0.5 * corr_term_gr_mu)
          post_mean_b[,h1,h2,it] <- exp(Bgrid %*% beta_b_old[,h1,h2] + 0.5 * corr_term_gr_b)
          
          for (i in 1:n_ind){
            post_ind_mu[,i,h1,h2,it] <- exp(Bgrid %*% beta_mu_old[,h1,h2] + Bgrid %*% beta_mu_u_old[i,])
            post_ind_b[,i,h1,h2,it] <- exp(Bgrid %*% beta_b_old[,h1,h2] + Bgrid %*% beta_b_u_old[i,])
          }
        }
      }
      post_ind_delta[,,it] <- delta_old
      
      
      # Sample from the predictive distributions of response category and 
      # response times (so we can use predictive checks to evaluate goodness of fit)
      for (t in 1:T_max){
        for (x_temp in 1:d_j[1]){
          mu_temp <- as.numeric(exp(B_basis(t, knots) %*% beta_mu_old[,x_temp,] + 
                                      0.5 * corr_term_mu[t]))
          b_temp <- as.numeric(exp(B_basis(t, knots) %*% beta_b_old[,x_temp,] + 
                                     0.5 * corr_term_b[t]))
          delta_temp <- mean(delta_old[x_temp,])
          pred_temp <- delta_temp + rinvgaussian(d_j[2], b_temp/mu_temp, 
                                                 b_temp^2)
          pred_ans[t,x_temp,it] <- which.min(pred_temp)
          pred_time[t,x_temp,it] <- min(pred_temp)
          
          for (i in 1:n_ind){
            mu_temp <- as.numeric(exp(B_basis(t, knots) %*% (beta_mu_old[,x_temp,] + 
                                                               beta_mu_u_old[i,])))
            b_temp <- as.numeric(exp(B_basis(t, knots) %*% (beta_b_old[,x_temp,] + 
                                                              beta_b_u_old[i,])))
            delta_temp <- delta_old[x_temp,i]
            pred_temp <- delta_temp + rinvgaussian(d_j[2], b_temp/mu_temp, 
                                                   b_temp^2)
            pred_ans_ind[i,t,x_temp,it] <- which.min(pred_temp)
            pred_time_ind[i,t,x_temp,it] <- min(pred_temp)
          }
        }
      }
      
      # Save MCMC objects
      post_mean_delta[it,] <- rowMeans(delta_old)
      sigma2_mu_us[it] <- sigma2_mu_us_old
      sigma2_mu_ua[it] <- sigma2_mu_ua_old
      sigma2_b_us[it] <- sigma2_b_us_old
      sigma2_b_ua[it] <- sigma2_b_ua_old
      sigma2_1_mu[it] <- sigma2_1_mu_old
      sigma2_1_b[it] <- sigma2_1_b_old
      z[,,it] <- z_temp		
      
      it <- it + 1
    }
    setTxtProgressBar(pb, iter/Niter)
  }
  
  return(list('Z' = z, 
              'post_mean_delta' = post_mean_delta, 
              'post_mean_mu' = post_mean_mu,
              'post_mean_b' = post_mean_b,
              'post_ind_delta' = post_ind_delta,
              'post_ind_mu' = post_ind_mu,
              'post_ind_b' = post_ind_b,
              'sigma2_mu_us' = sigma2_mu_us, 
              'sigma2_mu_ua' = sigma2_mu_ua,
              'sigma2_b_us' = sigma2_b_us, 
              'sigma2_b_ua' = sigma2_b_ua,
              'sigma2_1_mu' = sigma2_1_mu, 
              'sigma2_1_b' = sigma2_1_b, 
              'pred_ans' = pred_ans, 
              'pred_time' = pred_time,
              'pred_ans_ind' = pred_ans_ind, 
              'pred_time_ind' = pred_time_ind
  ))
}
