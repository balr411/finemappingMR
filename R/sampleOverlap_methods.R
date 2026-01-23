#' @title Relevant functions for running finemappingMR with sample overlap
#'
#' @keywords internal

negloglik_sampleOverlap_alpha <- function(sl2, RRinvR, Z_star_l, n_y, rho){
  omega_j <- (n_y/(1 - rho^2)) * RRinvR + 1/sl2

  return(-matrixStats::logSumExp(-0.5*log(sl2*omega_j) + Z_star_l^2/(2 * omega_j)))
}

negloglik_sampleOverlap_b <- function(sl2, RRinvR, Z_star_l, n_x, n_y, Egamma, Egamma2, rho){
  omega_j <- ((n_x + 2*rho*sqrt(n_x*n_y)*Egamma + n_y*Egamma2)/(1 - rho^2)) * RRinvR + 1/sl2

  return(-matrixStats::logSumExp(-0.5*log(sl2*omega_j) + Z_star_l^2/(2 * omega_j)))
}

negloglik_sampleOverlap_alpha_trans <- function(theta, RRinvR, Z_star_l, n_y, rho){
  sl2 <- exp(theta)

  omega_j <- (n_y/(1 - rho^2)) * RRinvR + 1/sl2

  return(-matrixStats::logSumExp(-0.5*log(sl2*omega_j) + Z_star_l^2/(2 * omega_j)))
}

negloglik_sampleOverlap_b_trans <- function(theta, RRinvR, Z_star_l, n_x, n_y, Egamma, Egamma2, rho){
  sl2 <- exp(theta)

  omega_j <- ((n_x + 2*rho*sqrt(n_x*n_y)*Egamma + n_y*Egamma2)/(1 - rho^2)) * RRinvR + 1/sl2

  return(-matrixStats::logSumExp(-0.5*log(sl2*omega_j) + Z_star_l^2/(2 * omega_j)))
}


elbo_sampleOverlap <- function(ZxRinvZx, ZyRinvZy, ZxRinvZy, Z_x, Z_y, mu_b, mu2_b, alpha_b,
                               mu_a, mu2_a, alpha_a, n_x, n_y, mu_gamma, mu2_gamma, R, Rinv,
                               kl_a, kl_b, kl_gamma, rho){

  a_post <- lapply(Map("*", mu_a, alpha_a), colSums)
  b_post <- lapply(Map("*", mu_b, alpha_b), colSums)

  E_btRb <- b_t_b(mu_b, mu2_b, alpha_b, R, n = 1)

  #These parts shared with the regular elbo_rss
  part1 <- ZxRinvZx + ZyRinvZy
  part2 <- -2 * crossprod((sqrt(n_x) * Z_x + sqrt(n_y) * mu_gamma * Z_y), b_post[[1]])
  part3 <- -2 * crossprod((sqrt(n_y) * Z_y), a_post[[1]])
  part4 <- (n_x + n_y*mu2_gamma) * E_btRb
  part5 <- b_t_b(mu_a, mu2_a, alpha_a, R, n = n_y)
  part6 <- 2 * n_y * mu_gamma * crossprod(a_post[[1]], R) %*% b_post[[1]]

  #These parts from the cross term
  part7 <- ZxRinvZy
  part8 <- -crossprod((sqrt(n_x) * Z_y + sqrt(n_y) * mu_gamma * Z_x), b_post[[1]])
  part9 <- -crossprod((sqrt(n_y) * Z_x), a_post[[1]])
  part10 <- sqrt(n_x*n_y) * mu_gamma * E_btRb
  part11 <- sqrt(n_x*n_y) * crossprod(a_post[[1]], R) %*% b_post[[1]]

  Eloglik <- (-0.5/(1-rho^2)) * (part1 + part2 + part3 + part4 + part5 + part6 + 2*rho*(part7 + part8 + part9 + part10 + part11))

  elbo <- Eloglik + sum(unlist(kl_a)) + sum(unlist(kl_b)) + kl_gamma

  return(elbo)
}

sampleOverlap_total_var_rss <- function(V_x, V_y, mu_b, mu2_b, alpha_b, mu_a, mu2_a,
                                        alpha_a, R, Z_y, num_samples,
                                        sigma2_gamma_prior, rho){

  M <- length(mu_b)

  idx_b_keep_full <- apply(V_x, 2, function(x) which(x > 0))
  idx_a_keep_full <- apply(V_y, 2, function(x) which(x > 0))

  E_var_full <- c()
  Var_E_full <- c()
  Var_E2_full <- c()

  for(j in 1:num_samples){
    if(j %% 50 == 0){
      print(j)
    }

    #For each region, iterate over each of the effects, sample one, and generate the effects
    b_full <- list()
    a_full <- list()

    for(m in 1:M){

      if(length(idx_b_keep_full) > 0){
        idx_b_keep <- idx_b_keep_full[[m]]
      }else{
        idx_b_keep <- integer(0)
      }

      if(length(idx_a_keep_full) > 0){
        idx_a_keep <- idx_a_keep_full[[m]]
      }else{
        idx_a_keep <- integer(0)
      }


      if(length(idx_b_keep) > 0){
        N_col <- ncol(mu_b[[m]])
        b_curr <- matrix(0, nrow = length(idx_b_keep), ncol = N_col)
        for(l in 1:length(idx_b_keep)){
          idx_curr <- idx_b_keep[l]

          #Sample the effect variant
          b_eff_var <- sample(1:N_col, 1, prob = alpha_b[[m]][idx_curr,])

          #Now sample the effect
          b_eff_size <- rnorm(1, mean = mu_b[[m]][idx_curr, b_eff_var], sd = sqrt(mu2_b[[m]][idx_curr,b_eff_var] - mu_b[[m]][idx_curr,b_eff_var]^2))

          vec_curr <- rep(0, N_col)
          vec_curr[b_eff_var] <- b_eff_size
          b_curr[l,] <- vec_curr
        }

        b_full[[m]] <- as.matrix(colSums(b_curr))
      }else{
        b_full[[m]] <- as.matrix(rep(0, ncol(mu_b[[m]])))
      }

      if(length(idx_a_keep) > 0){
        N_col <- ncol(mu_a[[m]])
        a_curr <- matrix(0, nrow = length(idx_a_keep), ncol = N_col)
        for(l in 1:length(idx_a_keep)){
          idx_curr <- idx_a_keep[l]

          #Sample the effect variant
          a_eff_var <- sample(1:N_col, 1, prob = alpha_a[[m]][idx_curr,])

          #Now sample the effect
          a_eff_size <- rnorm(1, mean = mu_a[[m]][idx_curr,a_eff_var], sd = sqrt(mu2_a[[m]][idx_curr,a_eff_var] - mu_a[[m]][idx_curr,a_eff_var]^2))

          vec_curr <- rep(0, N_col)
          vec_curr[a_eff_var] <- a_eff_size
          a_curr[l,] <- vec_curr
        }
        a_full[[m]] <- as.matrix(colSums(a_curr))
      }else{
        a_full[[m]] <- as.matrix(rep(0, ncol(mu_a[[m]])))
      }
    }


    den_curr <- (1 - rho^2 + sigma2_gamma_prior * n_y * sum(unlist(crossprod(b_full[[1]], R) %*% b_full[[1]])))
    E_var_curr <- (sigma2_gamma_prior * (1-rho^2))/den_curr

    #Now generate the estimates of the expectation and variance
    #Note when Z_y is made to be a list the following will work
    #den_curr <- (1 + sigma2_gamma_prior*sum(unlist(Map("%*%", Map(crossprod, b_full, Gy_t_Gy), b_full))))
    #E_var_curr <- sigma2_gamma_prior/den_curr

    #Again note when Z_y is made to be a list the following will work
    #num_curr <- sum(unlist(Map(crossprod, b_full, Z_y))) - sqrt(n_y)*sum(unlist(Map("%*%", Map(crossprod, b_full, R), a_full)))
    num_curr_part1 <- sqrt(n_y) * (sum(unlist(Z_y %*% b_full[[1]])) - sqrt(n_y)*sum(unlist(crossprod(b_full[[1]], R) %*% a_full[[1]])))
    num_curr_part2 <- rho*sqrt(n_y)*sum(unlist(Z_x %*% b_full[[1]]))
    num_curr_part3 <- -rho * sqrt(n_x*n_y) * sum(unlist(crossprod(b_full[[1]], R) %*% b_full[[1]]))

    Var_E_curr <- (sigma2_gamma_prior * (num_curr_part1 + num_curr_part2 + num_curr_part3))/den_curr
    Var_E2_curr <- ((sigma2_gamma_prior * (num_curr_part1 + num_curr_part2 + num_curr_part3))/den_curr)^2

    E_var_full <- c(E_var_full, E_var_curr)
    Var_E_full <- c(Var_E_full, Var_E_curr)
    Var_E2_full <- c(Var_E2_full, Var_E2_curr)
  }

  gamma_var_total_var <- mean(E_var_full) + mean(Var_E2_full) - mean(Var_E_full)^2

  return(gamma_var_total_var)

}

