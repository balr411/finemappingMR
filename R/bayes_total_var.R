#' @title Calculate the variance of gamma using the law of total variance
#'
#' @keywords internal

bayes_total_var <- function(V_x, V_y, mu_b, mu2_b, alpha_b, mu_a, mu2_a, alpha_a,
                            Gy_t_Gy, Gy_t_y, num_samples, sigma2_y, sigma2_gamma_prior){

  M <- length(mu_b)

  idx_b_keep_full <- apply(V_x, 2, function(x) which(x > 0))
  idx_a_keep_full <- apply(V_y, 2, function(x) which(x > 0))

  E_var_full <- c()
  Var_E_full <- c()
  Var_E2_full <- c()

  for(j in 1:1000){
    if(j %% 50 == 0){
      print(j)
    }

    #For each region, iterate over each of the effects, sample one, and generate the effects
    b_full <- list()
    a_full <- list()

    for(m in 1:M){
      idx_b_keep <- idx_b_keep_full[[m]]

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
        a_full[[m]] <- as.matrix(rep(0, N_col))
      }
    }

    #Now generate the estimates of the expectation and variance
    den_curr <- (sigma2_y + sigma2_gamma_prior*sum(unlist(Map("%*%", Map(crossprod, b_full, Gy_t_Gy), b_full))))
    E_var_curr <- sigma2_y*sigma2_gamma_prior/den_curr

    num_curr <- sum(unlist(Map(crossprod, b_full, Gy_t_y))) - sum(unlist(Map("%*%", Map(crossprod, b_full, Gy_t_Gy), a_full)))

    Var_E_curr <- (sigma2_gamma_prior*num_curr)/den_curr
    Var_E2_curr <- ((sigma2_gamma_prior*num_curr)/den_curr)^2

    E_var_full <- c(E_var_full, E_var_curr)
    Var_E_full <- c(Var_E_full, Var_E_curr)
    Var_E2_full <- c(Var_E2_full, Var_E2_curr)
  }

  gamma_var_total_var <- mean(E_var_full) + mean(Var_E2_full) - mean(Var_E_full)^2

  return(gamma_var_total_var)
}






