#' @title Relevant functions for running finemappingMR with Z-score derivation
#'
#' @keywords internal

negloglik_rss_alpha <- function(sl2, RRinvR, Z_star_l, n_y){
  omega_j <- n_y * RRinvR + 1/sl2

  return(-matrixStats::logSumExp(-0.5*log(sl2*omega_j) + Z_star_l^2/(2 * omega_j)))
}

negloglik_rss_b <- function(sl2, RRinvR, Z_star_l, n_x, n_y, Egamma2){

  omega_j <- (n_x + n_y*Egamma2) * RRinvR + 1/sl2

  return(-matrixStats::logSumExp(-0.5*log(sl2*omega_j) + Z_star_l^2/(2 * omega_j)))
}

get_ser_neg_KL_divergence_rss <- function(pip, prior_var, post_mean, post_var){
  #Note this this is the negative KL-divergence
  #Note that in the future if we wish to implement different prior probability
  #of being causal for each variant, it will appear in the sum(pip * log(1/(p*pip))) term

  #Note that sometimes pip = 0, which makes the const term below become NaN.
  #Other times the pip is so small that 1/(p*pip) overflows to infinity
  #To fix this, remove all variants with pip < 1e-200 since they don't contribute to the
  #KL-divergence anyways
  p <- length(pip)

  idx_nonZero <- which(pip > 1e-200)
  pip <- pip[idx_nonZero]
  post_mean <- post_mean[idx_nonZero]
  post_var <- post_var[idx_nonZero]

  kl_normal <- sum((pip/2) * (1 + log(post_var/prior_var) - (post_mean^2 + post_var)/prior_var))
  const <- sum(pip * log(1/(p * pip)))

  return(kl_normal + const)
}

b_t_b <- function(mu_b, mu2_b, alpha_b, R, n){
  #Note n can be whatever we want multiplied by R (I think)
  #Also works for alpha
  vec_den <- vector(length = 1)
  for(m in 1:1){
    B <- alpha_b[[m]] * mu_b[[m]]
    XB2 <- n * sum((B %*% R) * B)
    betabar <- colSums(B)
    d <- n
    postb2 <- alpha_b[[m]] * mu2_b[[m]]
    vec_den[m] <- sum(betabar * (n*R %*% betabar)) - XB2 + sum(d * t(postb2))
  }

  return(vec_den)
}

elbo_rss <- function(ZxRinvZx, ZyRinvZy, Z_x, Z_y, mu_b, mu2_b, alpha_b,
                     mu_a, mu2_a, alpha_a, n_x, n_y, mu_gamma, mu2_gamma, R, Rinv,
                     kl_a, kl_b, kl_gamma){

  a_post <- lapply(Map("*", mu_a, alpha_a), colSums)
  b_post <- lapply(Map("*", mu_b, alpha_b), colSums)

  part1 <- ZxRinvZx + ZyRinvZy
  part2 <- -2 * crossprod((sqrt(n_x) * Z_x + sqrt(n_y) * mu_gamma * Z_y), b_post[[1]])
  part3 <- -2 * crossprod((sqrt(n_y) * Z_y), a_post[[1]])
  part4 <- b_t_b(mu_b, mu2_b, alpha_b, R, n = n_x + n_y*mu2_gamma)
  part5 <- b_t_b(mu_a, mu2_a, alpha_a, R, n = n_y)
  part6 <- 2 * n_y * mu_gamma * crossprod(a_post[[1]], R) %*% b_post[[1]]

  Eloglik <- -0.5 * (part1 + part2 + part3 + part4 + part5 + part6)

  elbo <- Eloglik + sum(unlist(kl_a)) + sum(unlist(kl_b)) + kl_gamma

  return(elbo)
}




