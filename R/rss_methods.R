#' @title Relevant functions for running finemappingMR with Z-score derivation
#'
#' @keywords internal

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
  #R is a list of correlation matrices (one per region); the result is summed over regions
  M <- length(R)
  vec_den <- vector(length = M)
  for(m in 1:M){
    B <- alpha_b[[m]] * mu_b[[m]]
    XB2 <- n * sum((B %*% R[[m]]) * B)
    betabar <- colSums(B)
    d <- n
    postb2 <- alpha_b[[m]] * mu2_b[[m]]
    vec_den[m] <- sum(betabar * (n*R[[m]] %*% betabar)) - XB2 + sum(d * t(postb2))
  }

  return(sum(vec_den))
}
