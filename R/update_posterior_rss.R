#' @title Update the posteriors for b and alpha using Z-score derivation
#'
#' @keywords internal

update_posterior_rss <- function(Z_star_l, omega_j){
  #First update the PIPs
  #Note that the terms in the numerator can become very large, causing overflow issues
  #with the exponential, so center at the maximum
  #exp_num <- Z_star_l^2/(2*omega_j)
  #exp_num <- exp_num - max(exp_num)

  #num_pip <- omega_j^(-1/2) * exp(exp_num)
  #pip_out <- num_pip/sum(num_pip)

  #num_pip <- omega_j^(-1/2) * exp(Z_star_l^2/(2*omega_j))
  #pip_out <- num_pip/sum(num_pip)

  ####################################################
  log_num_pip <- -0.5 * log(omega_j) + Z_star_l^2 / (2 * omega_j)
  log_num_pip <- log_num_pip - max(log_num_pip)
  pip_out <- exp(log_num_pip) / sum(exp(log_num_pip))

  ####################################################

  #Now update the normal distribution parameters
  post_mean <- (1/omega_j)*Z_star_l
  post_var <- 1/omega_j
  post_mean2 <- post_var + post_mean^2

  return(list(alpha = pip_out, mu = post_mean, mu2 = post_mean2))
}
