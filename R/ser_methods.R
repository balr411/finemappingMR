#' @title Single effect regression that does not depend on susieR
#'
#' @description Fits the single effect regression (SER) for one effect given
#' only the score \eqn{X'r} and the diagonal \eqn{d = diag(X'X)} of the
#' (possibly stacked) regression with residual variance \code{residual_variance}.
#' This is all that is needed whether the data is individual-level (then
#' \eqn{X'r} and \eqn{d} are computed from the genotypes) or sufficient
#' statistics. It follows the implementation used for the Z-score method:
#' the prior variance is updated with \code{optim} (Brent) on the log scale,
#' the previous value is kept if it is better, and the prior variance is set to
#' 0 if the null model is better, in which case the posterior is set to 0.
#'
#' @param Xty Vector \eqn{X'r}.
#'
#' @param d Vector of the diagonal of \eqn{X'X}.
#'
#' @param V Current prior variance of the effect.
#'
#' @param residual_variance Residual variance.
#'
#' @param lower Lower bound of the log prior variance in the optimization.
#'
#' @param upper Upper bound of the log prior variance in the optimization.
#'
#' @return A list with the PIPs (\code{alpha}), posterior means (\code{mu}) and
#' second moments (\code{mu2}), the updated prior variance (\code{V}), and the
#' KL divergence of the posterior from the prior (\code{kl}, positive).
#'
#' @keywords internal

fmr_ser <- function(Xty, d, V, residual_variance = 1, lower = -30, upper = 15){
  p <- length(d)
  Z_star <- Xty/residual_variance

  #Negative log of the (unnormalized) Bayes factor sum, as a function of log(V)
  neg_loglik <- function(theta){
    omega_j <- d/residual_variance + exp(-theta)
    x <- -0.5*(theta + log(omega_j)) + Z_star^2/(2*omega_j)
    x_max <- max(x)
    -(x_max + log(sum(exp(x - x_max))))
  }

  #Update the prior variance
  res <- optim(par = ifelse(V > 0, log(V), lower),
               fn = neg_loglik,
               method = "Brent",
               lower = lower,
               upper = upper)

  if(!is.null(res$par) && res$convergence == 0){
    #First check if the new parameter beats the old
    if(V > 0 && neg_loglik(res$par) > neg_loglik(log(V))){
      V_new <- V
    }else{
      V_new <- exp(res$par)
    }

    #Now check if the parameter beats 0 (the null model)
    if(log(p) > -neg_loglik(log(V_new))){
      V_new <- 0
    }
  }else{
    warning("Prior variance update did not converge; keeping previous value")
    V_new <- V
  }

  if(V_new > 0){
    omega_j <- d/residual_variance + 1/V_new

    #Update the PIPs and posterior means
    post_update <- update_posterior_rss(Z_star, omega_j)

    #get_ser_neg_KL_divergence_rss gives the negative KL divergence
    kl <- -get_ser_neg_KL_divergence_rss(pip = post_update$alpha,
                                         prior_var = V_new,
                                         post_mean = post_update$mu,
                                         post_var = 1/omega_j)

    return(list(alpha = post_update$alpha, mu = post_update$mu,
                mu2 = post_update$mu2, V = V_new, kl = kl))
  }else{
    return(list(alpha = rep(0, p), mu = rep(0, p), mu2 = rep(0, p), V = 0, kl = 0))
  }
}
