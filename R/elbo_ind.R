#' @title Get the Bayesian ELBO when individual-level data is available
#'
#' @description Individual-level analogue of \code{bayes_elbo_ss}. Everything is
#' passed as lists with one element per locus, where \code{Xa} and \code{Xb_y}
#' are the lists of fitted values \eqn{G_{y,m}\bar{a}_m} and
#' \eqn{G_{y,m}\bar{b}_m} (the latter without \eqn{\gamma}), and \code{G_x},
#' \code{G_y} are centered and carry the \code{"d"} attribute.
#'
#' Unlike \code{bayes_elbo_ss}, the returned ELBO includes the \eqn{y'y} and
#' \eqn{x'x} terms. Hence it is equal to the sufficient statistic ELBO minus the
#' constant \eqn{y'y/(2\sigma_y^2) + x'x/(2\sigma_x^2)}.
#'
#' @keywords internal

bayes_elbo_ind <- function(sigma2_y, sigma2_x, x, y, G_x, G_y, Xa, Xb_y,
                           alpha_b, mu_b, mu2_b, alpha_a, mu_a, mu2_a,
                           mu_gamma, mu2_gamma, kl_a, kl_b, kl_gamma){
  n_x <- length(x)
  n_y <- length(y)
  M <- length(G_x)

  xtx <- sum(x^2)
  yty <- sum(y^2)

  part1 <- -((n_x + n_y)/2) * log(2*pi) - (n_y/2)*log(sigma2_y) - (n_x/2)*log(sigma2_x)

  #Loci are treated as separate (as in the sufficient statistic version), so
  #y'y and x'x enter once and the locus-specific terms are summed
  part2 <- yty
  part3 <- xtx
  for(m in 1:M){
    part2 <- part2 - 2*sum(y*Xa[[m]]) - 2*mu_gamma*sum(y*Xb_y[[m]]) +
      2*mu_gamma*sum(Xa[[m]]*Xb_y[[m]]) +
      mu2_gamma*test_ER(G_y[[m]], y, alpha_b[[m]], mu_b[[m]], mu2_b[[m]]) +
      test_ER(G_y[[m]], y, alpha_a[[m]], mu_a[[m]], mu2_a[[m]])

    part3 <- part3 + get_ER2(G_x[[m]], x, alpha_b[[m]], mu_b[[m]], mu2_b[[m]]) - xtx
  }

  part4 <- kl_gamma - sum(unlist(kl_a)) - sum(unlist(kl_b))

  elbo_curr <- part1 - (1/(2*sigma2_y))*part2 - (1/(2*sigma2_x))*part3 + part4

  return(list(elbo_curr, part2, part3))
}
