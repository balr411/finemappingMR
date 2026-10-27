#' @title Get the ELBO when there is only sufficient statistics
#'
#' @keywords internal

#Must be a better way than this for loop but leaving it for now
get_ER2_ss_list <- function (XtX, Xty, alpha, mu, mu2) {
  len <- length(XtX)
  vec <- vector(length = len)
  for(i in 1:len){
    B = alpha[[i]] * mu[[i]]
    XB2 = sum((B %*% XtX[[i]]) * B)
    betabar = colSums(B)
    d = attr(XtX[[i]],"d")
    postb2 = alpha[[i]] * mu2[[i]] # Posterior second moment.
    vec[i] <- -2*sum(betabar * Xty[[i]]) + sum(betabar * (XtX[[i]] %*% betabar)) -
      XB2 + sum(d * t(postb2))
  }
  return(vec)
}

#Must be a better way than this for loop but leaving it for now
test_ER_ss_list <- function (XtX, alpha, mu, mu2){
  len <- length(XtX)
  vec <- vector(length = len)
  for(i in 1:len){
    B <- alpha[[i]] * mu[[i]]
    XB2 <- sum((B %*% XtX[[i]]) * B)
    betabar <- colSums(B)
    d <- attr(XtX[[i]],"d")
    postb2 <- alpha[[i]] * mu2[[i]]
    vec[i] <- sum(betabar * (XtX[[i]] %*% betabar)) - XB2 + sum(d * t(postb2))
  }
  return(vec)
}

#Function that computes the Bayesian ELBO
bayes_elbo_ss <- function(sigma2_y, sigma2_x, GxtGx, GytGy, GxtX, GytY, alpha_b, mu_b, mu2_b,
                          alpha_a, mu_a, mu2_a, mu_gamma, mu2_gamma, kl_a, kl_b, kl_gamma, n_x, n_y){

  part1 <- -((n_x + n_y)/2) * log(2*pi) - (n_y/2)*log(sigma2_y) - (n_x/2)*log(sigma2_x)

  b_post <- lapply(Map("*", mu_b, alpha_b), colSums)
  a_post <- lapply(Map("*", mu_a, alpha_a), colSums)

  part2 <- -2*sum(unlist(Map("%*%", a_post, GytY))) - 2*mu_gamma*sum(unlist(Map("%*%", b_post, GytY))) +
    2*mu_gamma*sum(unlist(Map("%*%", Map("%*%", b_post, GytGy), a_post))) +
    mu2_gamma*sum(test_ER_ss_list(GytGy, alpha_b, mu_b, mu2_b)) +
    sum(test_ER_ss_list(GytGy, alpha_a, mu_a, mu2_a))

  part3 <- sum(get_ER2_ss_list(GxtGx, GxtX, alpha_b, mu_b, mu2_b))

  part4 <- kl_gamma - sum(unlist(kl_a)) - sum(unlist(kl_b))

  elbo_curr <- part1 - (1/(2*sigma2_y))*part2 - (1/(2*sigma2_x))*part3 + part4

  return(list(elbo_curr, part2, part3))
}




