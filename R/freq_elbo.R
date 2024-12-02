#' @title Get the ELBO estimate
#'
#' @keywords internal

# Expected squared residuals.
# Xr is column sum of Xr_L
#Copying and modifying this function since we don't create a SuSiE object
get_ER2 <- function (X, Y, alpha, mu, mu2) {
  Xr_L = susieR:::compute_MXt(alpha * mu,X) # L by N matrix
  Xr = colSums(Xr_L)
  postb2 = alpha * mu2 # Posterior second moment.
  return(sum((Y - Xr)^2) - sum(Xr_L^2) + sum(attr(X,"d") * t(postb2)))
}

test_ER <- function (X, Y, alpha, mu, mu2) {
  Xr_L = susieR:::compute_MXt(alpha * mu,X) # L by N matrix
  Xr = colSums(Xr_L)
  postb2 = alpha * mu2 # Posterior second moment.
  return(sum(Xr^2)- sum(Xr_L^2) + sum(attr(X,"d") * t(postb2)))
}

freq_elbo <- function(sigma2_y, sigma2_x, y, x, Xr_a, Xr_b_y, G_x, G_y,
                      alpha_b, mu_b, mu2_b, alpha_a, mu_a, mu2_a, kl_a,
                      kl_b, mu_gamma){
  n_x <- nrow(G_x)
  n_y <- nrow(G_y)

  part1 <- -((n_x + n_y)/2) * log(2*pi) - (n_y/2)*log(sigma2_y) - (n_x/2)*log(sigma2_x)

  part2 <- crossprod(y) - 2*sum(y*Xr_a) - 2*sum(y*mu_gamma*Xr_b_y) +
    2*mu_gamma*sum(Xr_a*Xr_b_y) +
    mu_gamma^2*test_ER(G_y, y, alpha_b, mu_b, mu2_b) +
    test_ER(G_y, y, alpha_a, mu_a, mu2_a)

  part3 <- get_ER2(G_x, x, alpha_b, mu_b, mu2_b)

  part4 <- - sum(kl_a) - sum(kl_b)

  elbo_curr <- part1 - (1/(2*sigma2_y))*part2 - (1/(2*sigma2_x))*part3 + part4

  return(elbo_curr)
}

