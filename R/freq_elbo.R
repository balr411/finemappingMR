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

