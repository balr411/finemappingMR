#' @title Bayesian finemappingMR-inf estimate using sufficient statistics for exposure and outcome.
#'
#' @description Function that performs our Bayesian MR estimation method using
#' sufficient statistics for exposure and outcome while incorporating an infinitesimal
#' term for both the exposure and outcome.
#'
#' @param Gx_t_Gx A list of matrices \eqn{G_x'G_x} in which the columns of \eqn{G_x}
#'   are centered to have mean zero. Exposure data.
#'
#' @param Gx_t_x A list of matrices \eqn{G_x'x} in which the columns of \eqn{G_x}
#'   are centered to have mean zero. Exposure data.
#'
#' @param xtx A scalar \eqn{x'x} in which x is centered to have mean 0. Exposure data.
#'
#' @param Gy_t_Gy A list of matrices \eqn{G_y'G_y} in which the columns of \eqn{G_y}
#'   are centered to have mean zero. Outcome data.
#'
#' @param Gy_t_y A list of matrices \eqn{G_y'y} in which the columns of \eqn{G_y}
#'   are centered to have mean zero. Outcome data.
#'
#' @param yty A scalar \eqn{y'y} in which y is centered to have mean 0. Outcome data.
#'
#' @param L_x Maximum number of non-zero effects for the exposure.
#'
#' @param L_y Maximum number of non-zero effects for the outcome.
#'
#' @param n_x The exposure sample size.
#'
#' @param n_y The outcome sample size.
#'
#' @param scaled_prior_variance_x The prior variance, divided by
#'   \code{var(x)}; that is, the prior variance of each
#'   non-zero element of b is \code{var(x) * scaled_prior_variance_x}. The
#'   value provided should be either a scalar or a vector of length
#'   \code{L_x}. If \code{estimate_prior_variance_x = TRUE}, this provides
#'   initial estimates of the prior variances.
#'
#' @param scaled_prior_variance_y The prior variance, divided by
#'   \code{var(y)}; that is, the prior variance of each
#'   non-zero element of alpha is \code{var(y) * scaled_prior_variance_y}. The
#'   value provided should be either a scalar or a vector of length
#'   \code{L_x}. If \code{estimate_prior_variance_y = TRUE}, this provides
#'   initial estimates of the prior variances.
#'
#' @param estimate_prior_variance_x If \code{estimate_prior_variance_x =
#'   TRUE}, the prior variance is estimated (this is a separate
#'   parameter for each of the L_x effects). If provided,
#'   \code{scaled_prior_variance_x} is then used as an initial value for
#'   the optimization which is done using \code{optim}.
#'   When \code{estimate_prior_variance_x = FALSE}, the
#'   prior variance for each of the L_x effects is determined by the
#'   value supplied to \code{scaled_prior_variance_x}.
#'
#' @param estimate_prior_variance_y If \code{estimate_prior_variance_y =
#'   TRUE}, the prior variance is estimated (this is a separate
#'   parameter for each of the L_y effects). If provided,
#'   \code{scaled_prior_variance_y} is then used as an initial value for
#'   the optimization which is done using \code{optim}.
#'   When \code{estimate_prior_variance_y = FALSE}, the
#'   prior variance for each of the L_y effects is determined by the
#'   value supplied to \code{scaled_prior_variance_y}. Currently has
#'   no functionality.
#'
#' @param residual_variance_x The variance of the exposure residual. If \code{estimate_residual_variance_x = TRUE},
#' this is used as an initial value.
#'
#' @param residual_variance_y The variance of the outcome residual. If \code{estimate_residual_variance_y = TRUE},
#' this is used as an initial value.
#'
#' @param estimate_residual_variance_x If \code{TRUE}, exposure residual variance is estimated
#' using initial value \code{residual_variance_x}. Currently has
#' no functionality.
#'
#' @param estimate_residual_variance_y If \code{TRUE}, outcome residual variance is estimated
#' using initial value \code{residual_variance_y}. Currently has
#' no functionality.
#'
#' @param tol A small, non-negative number specifying the convergence
#'   tolerance for the variational Bayes procedure. The fitting procedure
#'   will halt when the difference in the variational lower bound, or
#'   \dQuote{ELBO} (the objective function to be maximized), is
#'   less than \code{tol}.
#'
#' @param max_iter Maximum number of iterations before stopping estimation procedure.
#'
#' @param calc_cs_x Calculate credible sets for the exposure? Default = FALSE.
#' Currently has no function.
#'
#' @param calc_cs_y Calculate credible sets for the outcome? Default = FALSE.
#' Currently has no function.
#'
#' @param verbose Output progress? Default = FALSE.
#'
#' @param susie_init A list of SuSiE objects to initialize the exposure to. Default = NULL.
#' Currently has no function.
#'
#' @param susie_init_y A list of SuSiE objects to initialize the outcome to. Default = NULL.
#' Currently has no function.
#'
#' @param mu_gamma_init Value to initialize the posterior mean of gamma to.
#'  Default = 0.
#'
#' @param mu2_gamma_init Value to initialize the posterior second moment of gamma
#'  to. Default = 0.
#'
#' @param sigma2_gamma_prior Prior variance of gamma to use. Default = 10^6.
#'
#' @param num_samples Number of samples to use in the Monte Carlo sampling for
#'  the variance. Default = 1,000.
#'
#' @param tausq_x The variance of the infinitesimal term for the exposure. Later
#' will add functionality to estimate this value, and this variable will be used
#' as an initial estimate.
#'
#' @param tausq_y The variance of the infinitesimal term for the outcome. Later
#' will add functionality to estimate this value, and this variable will be used
#' as an initial estimate.
#'
#' @return A list containing various results from the estimation procedure,
#' including a data frame containing the gamma estimation results, as well as
#' the posterior first and second moments for b and alpha, the prior non-zero
#' effect variance estimates, and the PIPs for each non-zero effect, the
#' vector of ELBO iterations, the vector of likelihood components, and the KL
#' divergence terms for b and alpha.
#'
#' @importFrom stats var
#'
#' @export

finemappingMR_inf_ss <- function(Gx_t_Gx, Gx_t_x, xtx,
                                 Gy_t_Gy, Gy_t_y, yty,
                                 n_x, n_y,
                                 L_x = 10, L_y = 10,
                                 scaled_prior_variance_x = 0.2,
                                 scaled_prior_variance_y = 0.2,
                                 estimate_prior_variance_x = TRUE,
                                 estimate_prior_variance_y = TRUE,
                                 residual_variance_x = NULL,
                                 residual_variance_y = NULL,
                                 estimate_residual_variance_x = FALSE,
                                 estimate_residual_variance_y = FALSE,
                                 tol = 1e-4,
                                 max_iter = 1000,
                                 calc_cs_x = FALSE,
                                 calc_cs_y = FALSE,
                                 verbose = FALSE,
                                 susie_init = NULL,
                                 susie_init_y = NULL,
                                 mu_gamma_init = 0,
                                 mu2_gamma_init = 0,
                                 sigma2_gamma_prior = 10^6,
                                 num_samples = 1000,
                                 tausq_x = 1,
                                 tausq_y = 1){

  varX <- xtx/(n_x - 1) #Note need to think about changing this to add functionality for using summary statistics where xtx/yty are unknown
  varY <- yty/(n_y - 1)

  if(is.null(residual_variance_x)){
    residual_variance_x <- varX
  }

  if(is.null(residual_variance_y)){
    residual_variance_y <- varY
  }

  sigma2_x <- residual_variance_x
  sigma2_y <- residual_variance_y

  #Initialize all of the estimates
  M <- length(Gx_t_Gx) #Number of loci - should later add some checks that the same amount of regions were given

  #First perform the Eigen decompositions for each of the matrices
  #Note that in Alex Mccreight's code he used LD = GtG/n in the eigen value computation but I don't know why he would do that?
  #Also note that for PCSK9 on my first attempt at this, 894/5746 eigenvalues were negative,
  #which shouldn't happen
  #Also note that right now this will only work if M = 1
  for(i in 1:M){
    #Gx_t_Gx_nearPD <- Matrix::nearPD(Gx_t_Gx[[i]]) - it looked like for the first few entries that this
    #didn't change the actual matrix very much, and the largest eigenvalues were identical to the original ones,
    #but now the smaller ones are non-negative
    #Will talk to Hyun/William/Jean about this

    #Do for exposure first
    eig_x <- eigen(Gx_t_Gx[[i]], symmetric = TRUE)

    #Order decreasing - what is the consequence of doing this?
    idx_x <- order(eig_x$values, decreasing = TRUE)
    eig_x$values <- eig_x$values[idx_x]
    eig_x$vectors <- eig_x$vectors[, idx_x]
    V_eigen_x <- eig_x$vectors

    Dsq_x <- pmax(eig_x$values, 0) #Note that before the code was pmax(eig$values * n, 0), which
    #simultaneously puts the matrix back onto the correct scale and sets all of the negative
    #eigenvalues to 0. Since I decomposed the matrix without dividing by n, I do not
    #need to multiple back by n here (should check this though)

    #Define the other matrices that will be of use
    u_x <- t(V_eigen_x) %*% Gx_t_x[[i]]
    tauDx_sigma2x <- tausq_x*Dsq_x + sigma2_x
    diag_Gxt_Omega_Gx <- rowSums(sweep(V_eigen_x^2, 2, (Dsq_x / tauDx_sigma2x), `*`)) #Diagonal of the VD(tau^2D + sigma^2I)^-1 V^T matrix
    Gx_t_Omega_x <- V_eigen_x %*% (u_x / tauDx_sigma2x) #The Gx^T Omega_x x matrix - note that in the paper it is written V (tau^2D + sigma^2I)^-1 u but this should work

    #Now do for outcome
    eig_y <- eigen(Gy_t_Gy[[i]], symmetric = TRUE)

    #Order decreasing - what is the consequence of doing this?
    idx_y <- order(eig_y$values, decreasing = TRUE)
    eig_y$values <- eig_y$values[idx_y]
    eig_y$vectors <- eig_y$vectors[, idx_y]
    V_eigen_y <- eig_y$vectors

    Dsq_y <- pmax(eig_y$values, 0) #Note that before the code was pmax(eig$values * n, 0), which
    #simultaneously puts the matrix back onto the correct scale and sets all of the negative
    #eigenvalues to 0. Since I decomposed the matrix without dividing by n, I do not
    #need to multiple back by n here (should check this though)

    #Define the other matrices that will be of use
    u_y <- t(V_eigen_y) %*% Gy_t_y[[i]]
    tauDy_sigma2y <- tausq_y*Dsq_y + sigma2_y
    diag_Gyt_Omega_Gy <- rowSums(sweep(V_eigen_y^2, 2, (Dsq_y / tauDy_sigma2y), `*`)) #Diagonal of the VD(tau^2D + sigma^2I)^-1 V^T matrix
    Gy_t_Omega_y <- V_eigen_y %*% (u_y / tauDy_sigma2y) #The Gx^T Omega_x x matrix - note that in the paper it is written V (tau^2D + sigma^2I)^-1 u but this should work

  }

  #Now initialize the prior variances and other things
  #Note that later I will need to put these in lists/matrices for M > 1
  p <- nrow(Gx_t_Gx[[1]])
  ssq_b <- rep(scaled_prior_variance_x * varX, L_x)
  PIP_b <- matrix(1/p, nrow = p, ncol = L_x) #In other code this is alpha_b so probably change it after
  mu_b <- matrix(0, nrow = p, ncol = L_x)

  ssq_a <- rep(scaled_prior_variance_y * varY, L_y)
  PIP_a <- matrix(1/p, nrow = p, ncol = L_y) #In other code this is alpha_a so probably change it after
  mu_a <- matrix(0, nrow = p, ncol = L_y)

  #################################################
  #Stopped here; line 248 of Alex's code; Uses lbf_variable_b to update the pip
  #lbf_variable_b <- matrix(0, nrow = p, ncol = L_x)
  #lbf_b <- rep(0, L_x)
  #omega_x <- matrix(diagXtOmegaX, nrow = p, ncol = L_x) + 1 / ssq
  ################################################

  #Start loop
  conv <- FALSE
  iter <- 1

  while(!conv & iter < max_iter){
    #Update alpha first
    #Easiest to make a Zl and omega vector for each l




  }


}



