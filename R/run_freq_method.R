#' @title Frequentist MR estimate
#'
#' @description Function that performs our frequentist MR estimation method.
#'
#' @param x Exposure individual-level phenotype data.
#'
#' @param y Outcome individual-level phenotype data.
#'
#' @param G_x Exposure genotype matrix.
#'
#' @param G_y Outcome genotype matrix.
#'
#' @param L_x Maximum number of non-zero effects for the exposure.
#'
#' @param L_y Maximum number of non-zero effects for the outcome.
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
#' @return A list containing various results from the estimation procedure,
#' including a data frame containing the gamma estimation results, as well as
#' the posterior first and second moments for b and alpha, the prior non-zero
#' effect variance estimates, and the PIPs for each non-zero effect, and the
#' vector of ELBO iterations.
#'
#' @export

run_freq_method <- function(x, y, G_x, G_y, L_x = 10, L_y = 10,
                            scaled_prior_variance_x = 0.2,
                            scaled_prior_variance_y = 0.2,
                            estimate_prior_variance_x = TRUE,
                            estimate_prior_variance_y = TRUE,
                            residual_variance_x = NULL,
                            residual_variance_y = NULL,
                            estimate_residual_variance_x = FALSE,
                            estimate_residual_variance_y = FALSE,
                            tol = 1e-4,
                            max_iter = 1000){

  if(is.null(residual_variance_x)){
    residual_variance_x <- var(x)
  }

  if(is.null(residual_variance_y)){
    residual_variance_y <- var(y)
  }

  sigma2_x <- residual_variance_x
  sigma2_y <- residual_variance_y

  #Initialize variance estimates
  V_x_init <- scaled_prior_variance_x*var(x)
  V_y_init <- scaled_prior_variance_y*var(y)

  V_x <- rep(V_x_init, L_x)
  V_y <- rep(V_y_init, L_y)

  #Center x and y - note in original SuSiE code this is equivalent to intercept = TRUE. Going to force this to be true for now.
  x <- x - mean(x)
  y <- y - mean(y)

  #Center G_x, G_y
  out_x = susieR:::compute_colstats(G_x, center = TRUE, scale = FALSE)
  G_x <- sweep(G_x, 2, out_x$cm) #This assumes intercept = TRUE in original SuSiE code. Going to force this to be true for now.
  #G_x <- sweep(G_x, 2, out_x$csd, "/") - Note: this would scale the columns of G_x (i.e. standardize = TRUE in original SuSiE code)

  out_y = susieR:::compute_colstats(G_y, center = TRUE, scale = FALSE)
  G_y <- sweep(G_y, 2, out_y$cm) #This assumes intercept = TRUE in original SuSiE code. Going to force this to be true for now.
  #G_y <- sweep(G_y, 2, out_y$csd, "/") - Note: this would scale the columns of G_y (i.e. standardize = TRUE in original SuSiE code)

  #Fitted values for alpha estimation
  Xr_a <- rep(0, nrow(G_y))

  #b estimates
  mu_b <- matrix(0, L_x, ncol(G_x))
  mu2_b <- matrix(0, L_x, ncol(G_x))
  #Note that alpha_b is the PIP for the exposure variants, not the alpha parameter in our model. Should eventually change notation.
  alpha_b <- matrix(1/ncol(G_x), L_x, ncol(G_x)) #this is equivalent to setting prior_weights = 1/number of variants in SuSiE
  kl_b <- vector(length = L_x)

  #alpha (a) estimates
  mu_a <- matrix(0, L_y, ncol(G_y))
  mu2_a <- matrix(0, L_y, ncol(G_y))
  #Note that alpha_a is the PIP for the outcome variants, not the alpha parameter in our model. Should eventually change notation.
  alpha_a <- matrix(1/ncol(G_y), L_y, ncol(G_y)) #this is equivalent to setting prior_weights = 1/number of variants in SuSiE
  kl_a <- vector(length = L_y)

  #gamma estimate
  sigma2_gamma_curr <- 0
  mu_gamma <- 0

  #Note that the columns of G_x and G_y have been manually centered. Hence we
  #set intercept = FALSE here.
  standardize = FALSE
  intercept = FALSE
  out_x = susieR:::compute_colstats(G_x, center = intercept, scale = standardize)

  attr(G_x,"scaled:center") = out_x$cm
  attr(G_x,"scaled:scale") = out_x$csd
  attr(G_x,"d") = out_x$d

  out_y = susieR:::compute_colstats(G_y, center = intercept, scale = standardize)

  attr(G_y,"scaled:center") = out_y$cm
  attr(G_y,"scaled:scale") = out_y$csd
  attr(G_y,"d") = out_y$d

  elbo_conv_vec <- c()

  iter <- 1
  conv <- FALSE

  #G_star matrix for b updates
  G_star_x <- (1/(sqrt(sigma2_x)))*G_x

  while(!conv & iter < max_iter){
    #Re-initialize beta part of the residuals with updated gamma and beta values
    Xr_b_y <- mu_gamma*colSums(susieR:::compute_MXt(alpha_b*mu_b, G_y))

    #Update alpha first
    for(l in 1:L_y){
      # Remove lth effect from fitted values.
      Xr_a = Xr_a - susieR:::compute_Xb(G_y, alpha_a[l,] * mu_a[l,])

      R = y - Xr_a - Xr_b_y
      test_curr <- susieR:::single_effect_regression(R, G_y, V_y[l], residual_variance = sigma2_y, optimize_V = "optim")

      mu_a[l,] <- test_curr$mu
      mu2_a[l,] <- test_curr$mu2
      alpha_a[l,] <- test_curr$alpha
      V_y[l] <- test_curr$V
      kl_a[l] <- -test_curr$loglik +
        susieR:::SER_posterior_e_loglik(G_y, R, sigma2_y, test_curr$alpha * test_curr$mu,
                                        test_curr$alpha * test_curr$mu2)

      Xr_a = Xr_a + susieR:::compute_Xb(G_y, alpha_a[l,] * mu_a[l,])
    }

    #Update b next
    for(l in 1:10){
      # Update beta part of residuals
      Xr_b_x <- colSums(susieR:::compute_MXt(alpha_b[-l,]*mu_b[-l,], G_x))
      Xr_b_y <- colSums(susieR:::compute_MXt(alpha_b[-l,]*mu_b[-l,], G_y))

      R_1 <- (1/(sqrt(sigma2_x)))*(x - Xr_b_x)
      R_2 <- (1/(sqrt(sigma2_y)))*(y - Xr_a - mu_gamma * Xr_b_y)

      R <- c(R_1, R_2)

      G_star <- rbind(G_star_x, (mu_gamma/sqrt(sigma2_y))*G_y)

      out_star = susieR:::compute_colstats(G_star, center = intercept, scale = standardize)

      attr(G_star,"scaled:center") = out_star$cm
      attr(G_star,"scaled:scale") = out_star$csd
      attr(G_star,"d") = out_star$d

      resid_var_curr <- 1

      test_curr <- susieR:::single_effect_regression(R, G_star, V_x[l], residual_variance = resid_var_curr, optimize_V = "optim")

      mu_b[l,] <- test_curr$mu
      mu2_b[l,] <- test_curr$mu2
      alpha_b[l,] <- test_curr$alpha
      V_x[l] <- test_curr$V
      kl_b[l] <- -test_curr$loglik +
        susieR:::SER_posterior_e_loglik(G_star, R, resid_var_curr, test_curr$alpha * test_curr$mu,
                                        test_curr$alpha * test_curr$mu2)
    }

    elbo_curr <- freq_elbo(sigma2_y, sigma2_x, y, x, Xr_a,
                           Xr_b_y = colSums(susieR:::compute_MXt(alpha_b*mu_b, G_y)),
                           G_x, G_y, alpha_b, mu_b, mu2_b, alpha_a, mu_a, mu2_a, kl_a,
                           kl_b, mu_gamma) #Do this before updating gamma since previous gamma was used in the update for kl(a) and kl(b)

    elbo_conv_vec <- c(elbo_conv_vec, elbo_curr)

    #Now update gamma
    #Update b part of residuals
    Xr_b_y <- colSums(susieR:::compute_MXt(alpha_b*mu_b, G_y))

    #Should give E[sum Gbl^T sum Gbl]
    sum_transpose <- get_ER2(G_y, y, alpha_b, mu_b, mu2_b) - crossprod(y) + 2*sum(y*Xr_b_y)

    R_gamma <- y - Xr_a
    R_Eb_t <- sum(R_gamma*Xr_b_y) # Should be equal to t(y-Xr_a) %*% Xr_b_y

    mu_gamma_num <- R_Eb_t
    mu_gamma_den <- sum_transpose

    mu_gamma <- as.numeric(mu_gamma_num/mu_gamma_den)

    sigma2_gamma_curr <- as.numeric((sigma2_y)/mu_gamma_den)

    if(iter > 1){
      conv <- abs(elbo_conv_vec[iter] - elbo_conv_vec[iter-1]) < 1e-4
    }else{
      conv <- FALSE
    }

    iter <- iter + 1
  }

  if(iter >= max_iter){
    warning("Maximum iteration number reached, convergence not attained")
  }

  #Return a list with the various components
  to_return <- list()

  #Gamma estimation
  to_return$res <- data.frame(gamma = mu_gamma,
                              gamma_var = sigma2_gamma_curr,
                              iter = iter)

  #b estimation
  to_return$mu_b <- mu_b
  to_return$mu2_b <- mu2_b
  to_return$alpha_b <- alpha_b
  to_return$V_x <- V_x

  #alpha estimation
  to_return$mu_a <- mu_a
  to_return$mu2_a <- mu2_a
  to_return$alpha_a <- alpha_a
  to_return$V_y <- V_y

  #ELBO
  to_return$elbo <- elbo_conv_vec

}



