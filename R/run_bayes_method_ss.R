#' @title Bayesian MR estimate using sufficient statistics for exposure and outcome.
#'
#' @description Function that performs our Bayesian MR estimation method using
#' sufficient statistics for exposure and outcome.
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
#'
#' @param calc_cs_y Calculate credible sets for the outcome? Default = FALSE.
#'
#' @param verbose Output progress? Default = FALSE.
#'
#' @param susie_init A list of SuSiE objects to initialize the exposure to. Default = NULL.
#'
#' @param susie_init_y A list of SuSiE objects to initialize the outcome to. Default = NULL.
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
#' @param beta_gamma_alpha Change the order of estimation from alpha - beta - gamma
#' to beta - gamma - alpha? Default = FALSE.
#'
#' @param include_pleiotropy Include a pleiotropy term in the estimation of gamma?
#' Default = TRUE.
#'
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

run_bayes_method_ss <- function(Gx_t_Gx, Gx_t_x, xtx,
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
                                beta_gamma_alpha = FALSE,
                                include_pleiotropy = TRUE){



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

  est_init <- fmr_init(M, L_x, L_y, susie_init, susie_init_y,
                       scaled_prior_variance_x, scaled_prior_variance_y, varX,
                       varY, Gx_t_Gx, Gy_t_Gy, mu_gamma_init, mu2_gamma_init)

  V_x <- est_init$V_x
  mu_b <- est_init$mu_b
  mu2_b <- est_init$mu2_b
  alpha_b <- est_init$alpha_b
  kl_b <- est_init$kl_b

  V_y <- est_init$V_y
  mu_a <- est_init$mu_a
  mu2_a <- est_init$mu2_a
  alpha_a <- est_init$alpha_a
  kl_a <- est_init$kl_a

  mu_gamma <- est_init$mu_gamma
  mu2_gamma <- est_init$mu2_gamma
  sigma2_gamma_curr <- est_init$sigma2_gamma_curr
  kl_gamma <- est_init$kl_gamma

  elbo_conv_vec <- c()
  lik_x <- c()
  lik_y <- c()

  #Initialize the d and scale for each of the GtG matrices
  for(i in 1:length(Gx_t_Gx)){
    csd = rep(1, ncol(Gx_t_Gx[[i]]))
    attr(Gx_t_Gx[[i]], "d") = diag(Gx_t_Gx[[i]])
    attr(Gx_t_Gx[[i]], "scaled:scale") = csd

    csd = rep(1, ncol(Gy_t_Gy[[i]]))
    attr(Gy_t_Gy[[i]], "d") = diag(Gy_t_Gy[[i]])
    attr(Gy_t_Gy[[i]], "scaled:scale") = csd
  }

  conv <- FALSE
  iter <- 1

  if(!beta_gamma_alpha & include_pleiotropy){
    while(!conv & iter < max_iter){
      #Update alpha first
      for(m in 1:M){
        for(l in 1:L_y){
          # Remove lth effect from fitted values.
          # Note that we shouldn't actually have to include each of the m loci since we are assuming that Gj^TGk = 0
          GytR <- Gy_t_y[[m]] - Gy_t_Gy[[m]] %*% (colSums(mu_a[[m]][-l,] * alpha_a[[m]][-l,]) + mu_gamma*colSums(mu_b[[m]] * alpha_b[[m]]))

          test_curr_ss <- susieR:::single_effect_regression_ss(as.matrix(GytR), attr(Gy_t_Gy[[m]],"d"), V_y[l, m], residual_variance = sigma2_y, optimize_V = "optim")

          mu_a[[m]][l,] <- test_curr_ss$mu
          mu2_a[[m]][l,] <- test_curr_ss$mu2
          alpha_a[[m]][l,] <- test_curr_ss$alpha
          V_y[l,m] <- test_curr_ss$V
          kl_a[[m]][l] <- -test_curr_ss$lbf_model +
            susieR:::SER_posterior_e_loglik_ss(attr(Gy_t_Gy[[m]],"d"), GytR, sigma2_y, test_curr_ss$alpha * test_curr_ss$mu,
                                               test_curr_ss$alpha * test_curr_ss$mu2)
        }
      }

      #Update b next
      for(m in 1:M){
        for(l in 1:L_x){
          #Get GtG
          Gstar_t_Gstar <- (1/sigma2_x)*Gx_t_Gx[[m]] + (mu2_gamma/sigma2_y)*Gy_t_Gy[[m]]
          csd = rep(1, length = ncol(Gx_t_Gx[[m]]))
          attr(Gstar_t_Gstar, "d") = diag(Gstar_t_Gstar)
          attr(Gstar_t_Gstar, "scaled:scale") = rep(1, length = nrow(Gstar_t_Gstar))

          Gstar_t_R <- (1/sigma2_x) * (Gx_t_x[[m]] -  Gx_t_Gx[[m]] %*% (colSums(mu_b[[m]][-l,] * alpha_b[[m]][-l,]))) + (1/sigma2_y)*(mu_gamma*Gy_t_y[[m]] - mu_gamma*Gy_t_Gy[[m]] %*% (colSums(mu_a[[m]] * alpha_a[[m]])) - mu2_gamma*Gy_t_Gy[[m]] %*% (colSums(mu_b[[m]][-l,] * alpha_b[[m]][-l,])))

          test_curr_ss <- susieR:::single_effect_regression_ss(as.matrix(Gstar_t_R), attr(Gstar_t_Gstar,"d"), V_x[l, m], residual_variance = 1, optimize_V = "optim")

          mu_b[[m]][l,] <- test_curr_ss$mu
          mu2_b[[m]][l,] <- test_curr_ss$mu2
          alpha_b[[m]][l,] <- test_curr_ss$alpha
          V_x[l,m] <- test_curr_ss$V
          kl_b[[m]][l] <- -test_curr_ss$lbf_model +
            susieR:::SER_posterior_e_loglik_ss(attr(Gstar_t_Gstar,"d"), Gstar_t_R, 1, test_curr_ss$alpha * test_curr_ss$mu,
                                               test_curr_ss$alpha * test_curr_ss$mu2)
        }
      }

      #Now update gamma
      b_post <- lapply(Map("*", mu_b, alpha_b), colSums)
      a_post <- lapply(Map("*", mu_a, alpha_a), colSums)

      mu_gamma_num <- sigma2_gamma_prior*(sum(unlist(Map("%*%", b_post, Gy_t_y))) - sum(unlist(Map("%*%", Map("%*%", b_post, Gy_t_Gy), a_post))))

      vec_den <- vector(length = M)
      for(m in 1:M){
        B <- alpha_b[[m]] * mu_b[[m]]
        XB2 <- sum((B %*% Gy_t_Gy[[m]]) * B)
        betabar <- colSums(B)
        d <- attr(Gy_t_Gy[[m]],"d")
        postb2 <- alpha_b[[m]] * mu2_b[[m]]
        vec_den[m] <- sum(betabar * (Gy_t_Gy[[m]] %*% betabar)) - XB2 + sum(d * t(postb2))
      }

      mu_gamma_den <- sigma2_y + sigma2_gamma_prior*(sum(vec_den))

      mu_gamma <- as.numeric(mu_gamma_num/mu_gamma_den)

      sigma2_gamma_curr <- as.numeric((sigma2_gamma_prior*sigma2_y)/mu_gamma_den)

      mu2_gamma <- as.numeric(mu_gamma^2 + sigma2_gamma_curr)

      #Get KL-divergence of gamma (taken from https://stats.stackexchange.com/questions/7440/kl-divergence-between-two-univariate-gaussians):
      kl_gamma <- (1/2) - log(sqrt(sigma2_gamma_prior/sigma2_gamma_curr)) - (sigma2_gamma_curr + mu_gamma^2)/(2*sigma2_gamma_prior) #Should always be <0 (its negative kl divergence)

      elbo_full <- bayes_elbo_ss(sigma2_y, sigma2_x, Gx_t_Gx, Gy_t_Gy, Gx_t_x, Gy_t_y, alpha_b, mu_b, mu2_b,
                                 alpha_a, mu_a, mu2_a, mu_gamma, mu2_gamma, kl_a, kl_b, kl_gamma, n_x, n_y)

      elbo_curr <- elbo_full[[1]]
      lik_y_curr <- elbo_full[[2]]
      lik_x_curr <- elbo_full[[3]]

      elbo_conv_vec <- c(elbo_conv_vec, elbo_curr)
      lik_x <- c(lik_x, lik_x_curr)
      lik_y <- c(lik_y, lik_y_curr)


      if(iter > 1){
        conv <- abs(elbo_conv_vec[iter] - elbo_conv_vec[iter-1]) < 1e-4
      }else{
        conv <- FALSE
      }

      iter <- iter + 1

      if(verbose){
        print(str_glue("Starting iteration {iter}, gamma_curr = {mu_gamma}"))
      }

    }

  }else if(beta_gamma_alpha & include_pleiotropy){
    while(!conv & iter < max_iter){
      #Update b first
      for(m in 1:M){
        for(l in 1:L_x){
          #Get GtG
          Gstar_t_Gstar <- (1/sigma2_x)*Gx_t_Gx[[m]] + (mu2_gamma/sigma2_y)*Gy_t_Gy[[m]]
          csd = rep(1, length = ncol(Gx_t_Gx[[m]]))
          attr(Gstar_t_Gstar, "d") = diag(Gstar_t_Gstar)
          attr(Gstar_t_Gstar, "scaled:scale") = rep(1, length = nrow(Gstar_t_Gstar))

          Gstar_t_R <- (1/sigma2_x) * (Gx_t_x[[m]] -  Gx_t_Gx[[m]] %*% (colSums(mu_b[[m]][-l,] * alpha_b[[m]][-l,]))) + (1/sigma2_y)*(mu_gamma*Gy_t_y[[m]] - mu_gamma*Gy_t_Gy[[m]] %*% (colSums(mu_a[[m]] * alpha_a[[m]])) - mu2_gamma*Gy_t_Gy[[m]] %*% (colSums(mu_b[[m]][-l,] * alpha_b[[m]][-l,])))

          test_curr_ss <- susieR:::single_effect_regression_ss(as.matrix(Gstar_t_R), attr(Gstar_t_Gstar,"d"), V_x[l, m], residual_variance = 1, optimize_V = "optim")

          mu_b[[m]][l,] <- test_curr_ss$mu
          mu2_b[[m]][l,] <- test_curr_ss$mu2
          alpha_b[[m]][l,] <- test_curr_ss$alpha
          V_x[l,m] <- test_curr_ss$V
          kl_b[[m]][l] <- -test_curr_ss$lbf_model +
            susieR:::SER_posterior_e_loglik_ss(attr(Gstar_t_Gstar,"d"), Gstar_t_R, 1, test_curr_ss$alpha * test_curr_ss$mu,
                                               test_curr_ss$alpha * test_curr_ss$mu2)
        }
      }

      #Now update gamma
      b_post <- lapply(Map("*", mu_b, alpha_b), colSums)
      a_post <- lapply(Map("*", mu_a, alpha_a), colSums)

      mu_gamma_num <- sigma2_gamma_prior*(sum(unlist(Map("%*%", b_post, Gy_t_y))) - sum(unlist(Map("%*%", Map("%*%", b_post, Gy_t_Gy), a_post))))

      vec_den <- vector(length = M)
      for(m in 1:M){
        B <- alpha_b[[m]] * mu_b[[m]]
        XB2 <- sum((B %*% Gy_t_Gy[[m]]) * B)
        betabar <- colSums(B)
        d <- attr(Gy_t_Gy[[m]],"d")
        postb2 <- alpha_b[[m]] * mu2_b[[m]]
        vec_den[m] <- sum(betabar * (Gy_t_Gy[[m]] %*% betabar)) - XB2 + sum(d * t(postb2))
      }

      mu_gamma_den <- sigma2_y + sigma2_gamma_prior*(sum(vec_den))

      mu_gamma <- as.numeric(mu_gamma_num/mu_gamma_den)

      sigma2_gamma_curr <- as.numeric((sigma2_gamma_prior*sigma2_y)/mu_gamma_den)

      mu2_gamma <- as.numeric(mu_gamma^2 + sigma2_gamma_curr)

      #Get KL-divergence of gamma (taken from https://stats.stackexchange.com/questions/7440/kl-divergence-between-two-univariate-gaussians):
      kl_gamma <- (1/2) - log(sqrt(sigma2_gamma_prior/sigma2_gamma_curr)) - (sigma2_gamma_curr + mu_gamma^2)/(2*sigma2_gamma_prior) #Should always be <0 (its negative kl divergence)

      #Update alpha
      for(m in 1:M){
        for(l in 1:L_y){
          # Remove lth effect from fitted values.
          # Note that we shouldn't actually have to include each of the m loci since we are assuming that Gj^TGk = 0
          GytR <- Gy_t_y[[m]] - Gy_t_Gy[[m]] %*% (colSums(mu_a[[m]][-l,] * alpha_a[[m]][-l,]) + mu_gamma*colSums(mu_b[[m]] * alpha_b[[m]]))

          test_curr_ss <- susieR:::single_effect_regression_ss(as.matrix(GytR), attr(Gy_t_Gy[[m]],"d"), V_y[l, m], residual_variance = sigma2_y, optimize_V = "optim")

          mu_a[[m]][l,] <- test_curr_ss$mu
          mu2_a[[m]][l,] <- test_curr_ss$mu2
          alpha_a[[m]][l,] <- test_curr_ss$alpha
          V_y[l,m] <- test_curr_ss$V
          kl_a[[m]][l] <- -test_curr_ss$lbf_model +
            susieR:::SER_posterior_e_loglik_ss(attr(Gy_t_Gy[[m]],"d"), GytR, sigma2_y, test_curr_ss$alpha * test_curr_ss$mu,
                                               test_curr_ss$alpha * test_curr_ss$mu2)
        }
      }

      elbo_full <- bayes_elbo_ss(sigma2_y, sigma2_x, Gx_t_Gx, Gy_t_Gy, Gx_t_x, Gy_t_y, alpha_b, mu_b, mu2_b,
                                 alpha_a, mu_a, mu2_a, mu_gamma, mu2_gamma, kl_a, kl_b, kl_gamma, n_x, n_y)

      elbo_curr <- elbo_full[[1]]
      lik_y_curr <- elbo_full[[2]]
      lik_x_curr <- elbo_full[[3]]

      elbo_conv_vec <- c(elbo_conv_vec, elbo_curr)
      lik_x <- c(lik_x, lik_x_curr)
      lik_y <- c(lik_y, lik_y_curr)


      if(iter > 1){
        conv <- abs(elbo_conv_vec[iter] - elbo_conv_vec[iter-1]) < 1e-4
      }else{
        conv <- FALSE
      }

      iter <- iter + 1

      if(verbose){
        print(str_glue("Starting iteration {iter}, gamma_curr = {mu_gamma}"))
      }

    }
  }else{
    #This is the case where we don't estimate pleiotropy. Let us keep it as a fixed
    #value (i.e. 0) for now and not update it, since this should be equivalent to setting it to 0

    #I will also set all prior variance estimates to 0 so that the total variance
    #calculation will work
    V_y_init <- 0
    V_y <- matrix(rep(V_y_init, L_y*M), nrow = L_y, ncol = M)

    #Update b first
    for(m in 1:M){
      for(l in 1:L_x){
        #Get GtG
        Gstar_t_Gstar <- (1/sigma2_x)*Gx_t_Gx[[m]] + (mu2_gamma/sigma2_y)*Gy_t_Gy[[m]]
        csd = rep(1, length = ncol(Gx_t_Gx[[m]]))
        attr(Gstar_t_Gstar, "d") = diag(Gstar_t_Gstar)
        attr(Gstar_t_Gstar, "scaled:scale") = rep(1, length = nrow(Gstar_t_Gstar))

        Gstar_t_R <- (1/sigma2_x) * (Gx_t_x[[m]] -  Gx_t_Gx[[m]] %*% (colSums(mu_b[[m]][-l,] * alpha_b[[m]][-l,]))) + (1/sigma2_y)*(mu_gamma*Gy_t_y[[m]] - mu_gamma*Gy_t_Gy[[m]] %*% (colSums(mu_a[[m]] * alpha_a[[m]])) - mu2_gamma*Gy_t_Gy[[m]] %*% (colSums(mu_b[[m]][-l,] * alpha_b[[m]][-l,])))

        test_curr_ss <- susieR:::single_effect_regression_ss(as.matrix(Gstar_t_R), attr(Gstar_t_Gstar,"d"), V_x[l, m], residual_variance = 1, optimize_V = "optim")

        mu_b[[m]][l,] <- test_curr_ss$mu
        mu2_b[[m]][l,] <- test_curr_ss$mu2
        alpha_b[[m]][l,] <- test_curr_ss$alpha
        V_x[l,m] <- test_curr_ss$V
        kl_b[[m]][l] <- -test_curr_ss$lbf_model +
          susieR:::SER_posterior_e_loglik_ss(attr(Gstar_t_Gstar,"d"), Gstar_t_R, 1, test_curr_ss$alpha * test_curr_ss$mu,
                                             test_curr_ss$alpha * test_curr_ss$mu2)
      }
    }

    #Now update gamma
    b_post <- lapply(Map("*", mu_b, alpha_b), colSums)
    a_post <- lapply(Map("*", mu_a, alpha_a), colSums)

    mu_gamma_num <- sigma2_gamma_prior*(sum(unlist(Map("%*%", b_post, Gy_t_y))) - sum(unlist(Map("%*%", Map("%*%", b_post, Gy_t_Gy), a_post))))

    vec_den <- vector(length = M)
    for(m in 1:M){
      B <- alpha_b[[m]] * mu_b[[m]]
      XB2 <- sum((B %*% Gy_t_Gy[[m]]) * B)
      betabar <- colSums(B)
      d <- attr(Gy_t_Gy[[m]],"d")
      postb2 <- alpha_b[[m]] * mu2_b[[m]]
      vec_den[m] <- sum(betabar * (Gy_t_Gy[[m]] %*% betabar)) - XB2 + sum(d * t(postb2))
    }

    mu_gamma_den <- sigma2_y + sigma2_gamma_prior*(sum(vec_den))

    mu_gamma <- as.numeric(mu_gamma_num/mu_gamma_den)

    sigma2_gamma_curr <- as.numeric((sigma2_gamma_prior*sigma2_y)/mu_gamma_den)

    mu2_gamma <- as.numeric(mu_gamma^2 + sigma2_gamma_curr)

    #Get KL-divergence of gamma (taken from https://stats.stackexchange.com/questions/7440/kl-divergence-between-two-univariate-gaussians):
    kl_gamma <- (1/2) - log(sqrt(sigma2_gamma_prior/sigma2_gamma_curr)) - (sigma2_gamma_curr + mu_gamma^2)/(2*sigma2_gamma_prior) #Should always be <0 (its negative kl divergence)

    elbo_full <- bayes_elbo_ss(sigma2_y, sigma2_x, Gx_t_Gx, Gy_t_Gy, Gx_t_x, Gy_t_y, alpha_b, mu_b, mu2_b,
                               alpha_a, mu_a, mu2_a, mu_gamma, mu2_gamma, kl_a, kl_b, kl_gamma, n_x, n_y)

    elbo_curr <- elbo_full[[1]]
    lik_y_curr <- elbo_full[[2]]
    lik_x_curr <- elbo_full[[3]]

    elbo_conv_vec <- c(elbo_conv_vec, elbo_curr)
    lik_x <- c(lik_x, lik_x_curr)
    lik_y <- c(lik_y, lik_y_curr)


    if(iter > 1){
      conv <- abs(elbo_conv_vec[iter] - elbo_conv_vec[iter-1]) < 1e-4
    }else{
      conv <- FALSE
    }

    iter <- iter + 1

    if(verbose){
      print(str_glue("Starting iteration {iter}, gamma_curr = {mu_gamma}"))
    }

  }

  #Calculate the variance using the law of total variance
  total_var <- bayes_total_var(V_x, V_y, mu_b, mu2_b, alpha_b, mu_a, mu2_a,
                               alpha_a, Gy_t_Gy, Gy_t_y, num_samples, sigma2_y,
                               sigma2_gamma_prior)


  #Return a list with the various components
  to_return <- list()

  #Gamma estimation
  to_return$res <- data.frame(gamma = mu_gamma,
                              gamma_var = sigma2_gamma_curr,
                              gamma_total_var = total_var,
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
  to_return$lik_x <- lik_x
  to_return$lik_y <- lik_y

  #KL-divergence
  to_return$kl_a <- kl_a
  to_return$kl_b <- kl_b

  #Calculate and return credible sets if desired
  if(calc_cs_x){
    to_return$cs_x <- get_cs(V = V_x, alpha = alpha_b, G_t_G = Gx_t_Gx)
  }

  if(calc_cs_y){
    to_return$cs_y <- get_cs(V = V_y, alpha = alpha_a, G_t_G = Gy_t_Gy)
  }

  return(to_return)


}



