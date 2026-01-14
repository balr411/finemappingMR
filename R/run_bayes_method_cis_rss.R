#' @title Bayesian MR estimate using summary statistics for exposure and outcome
#'  using Z-scores instead of sufficient statistics.
#'
#' @description Function that performs our Bayesian MR estimation method using
#' summary statistics for exposure and outcome.
#'
#' @param Z_x A vector of Z-scores for the exposure X, assuming centered and
#'  standardized G and x.
#'
#' @param Z_y A vector of Z-scores for the outcome y, assuming centered and
#'  standardized G and y.
#'
#' @param R A correlation matrix for the Z-scores. Assume it is the same for
#'  exposure and outcome studies, and that it is a positive-definite matrix. If
#'  your matrix is not positive-definite, project it to a near positive-definite
#'  matrix using Matrix::nearPD().
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
#' this is used as an initial value. This is fixed to 1 for the Z-score implementation.
#'
#' @param residual_variance_y The variance of the outcome residual. If \code{estimate_residual_variance_y = TRUE},
#' this is used as an initial value. This is fixed to 1 for the Z-score implementation.
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
#' to beta - gamma - alpha? Default = FALSE. Currently no functionality.
#'
#' @param include_pleiotropy Include a pleiotropy term in the estimation of gamma?
#' Default = TRUE. Currently no functionality.
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
#' @importFrom stringr str_glue
#'
#' @export

run_bayes_method_cis_rss <- function(Z_x, Z_y, R,
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

  #Project R to near positive semi-definite
  #R_near_pd <- Matrix::nearPD(R)
  #R <- R_near_pd$mat

  #Initialize all of the estimates
  #Right now just slightly hacking the original function but it should come out fine
  #Note that the scale_prior_variances are probably not correct. The put those to
  #50 or something in the SuSiE code?
  est_init <- fmr_init(M = 1, L_x, L_y, susie_init = NULL, susie_init_y = NULL,
                       scaled_prior_variance_x, scaled_prior_variance_y, varX = 1,
                       varY = 1, Gx_t_Gx = list(R), Gy_t_Gy = list(R), mu_gamma_init,
                       mu2_gamma_init)

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

  #Make extra copy of initial estimates for convergence checks
  b_post_old <- lapply(1:length(mu_b), function(x) rep(0, ncol(mu_b[[x]])))
  a_post_old <- lapply(1:length(mu_a), function(x) rep(0, ncol(mu_a[[x]])))
  mu_gamma_old <- mu_gamma

  elbo_conv_vec <- c()
  elbo_full_vec <- c()

  #Now start the main iteration loop
  conv <- FALSE
  iter <- 1

  #Get the inverse of R
  eig <- eigen(R, symmetric = TRUE)
  Rinv <- eig$vectors %*% diag(1 / eig$values) %*% t(eig$vectors)
  #RinvR <- Rinv %*% R
  #RRinvR <- rowSums(R * t(RinvR))
  RRinvR <- rep(1, nrow(Rinv))
  ZxRinvZx <- crossprod(Z_x, Rinv) %*% Z_x
  ZyRinvZy <- crossprod(Z_y, Rinv) %*% Z_y

  while(!conv & iter < max_iter){
    #Update alpha
    for(l in 1:L_y){
      if(V_y[l, 1] > 0){
        resid_al <- Z_y - sqrt(n_y) * mu_gamma * R %*% colSums(mu_b[[1]] * alpha_b[[1]]) - sqrt(n_y) * R %*% colSums(mu_a[[1]][-l, ] * alpha_a[[1]][-l, ])
        Rinv_resid_al <- Rinv %*% resid_al
        #Z_star_l <- sqrt(n_y) * R %*% Rinv_resid_al #Note this only works because R is symmetric
        Z_star_l <- sqrt(n_y) * resid_al

        #Update the prior variance
        res <- optim(par = V_y[l, 1],
                     fn = function(x) negloglik_rss_alpha(x, RRinvR, Z_star_l, n_y),
                     method = "Brent",
                     lower = 0,
                     upper = 1)

        if (!is.null(res$par) && res$convergence == 0) {
          V_y[l, 1] <- res$par
          #if (verbose) {
          #  cat(sprintf("Update s^2 for alpha effect %d to %f\n", l, V_y[l, 1]))
          #}
        } else {
          cat(sprintf("WARNING: s^2 alpha update for iteration %d, effect %d failed to converge; keeping previous parameters\n", iter, l))
        }

        if(V_y[l,1] > 0){
          omega_j <- n_y * RRinvR + 1/V_y[l, 1]
          #Update the PIPs and posterior means
          post_update <- update_posterior_rss(Z_star_l, omega_j)

          mu_a[[1]][l,] <- post_update$mu
          mu2_a[[1]][l,] <- post_update$mu2
          alpha_a[[1]][l,] <- post_update$alpha

          #Update the KL-divergence
          kl_a[[1]][l] <- get_ser_neg_KL_divergence_rss(pip = alpha_a[[1]][l,],
                                                        prior_var = V_y[l, 1],
                                                        post_mean = mu_a[[1]][l,],
                                                        post_var = mu2_a[[1]][l,] - mu_a[[1]][l,]^2)

        }else{
          mu_a[[1]][l,] <- 0
          mu2_a[[1]][l,] <- 0
          alpha_a[[1]][l,] <- 0
          kl_a[[1]][l] <- 0
          #What to do with KL-divergence??
          #Only use those that don't go to 0 at the end??
        }

      }

      elbo_full_vec <- c(elbo_full_vec, elbo_rss(ZxRinvZx, ZyRinvZy, Z_x, Z_y, mu_b, mu2_b, alpha_b,
                                                 mu_a, mu2_a, alpha_a, n_x, n_y, mu_gamma, mu2_gamma, R, Rinv,
                                                 kl_a, kl_b, kl_gamma))
    }

    #Update bl - note if there is some failure here it might revert to using the alpha updates since many of the names overlap
    for(l in 1:L_x){
      if(V_x[l, 1] > 0) {
        rblx <- Z_x - sqrt(n_x)* R %*% colSums(mu_b[[1]][-l,] * alpha_b[[1]][-l,])
        rbly <- mu_gamma*Z_y - sqrt(n_y) * mu2_gamma * R %*% colSums(mu_b[[1]][-l,] * alpha_b[[1]][-l,]) - sqrt(n_y) * mu_gamma * R %*% colSums(mu_a[[1]] * alpha_a[[1]])
        #resid_bl <- (1/sqrt(n_x + n_y*mu2_gamma)) * (sqrt(n_x) * rblx + sqrt(n_y) * rbly)

        #Rinv_resid_bl <- Rinv %*% resid_bl
        #Z_star_l <- sqrt(n_x + n_y*mu2_gamma) * R %*% Rinv_resid_bl #Note this only works because R is symmetric

        Z_star_l <- (sqrt(n_x) * rblx + sqrt(n_y) * rbly)

        #Update the prior variance
        res_x <- optim(par = V_x[l, 1],
                       fn = function(x) negloglik_rss_b(x, RRinvR, Z_star_l, n_x, n_y, mu2_gamma),
                       method = "Brent",
                       lower = 0,
                       upper = 1)

        if (!is.null(res_x$par) && res_x$convergence == 0) {
          V_x[l, 1] <- res_x$par
          #if (verbose) {
          #  cat(sprintf("Update s^2 for b effect %d to %f\n", l, V_x[l, 1]))
          #}
        } else {
          cat(sprintf("WARNING: s^2 update for b iteration %d, effect %d failed to converge; keeping previous parameters\n", iter, l))
        }

        if(V_x[l,1] > 0){
          omega_j <- (n_x + n_y*mu2_gamma) * RRinvR + 1/V_x[l, 1]

          #Update the PIPs and posterior means
          post_update <- update_posterior_rss(Z_star_l, omega_j)

          #post_update <- update_posterior_rss(Z_star_l, rep((n_x + n_y*mu2_gamma), length(RRinvR)) + 1/V_x[l, 1])

          mu_b[[1]][l,] <- post_update$mu
          mu2_b[[1]][l,] <- post_update$mu2
          alpha_b[[1]][l,] <- post_update$alpha

          #Update the KL-divergence
          kl_b[[1]][l] <- get_ser_neg_KL_divergence_rss(pip = alpha_b[[1]][l,],
                                                        prior_var = V_x[l, 1],
                                                        post_mean = mu_b[[1]][l,],
                                                        post_var = mu2_b[[1]][l,] - mu_b[[1]][l,]^2)

        }else{
          mu_b[[1]][l,] <- 0
          mu2_b[[1]][l,] <- 0
          alpha_b[[1]][l,] <- 0
          kl_b[[1]][l] <- 0
        }
      }

      elbo_full_vec <- c(elbo_full_vec, elbo_rss(ZxRinvZx, ZyRinvZy, Z_x, Z_y, mu_b, mu2_b, alpha_b,
                                                 mu_a, mu2_a, alpha_a, n_x, n_y, mu_gamma, mu2_gamma, R, Rinv,
                                                 kl_a, kl_b, kl_gamma))
    }

    #Now update gamma
    b_post <- lapply(Map("*", mu_b, alpha_b), colSums)
    a_post <- lapply(Map("*", mu_a, alpha_a), colSums)

    #Later for multiple regions we will have Z_y be a list of vectors, and the commented
    #out code will work.
    #mu_gamma_num <- sigma2_gamma_prior*sqrt(n_y)*(sum(unlist(Map("%*%", b_post, Z_y))) - sqrt(n_y)*sum(unlist(Map("%*%", Map("%*%", b_post, R), a_post))))

    mu_gamma_num <- sigma2_gamma_prior*sqrt(n_y)*(sum(b_post[[1]] %*% Z_y) - sqrt(n_y)*sum(unlist(b_post[[1]] %*% R %*% a_post[[1]])))

    vec_den <- vector(length = 1)
    for(m in 1:1){
      B <- alpha_b[[m]] * mu_b[[m]]
      XB2 <- n_y * sum((B %*% R) * B)
      betabar <- colSums(B)
      d <- n_y
      postb2 <- alpha_b[[m]] * mu2_b[[m]]
      vec_den[m] <- sum(betabar * (n_y*R %*% betabar)) - XB2 + sum(d * t(postb2))
    }

    mu_gamma_den <- 1 + sigma2_gamma_prior*(sum(vec_den))

    mu_gamma <- as.numeric(mu_gamma_num/mu_gamma_den)

    sigma2_gamma_curr <- as.numeric((sigma2_gamma_prior)/mu_gamma_den)

    mu2_gamma <- as.numeric(mu_gamma^2 + sigma2_gamma_curr)

    #Get KL-divergence of gamma (taken from https://stats.stackexchange.com/questions/7440/kl-divergence-between-two-univariate-gaussians):
    kl_gamma <- (1/2) - log(sqrt(sigma2_gamma_prior/sigma2_gamma_curr)) - (sigma2_gamma_curr + mu_gamma^2)/(2*sigma2_gamma_prior) #Should always be <0 (its negative kl divergence)

    #Now check the convergence based on parameter values
    #conv <- crossprod(b_post[[1]] - b_post_old[[1]]) + crossprod(a_post[[1]] - a_post_old[[1]]) + (mu_gamma - mu_gamma_old)^2 < tol
    #b_post_old <- b_post
    #a_post_old <- a_post
    #mu_gamma_old <- mu_gamma

    #Now calculate the ELBO and check convergence
    elbo_curr <- elbo_rss(ZxRinvZx, ZyRinvZy, Z_x, Z_y, mu_b, mu2_b, alpha_b,
                          mu_a, mu2_a, alpha_a, n_x, n_y, mu_gamma, mu2_gamma, R, Rinv,
                          kl_a, kl_b, kl_gamma)

    elbo_conv_vec <- c(elbo_conv_vec, elbo_curr)
    elbo_full_vec <- c(elbo_full_vec, elbo_curr)

    if(iter > 1){
      conv <- abs(elbo_conv_vec[iter] - elbo_conv_vec[iter-1]) < tol
    }else{
      conv <- FALSE
    }

    iter <- iter + 1

    if(verbose){
      print(str_glue("Starting iteration {iter}, gamma_curr = {mu_gamma}"))
      print(str_glue("ELBO: {elbo_conv_vec[iter - 1]}"))
    }

  }

  #Calculate the variance using the law of total variance
  total_var <- bayes_total_var_rss(V_x, V_y, mu_b, mu2_b, alpha_b, mu_a, mu2_a,
                                   alpha_a, R, Z_y, num_samples,
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

  #Calculate and return credible sets if desired
  if(calc_cs_x){
    to_return$cs_x <- get_cs_rss(V = V_x, alpha = alpha_b, R = list(R))
  }

  if(calc_cs_y){
    to_return$cs_y <- get_csrss(V = V_y, alpha = alpha_a, R = list(R))
  }

  return(to_return)

}

