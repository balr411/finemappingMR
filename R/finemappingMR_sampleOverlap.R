#' @title Bayesian MR estimate using summary statistics for exposure and outcome
#'  using Z-scores, assuming there is correlation induced by sample overlap.
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
#'  exposure and outcome studies
#'
#' @param rho The correlation induced by sample overlap.
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
#'   less than \code{tol}. Currently we do not use the ELBO to check convergence
#'   but rather check the changes in the parameters
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

finemappingMR_sampleOverlap <- function(Z_x, Z_y, R, rho,
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

  p <- length(Z_x)

  elbo_conv_vec <- c()
  #elbo_full_vec <- c()

  #Now start the main iteration loop
  conv <- FALSE
  iter <- 1

  #Get the inverse of R
  eig <- eigen(R, symmetric = TRUE)
  Rinv <- eig$vectors %*% diag(1 / eig$values) %*% t(eig$vectors)
  RRinvR <- rep(1, nrow(Rinv))
  ZxRinvZx <- crossprod(Z_x, Rinv) %*% Z_x
  ZyRinvZy <- crossprod(Z_y, Rinv) %*% Z_y
  ZxRinvZy <- crossprod(Z_x, Rinv) %*% Z_y

  while(!conv & iter < max_iter){
    #Update alpha
    for(l in 1:L_y){
      resid_al <-  ((Z_y - sqrt(n_y) * mu_gamma * R %*% colSums(mu_b[[1]] * alpha_b[[1]]) - sqrt(n_y) * R %*% colSums(mu_a[[1]][-l, ] * alpha_a[[1]][-l, ])) + rho*Z_x - rho*sqrt(n_x) * R %*% colSums(mu_b[[1]] * alpha_b[[1]]))
      Z_star_l <- (sqrt(n_y)/(1 - rho^2)) * resid_al

      #Update the prior variance
      res <- optim(par = ifelse(V_y[l, 1] > 0, log(V_y[l, 1]), -30),
                   fn = function(x) negloglik_sampleOverlap_alpha_trans(x, RRinvR, Z_star_l, n_y, rho),
                   method = "Brent",
                   lower = -30,
                   upper = 0)

      if (!is.null(res$par) && res$convergence == 0) {
        #First check if the new parameter beats the old:
        if(V_y[l,1] > 0){
          if(negloglik_sampleOverlap_alpha_trans(res$par, RRinvR, Z_star_l, n_y, rho) > negloglik_sampleOverlap_alpha_trans(log(V_y[l,1]), RRinvR, Z_star_l, n_y, rho)){
            V_y_new <- V_y[l,1]
          }else{
            V_y_new <- exp(res$par)
          }
        }else{
          V_y_new <- exp(res$par)
        }


        #Now check if the parameter beats 0:
        if(log(p) > -negloglik_sampleOverlap_alpha_trans(log(V_y_new), RRinvR, Z_star_l, n_y, rho)){
          V_y[l, 1] <- 0
        }else{
          V_y[l, 1] <- V_y_new
        }
        #if (verbose) {
        #  cat(sprintf("Update s^2 for alpha effect %d to %f\n", l, V_y[l, 1]))
        #}
      } else {
        cat(sprintf("WARNING: s^2 alpha update for iteration %d, effect %d failed to converge; keeping previous parameters\n", iter, l))
      }

      if(V_y[l,1] > 0){
        omega_j <- (n_y/(1 - rho^2)) * RRinvR + 1/V_y[l, 1]

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


      #elbo_full_vec <- c(elbo_full_vec, elbo_sampleOverlap(ZxRinvZx, ZyRinvZy, ZxRinvZy, Z_x, Z_y, mu_b, mu2_b, alpha_b,
      #                                                     mu_a, mu2_a, alpha_a, n_x, n_y, mu_gamma, mu2_gamma, R, Rinv,
      #                                                     kl_a, kl_b, kl_gamma, rho))

    }

    #Update bl - note if there is some failure here it might revert to using the alpha updates since many of the names overlap
    for(l in 1:L_x){

      rblx <- Z_x - sqrt(n_x)* R %*% colSums(mu_b[[1]][-l,] * alpha_b[[1]][-l,])
      rbly <- Z_y - sqrt(n_y) * mu_gamma * R %*% colSums(mu_b[[1]][-l,] * alpha_b[[1]][-l,]) - sqrt(n_y) * R %*% colSums(mu_a[[1]] * alpha_a[[1]])

      r_gam_blx <- mu_gamma*Z_x - sqrt(n_x)* mu_gamma * R %*% colSums(mu_b[[1]][-l,] * alpha_b[[1]][-l,])
      r_gam_bly <- mu_gamma*Z_y - sqrt(n_y) * mu2_gamma * R %*% colSums(mu_b[[1]][-l,] * alpha_b[[1]][-l,]) - sqrt(n_y) * mu_gamma * R %*% colSums(mu_a[[1]] * alpha_a[[1]])

      resid_bl <- sqrt(n_x)*rblx + sqrt(n_x)*rho*rbly + sqrt(n_y)*rho*r_gam_blx + sqrt(n_y)*r_gam_bly

      Z_star_l <- (1/(1 - rho^2)) * resid_bl

      #Update the prior variance
      res_x <- optim(par = ifelse(V_x[l, 1] > 0, log(V_x[l, 1]), -30),
                     fn = function(x) negloglik_sampleOverlap_b_trans(x, RRinvR, Z_star_l, n_x, n_y, mu_gamma, mu2_gamma, rho),
                     method = "Brent",
                     lower = -30,
                     upper = 0)

      if (!is.null(res_x$par) && res_x$convergence == 0) {
        #First check if the new parameter beats the old:
        if(V_x[l,1] > 0){
          if(negloglik_sampleOverlap_b_trans(res_x$par, RRinvR, Z_star_l, n_x, n_y, mu_gamma, mu2_gamma, rho) > negloglik_sampleOverlap_b_trans(log(V_x[l,1]),  RRinvR, Z_star_l, n_x, n_y, mu_gamma, mu2_gamma, rho)){
            V_x_new <- V_x[l,1]
          }else{
            V_x_new <- exp(res_x$par)
          }
        }else{
          V_x_new <- exp(res_x$par)
        }

        #Now check if the parameter beats 0:
        if(log(p) > -negloglik_sampleOverlap_b_trans(log(V_x_new), RRinvR, Z_star_l, n_x, n_y, mu_gamma, mu2_gamma, rho)){
          V_x[l, 1] <- 0
        }else{
          V_x[l, 1] <- V_x_new
        }

        #if (verbose) {
        #  cat(sprintf("Update s^2 for b effect %d to %f\n", l, V_x[l, 1]))
        #}
      } else {
        cat(sprintf("WARNING: s^2 update for b iteration %d, effect %d failed to converge; keeping previous parameters\n", iter, l))
      }

      if(V_x[l,1] > 0){
        omega_j <- ((n_x + 2*rho*sqrt(n_x*n_y)*mu_gamma + n_y*mu2_gamma)/(1-rho^2)) * RRinvR + 1/V_x[l, 1]

        #Update the PIPs and posterior means
        post_update <- update_posterior_rss(Z_star_l, omega_j) #Getting NA in the alpha

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


      #elbo_full_vec <- c(elbo_full_vec, elbo_sampleOverlap(ZxRinvZx, ZyRinvZy, ZxRinvZy, Z_x, Z_y, mu_b, mu2_b, alpha_b,
      #                                                     mu_a, mu2_a, alpha_a, n_x, n_y, mu_gamma, mu2_gamma, R, Rinv,
      #                                                     kl_a, kl_b, kl_gamma, rho))
    }

    #Now update gamma
    b_post <- lapply(Map("*", mu_b, alpha_b), colSums)
    a_post <- lapply(Map("*", mu_a, alpha_a), colSums)

    #Later for multiple regions we will have Z_y be a list of vectors, and the commented
    #out code will work.
    #mu_gamma_num <- sigma2_gamma_prior*sqrt(n_y)*(sum(unlist(Map("%*%", b_post, Z_y))) - sqrt(n_y)*sum(unlist(Map("%*%", Map("%*%", b_post, R), a_post))))

    E_btRb <- b_t_b(mu_b, mu2_b, alpha_b, R, n = 1) #Will need to modify for multiple regions

    mu_gamma_den <- 1 - rho^2 + sigma2_gamma_prior*n_y*E_btRb

    mu_gamma_num_part1 <- sqrt(n_y)*(sum(b_post[[1]] %*% Z_y) - sqrt(n_y)*sum(unlist(b_post[[1]] %*% R %*% a_post[[1]])))
    mu_gamma_num_part2 <- rho * sqrt(n_y) * sum(b_post[[1]] %*% Z_x)
    mu_gamma_num_part3 <- -rho*sqrt(n_x*n_y)*E_btRb

    mu_gamma_num <- sigma2_gamma_prior * (mu_gamma_num_part1 + mu_gamma_num_part2 + mu_gamma_num_part3)

    mu_gamma <- as.numeric(mu_gamma_num/mu_gamma_den)

    sigma2_gamma_curr <- as.numeric((sigma2_gamma_prior * (1-rho^2))/mu_gamma_den)

    mu2_gamma <- as.numeric(mu_gamma^2 + sigma2_gamma_curr)

    #Get KL-divergence of gamma (taken from https://stats.stackexchange.com/questions/7440/kl-divergence-between-two-univariate-gaussians):
    kl_gamma <- (1/2) - log(sqrt(sigma2_gamma_prior/sigma2_gamma_curr)) - (sigma2_gamma_curr + mu_gamma^2)/(2*sigma2_gamma_prior) #Should always be <0 (its negative kl divergence)

    #Now calculate the ELBO and check convergence
    elbo_curr <- elbo_sampleOverlap(ZxRinvZx, ZyRinvZy, ZxRinvZy, Z_x, Z_y, mu_b, mu2_b, alpha_b,
                                    mu_a, mu2_a, alpha_a, n_x, n_y, mu_gamma, mu2_gamma, R, Rinv,
                                    kl_a, kl_b, kl_gamma, rho)

    elbo_conv_vec <- c(elbo_conv_vec, elbo_curr)
    #elbo_full_vec <- c(elbo_full_vec, elbo_curr)

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
  total_var <- sampleOverlap_total_var_rss(V_x, V_y, mu_b, mu2_b, alpha_b, mu_a, mu2_a,
                                           alpha_a, R, Z_y, Z_x, num_samples,
                                           sigma2_gamma_prior, rho, n_y)

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
    to_return$cs_y <- get_cs_rss(V = V_y, alpha = alpha_a, R = list(R))
  }

  return(to_return)
}
