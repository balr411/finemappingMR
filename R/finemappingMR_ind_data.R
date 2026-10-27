#' @title Bayesian MR estimate using individual-level data for exposure and outcome.
#'
#' @description Function that performs our Bayesian MR estimation method using
#' individual-level data for exposure and outcome. This is the individual-level
#' analogue of \code{\link{run_bayes_method_ss}}: it has the same options and
#' returns the same components, but the updates are computed from the
#' phenotypes and genotypes directly instead of from sufficient statistics.
#' A normal prior with variance \code{sigma2_gamma_prior} is placed on gamma.
#' The single effect regressions are computed with the package's own
#' implementation, so fitting does not depend on susieR (only the credible set
#' calculation, \code{calc_cs_x}/\code{calc_cs_y}, still does).
#'
#' Within each locus all quantities are computed exactly as in
#' \code{run_bayes_method_ss}; in particular, different loci are treated as
#' independent (\eqn{G_j'G_k = 0}).
#'
#' @param x Exposure individual-level phenotype data (vector). Will be centered.
#'
#' @param y Outcome individual-level phenotype data (vector). Will be centered.
#'
#' @param G_x Exposure genotype data. Either a single matrix (one locus) or a
#'   list of matrices (one per locus), each with \code{length(x)} rows.
#'   Columns are centered internally.
#'
#' @param G_y Outcome genotype data. Either a single matrix (one locus) or a
#'   list of matrices (one per locus), each with \code{length(y)} rows. Must
#'   have the same number of loci and variants per locus as \code{G_x}.
#'   Columns are centered internally.
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
#'   value supplied to \code{scaled_prior_variance_x}. Currently has
#'   no functionality.
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
#' this is used as an initial value. Default (\code{NULL}) is \code{var(x)}.
#'
#' @param residual_variance_y The variance of the outcome residual. If \code{estimate_residual_variance_y = TRUE},
#' this is used as an initial value. Default (\code{NULL}) is \code{var(y)}.
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
#'   less than \code{tol}. Note that this is used in all cases, whereas
#'   \code{run_bayes_method_ss} only uses it when \code{include_pleiotropy = FALSE}
#'   (and otherwise uses 1e-4).
#'
#' @param max_iter Maximum number of iterations before stopping estimation procedure.
#'
#' @param calc_cs_x Calculate credible sets for the exposure? Default = FALSE.
#'
#' @param calc_cs_y Calculate credible sets for the outcome? Default = FALSE.
#'
#' @param verbose Output progress? Default = FALSE.
#'
#' @param susie_init A list of SuSiE objects (one per locus) to initialize the
#'  exposure to. Default = NULL.
#'
#' @param susie_init_y A list of SuSiE objects (one per locus) to initialize the
#'  outcome to. Default = NULL.
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
#' @return A list containing various results from the estimation procedure,
#' including a data frame containing the gamma estimation results, as well as
#' the posterior first and second moments for b and alpha, the prior non-zero
#' effect variance estimates, and the PIPs for each non-zero effect, the
#' vector of ELBO iterations, the vector of likelihood components, and the KL
#' divergence terms for b and alpha. Note that the ELBO includes the
#' \eqn{y'y/(2\sigma_y^2) + x'x/(2\sigma_x^2)} constants, which are left out of
#' the ELBO returned by \code{run_bayes_method_ss}.
#'
#' @importFrom stats var
#' @importFrom stringr str_glue
#'
#' @export

finemappingMR_ind_data <- function(x, y, G_x, G_y,
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

  #A single matrix is a single locus
  if(!is.list(G_x) || is.data.frame(G_x)) G_x <- list(G_x)
  if(!is.list(G_y) || is.data.frame(G_y)) G_y <- list(G_y)

  M <- length(G_x) #Number of loci

  x <- as.numeric(x)
  y <- as.numeric(y)

  if(length(G_y) != M){
    stop("G_x and G_y must contain the same number of loci")
  }

  for(m in 1:M){
    G_x[[m]] <- matrix(as.numeric(G_x[[m]]), nrow = nrow(G_x[[m]]), ncol = ncol(G_x[[m]]))
    G_y[[m]] <- matrix(as.numeric(G_y[[m]]), nrow = nrow(G_y[[m]]), ncol = ncol(G_y[[m]]))

    if(nrow(G_x[[m]]) != length(x)){
      stop("The number of rows of each element of G_x must equal length(x)")
    }

    if(nrow(G_y[[m]]) != length(y)){
      stop("The number of rows of each element of G_y must equal length(y)")
    }

    if(ncol(G_x[[m]]) != ncol(G_y[[m]])){
      stop("G_x and G_y must have the same number of variants in each locus")
    }
  }

  n_x <- length(x)
  n_y <- length(y)

  varX <- var(x)
  varY <- var(y)

  if(is.null(residual_variance_x)){
    residual_variance_x <- varX
  }

  if(is.null(residual_variance_y)){
    residual_variance_y <- varY
  }

  sigma2_x <- residual_variance_x
  sigma2_y <- residual_variance_y

  #Center x and y, and the columns of G_x and G_y (intercept = TRUE in SuSiE).
  #Note that columns are not scaled (standardize = FALSE in SuSiE)
  x <- x - mean(x)
  y <- y - mean(y)

  #The "d" attribute is the diagonal of G'G, i.e. the column sums of squares of
  #the centered columns
  for(m in 1:M){
    G_x[[m]] <- sweep(G_x[[m]], 2, colMeans(G_x[[m]]))
    G_y[[m]] <- sweep(G_y[[m]], 2, colMeans(G_y[[m]]))

    attr(G_x[[m]], "d") <- colSums(G_x[[m]]^2)
    attr(G_y[[m]], "d") <- colSums(G_y[[m]]^2)
  }

  #Initialize all of the estimates
  est_init <- fmr_init(M, L_x, L_y, susie_init, susie_init_y,
                       scaled_prior_variance_x, scaled_prior_variance_y, varX,
                       varY, G_x, G_y, mu_gamma_init, mu2_gamma_init)

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

  #If pleiotropy is not estimated, set all prior variance estimates for the
  #outcome to 0 so that the total variance calculation will work. The
  #posterior for alpha is then never updated.
  if(!include_pleiotropy){
    V_y <- matrix(0, nrow = L_y, ncol = M)
  }

  #Order in which the parameters are updated
  if(!include_pleiotropy){
    update_order <- c("b", "gamma")
  }else if(beta_gamma_alpha){
    update_order <- c("b", "gamma", "alpha")
  }else{
    update_order <- c("alpha", "b", "gamma")
  }

  #Fitted values for a and b (without gamma) at the current posterior means
  fit_a <- function(m){
    as.numeric(G_y[[m]] %*% colSums(alpha_a[[m]] * mu_a[[m]]))
  }
  fit_b_x <- function(m){
    as.numeric(G_x[[m]] %*% colSums(alpha_b[[m]] * mu_b[[m]]))
  }
  fit_b_y <- function(m){
    as.numeric(G_y[[m]] %*% colSums(alpha_b[[m]] * mu_b[[m]]))
  }

  elbo_conv_vec <- c()
  lik_x <- c()
  lik_y <- c()

  conv <- FALSE
  iter <- 1

  while(!conv & iter < max_iter){
    for(step in update_order){

      if(step == "alpha"){
        #Update alpha
        for(m in 1:M){
          Xr_a <- fit_a(m)
          Xr_b_y <- mu_gamma * fit_b_y(m)

          for(l in 1:L_y){
            # Remove lth effect from fitted values.
            # Note that we shouldn't actually have to include each of the m loci since we are assuming that Gj^TGk = 0
            Xr_a <- Xr_a - as.numeric(G_y[[m]] %*% (alpha_a[[m]][l,] * mu_a[[m]][l,]))

            R <- y - Xr_a - Xr_b_y
            test_curr <- fmr_ser(as.numeric(crossprod(G_y[[m]], R)), attr(G_y[[m]], "d"), V_y[l, m], residual_variance = sigma2_y)

            mu_a[[m]][l,] <- test_curr$mu
            mu2_a[[m]][l,] <- test_curr$mu2
            alpha_a[[m]][l,] <- test_curr$alpha
            V_y[l, m] <- test_curr$V
            kl_a[[m]][l] <- test_curr$kl

            Xr_a <- Xr_a + as.numeric(G_y[[m]] %*% (alpha_a[[m]][l,] * mu_a[[m]][l,]))
          }
        }

      }else if(step == "b"){
        #Update b
        for(m in 1:M){
          Xr_a <- fit_a(m)
          Xr_b_x <- fit_b_x(m)
          Xr_b_y <- fit_b_y(m)

          #The b update is a regression on G_star = (G_x/sqrt(sigma2_x) \\ sqrt(mu2_gamma/sigma2_y) * G_y)
          #with residual variance 1. The single effect regression only needs
          #the diagonal of G_star'G_star (the same for every effect l) and
          #G_star'R, so G_star does not need to be built.
          d_star <- attr(G_x[[m]], "d")/sigma2_x + (mu2_gamma/sigma2_y) * attr(G_y[[m]], "d")

          resid_var_curr <- 1

          for(l in 1:L_x){
            # Remove lth effect from fitted values.
            b_l_old <- alpha_b[[m]][l,] * mu_b[[m]][l,]
            Xr_b_x <- Xr_b_x - as.numeric(G_x[[m]] %*% b_l_old)
            Xr_b_y <- Xr_b_y - as.numeric(G_y[[m]] %*% b_l_old)

            #G_star'R = G_x'(x - Xr_b_x)/sigma2_x + G_y'(mu_gamma*(y - Xr_a) - mu2_gamma*Xr_b_y)/sigma2_y
            Xty_star <- as.numeric(crossprod(G_x[[m]], x - Xr_b_x))/sigma2_x +
              as.numeric(crossprod(G_y[[m]], mu_gamma * (y - Xr_a) - mu2_gamma * Xr_b_y))/sigma2_y

            test_curr <- fmr_ser(Xty_star, d_star, V_x[l, m], residual_variance = resid_var_curr)
            #Note the residual variance comes from what is multiplied by the Eq term in the ELBO

            mu_b[[m]][l,] <- test_curr$mu
            mu2_b[[m]][l,] <- test_curr$mu2
            alpha_b[[m]][l,] <- test_curr$alpha
            V_x[l, m] <- test_curr$V
            kl_b[[m]][l] <- test_curr$kl

            b_l_new <- alpha_b[[m]][l,] * mu_b[[m]][l,]
            Xr_b_x <- Xr_b_x + as.numeric(G_x[[m]] %*% b_l_new)
            Xr_b_y <- Xr_b_y + as.numeric(G_y[[m]] %*% b_l_new)
          }
        }

      }else{
        #Now update gamma
        mu_gamma_num <- 0
        sum_transpose <- 0
        for(m in 1:M){
          Xr_a <- fit_a(m)
          Xr_b_y <- fit_b_y(m)

          # sum_m b_m' G_ym' (y - G_ym a_m)
          mu_gamma_num <- mu_gamma_num + sum(Xr_b_y * (y - Xr_a))

          # E[sum_m (G_ym b_m)' (G_ym b_m)]
          sum_transpose <- sum_transpose + test_ER(G_y[[m]], y, alpha_b[[m]], mu_b[[m]], mu2_b[[m]])
        }

        mu_gamma_num <- sigma2_gamma_prior * mu_gamma_num
        mu_gamma_den <- sigma2_y + sigma2_gamma_prior * sum_transpose

        mu_gamma <- as.numeric(mu_gamma_num/mu_gamma_den)

        sigma2_gamma_curr <- as.numeric((sigma2_gamma_prior*sigma2_y)/mu_gamma_den)

        mu2_gamma <- as.numeric(mu_gamma^2 + sigma2_gamma_curr)

        #Get KL-divergence of gamma (taken from https://stats.stackexchange.com/questions/7440/kl-divergence-between-two-univariate-gaussians):
        kl_gamma <- (1/2) - log(sqrt(sigma2_gamma_prior/sigma2_gamma_curr)) - (sigma2_gamma_curr + mu_gamma^2)/(2*sigma2_gamma_prior) #Should always be <0 (its negative kl divergence)
      }
    }

    Xa_list <- lapply(1:M, fit_a)
    Xb_y_list <- lapply(1:M, fit_b_y)

    elbo_full <- bayes_elbo_ind(sigma2_y, sigma2_x, x, y, G_x, G_y, Xa_list, Xb_y_list,
                                alpha_b, mu_b, mu2_b, alpha_a, mu_a, mu2_a,
                                mu_gamma, mu2_gamma, kl_a, kl_b, kl_gamma)

    elbo_curr <- elbo_full[[1]]
    lik_y_curr <- elbo_full[[2]]
    lik_x_curr <- elbo_full[[3]]

    elbo_conv_vec <- c(elbo_conv_vec, elbo_curr)
    lik_x <- c(lik_x, lik_x_curr)
    lik_y <- c(lik_y, lik_y_curr)

    if(iter > 1){
      conv <- abs(elbo_conv_vec[iter] - elbo_conv_vec[iter-1]) < tol
    }else{
      conv <- FALSE
    }

    iter <- iter + 1

    if(verbose){
      print(str_glue("Starting iteration {iter}, gamma_curr = {mu_gamma}"))
    }
  }

  if(!conv){
    warning("Maximum iteration number reached, convergence not attained")
  }

  #The total variance calculation and credible sets use G'G and G'y
  Gy_t_Gy <- lapply(G_y, function(G) crossprod(G))
  Gy_t_y <- lapply(G_y, function(G) crossprod(G, y))

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
    cs_x_tmp <- get_cs(V = V_x, alpha = alpha_b, G_t_G = lapply(G_x, function(G) crossprod(G)))
    to_return$cs_x <- cs_x_tmp$cs
    to_return$cs_x_idx <- cs_x_tmp$idx
  }

  if(calc_cs_y){
    cs_y_tmp <- get_cs(V = V_y, alpha = alpha_a, G_t_G = Gy_t_Gy)
    to_return$cs_y <- cs_y_tmp$cs
    to_return$cs_y_idx <- cs_y_tmp$idx
  }

  return(to_return)
}
