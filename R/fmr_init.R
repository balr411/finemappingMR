#' @title Initialize the posteriors for b and alpha
#'
#' @keywords internal

fmr_init <- function(M, L_x, L_y, susie_init, susie_init_y,
                     scaled_prior_variance_x, scaled_prior_variance_y, varX,
                     varY, Gx_t_Gx, Gy_t_Gy, mu_gamma_init, mu2_gamma_init){

  #Check if a SuSiE object was given for initialization for the exposure
  if(!is.null(susie_init)){

    V_x <- matrix(rep(0, L_x*M), nrow = L_x, ncol = M)
    mu_b <- list()
    mu2_b <- list()
    alpha_b <- list()

    for(i in 1:M){
      V_x[,i] <- susie_init[[i]]$V
      mu_b[[i]] <- susie_init[[i]]$mu
      mu2_b[[i]] <- susie_init[[i]]$mu2
      alpha_b[[i]] <- susie_init[[i]]$alpha
    }

    kl_b <- lapply(Gx_t_Gx, FUN = function(x) rep(0, L_x)) #Note not sure that this would be correct
  }else{
    V_x_init <- scaled_prior_variance_x * varX
    V_x <- matrix(rep(V_x_init, L_x*M), nrow = L_x, ncol = M)

    #b estimates
    mu_b <- lapply(Gx_t_Gx, function(x) matrix(0, nrow = L_x, ncol = ncol(x)))
    mu2_b <- lapply(Gx_t_Gx, function(x) matrix(0, nrow = L_x, ncol = ncol(x)))
    alpha_b <- lapply(Gx_t_Gx, function(x) matrix(1/ncol(x), nrow = L_x, ncol = ncol(x)))
    kl_b <- lapply(Gx_t_Gx, FUN = function(x) rep(0, L_x))
  }

  #Check if a SuSiE object was given for initialization for the outcome
  if(!is.null(susie_init_y)){
    V_y <- matrix(rep(0, L_y*M), nrow = L_y, ncol = M)
    mu_a <- list()
    mu2_a <- list()
    alpha_a <- list()

    for(i in 1:M){
      V_y[,i] <- susie_init_y[[i]]$V
      mu_a[[i]] <- susie_init_y[[i]]$mu
      mu2_a[[i]] <- susie_init_y[[i]]$mu2
      alpha_a[[i]] <- susie_init_y[[i]]$alpha
    }

    kl_a <- lapply(Gy_t_Gy, FUN = function(x) rep(0, L_y)) #Note not sure that this would be correct

  }else{
    #Initialize variance estimates
    V_y_init <- scaled_prior_variance_y * varY
    V_y <- matrix(rep(V_y_init, L_y*M), nrow = L_y, ncol = M)

    #alpha (a) estimates
    mu_a <- lapply(Gy_t_Gy, function(x) matrix(0, nrow = L_y, ncol = ncol(x)))
    mu2_a <- lapply(Gy_t_Gy, function(x) matrix(0, nrow = L_y, ncol = ncol(x)))
    alpha_a <- lapply(Gy_t_Gy, function(x) matrix(1/ncol(x), nrow = L_y, ncol = ncol(x)))
    kl_a <- lapply(Gy_t_Gy, FUN = function(x) rep(0, L_y))
  }

  #gamma estimate
  mu_gamma <- mu_gamma_init
  mu2_gamma <- mu2_gamma_init
  sigma2_gamma_curr <- mu2_gamma - mu_gamma^2
  kl_gamma <- 0

  return(list(V_x = V_x, mu_b = mu_b, mu2_b = mu2_b, alpha_b = alpha_b, kl_b = kl_b,
              V_y = V_y, mu_a = mu_a, mu2_a = mu2_a, alpha_a = alpha_a, kl_a = kl_a,
              mu_gamma = mu_gamma, mu2_gamma = mu2_gamma, sigma2_gamma_curr = sigma2_gamma_curr,
              kl_gamma = kl_gamma))
}
