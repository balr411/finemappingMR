# Checks that finemappingMR_ind_data() (individual-level data) gives the same
# answer as run_bayes_method_ss() (sufficient statistics computed from the same
# individual-level data), for each option combination, using the data saved by
# data-raw/sim_ind_data_gamma0.2.R.
#
# Usage (from the package root):
#   SUSIE_LIB=<dir containing susieR 0.12.x> Rscript data-raw/check_ind_data_vs_ss.R [data.rds]
#
# NOTE: the package no longer depends on susieR (single effect regressions and
# credible sets are implemented in the package). susieR is only used here to build
# the SuSiE fits used to test the susie_init option; any version can be used.

lib_loc <- "/net/mario/cluster/lib/R/site-library-jammy-40"
.libPaths(c(lib_loc, .libPaths()))

susie_lib <- Sys.getenv("SUSIE_LIB", unset = "")
if(nzchar(susie_lib)){
  library("susieR", lib.loc = susie_lib)
}else{
  library("susieR", lib.loc = lib_loc)
}
cat("susieR version:", as.character(packageVersion("susieR")), "\n")
library("stringr", lib.loc = lib_loc)
library("pkgload", lib.loc = lib_loc)
pkgload::load_all(".", quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
data_file <- if(length(args) > 0) args[1] else "data-raw/sim_ind_data_gamma0.2_N5000_run1.rds"
sim <- readRDS(data_file)

x <- sim$x; y <- sim$y
G_x <- matrix(as.numeric(sim$G_x), nrow(sim$G_x), ncol(sim$G_x))
G_y <- matrix(as.numeric(sim$G_y), nrow(sim$G_y), ncol(sim$G_y))
cat("N_x =", length(x), " N_y =", length(y), " p =", ncol(G_x), "\n")
cat("true gamma =", sim$gamma_true, "; X causal:", sim$idx_X_causal, "; Y causal:", sim$idx_Y_causal, "\n")

#Sufficient statistics from the same (centered) individual-level data
get_ss <- function(x, y, G_x_list, G_y_list){
  xc <- x - mean(x); yc <- y - mean(y)
  Gxc <- lapply(G_x_list, function(G) sweep(G, 2, colMeans(G)))
  Gyc <- lapply(G_y_list, function(G) sweep(G, 2, colMeans(G)))
  list(Gx_t_Gx = lapply(Gxc, crossprod),
       Gx_t_x = lapply(Gxc, function(G) crossprod(G, xc)),
       xtx = sum(xc^2),
       Gy_t_Gy = lapply(Gyc, crossprod),
       Gy_t_y = lapply(Gyc, function(G) crossprod(G, yc)),
       yty = sum(yc^2))
}

max_rel_diff <- function(a, b){
  a <- unlist(a); b <- unlist(b)
  max(abs(a - b)/pmax(1, abs(a), abs(b)))
}

compare <- function(label, G_x_in, G_y_in, ...){
  ss <- get_ss(x, y, if(is.list(G_x_in)) G_x_in else list(G_x_in),
               if(is.list(G_y_in)) G_y_in else list(G_y_in))
  M <- length(ss$Gx_t_Gx)

  set.seed(1)
  t_ss <- system.time(
    fit_ss <- run_bayes_method_ss(ss$Gx_t_Gx, ss$Gx_t_x, ss$xtx, ss$Gy_t_Gy, ss$Gy_t_y, ss$yty,
                                  length(x), length(y), ...))["elapsed"]
  set.seed(1)
  t_ind <- system.time(
    fit_ind <- finemappingMR_ind_data(x, y, G_x_in, G_y_in, ...))["elapsed"]

  #ELBO differs by a constant: y'y/(2 sigma2_y) + x'x/(2 sigma2_x)
  const <- ss$yty/(2*var(y)) + ss$xtx/(2*var(x))
  d <- c(gamma = max_rel_diff(fit_ss$res$gamma, fit_ind$res$gamma),
         gamma_var = max_rel_diff(fit_ss$res$gamma_var, fit_ind$res$gamma_var),
         iter = fit_ss$res$iter - fit_ind$res$iter,
         elbo = max_rel_diff(fit_ss$elbo, fit_ind$elbo + const),
         alpha_b = max_rel_diff(fit_ss$alpha_b, fit_ind$alpha_b),
         mu_b = max_rel_diff(fit_ss$mu_b, fit_ind$mu_b),
         alpha_a = max_rel_diff(fit_ss$alpha_a, fit_ind$alpha_a),
         mu_a = max_rel_diff(fit_ss$mu_a, fit_ind$mu_a),
         V_x = max_rel_diff(fit_ss$V_x, fit_ind$V_x),
         V_y = max_rel_diff(fit_ss$V_y, fit_ind$V_y),
         kl_b = max_rel_diff(fit_ss$kl_b, fit_ind$kl_b))
  cat("\n=====", label, "=====\n")
  cat(sprintf("SS:  gamma = %.6f  var = %.3e  total var = %.3e  iter = %d  (%.1fs)\n",
              fit_ss$res$gamma, fit_ss$res$gamma_var, fit_ss$res$gamma_total_var, fit_ss$res$iter, t_ss))
  cat(sprintf("IND: gamma = %.6f  var = %.3e  total var = %.3e  iter = %d  (%.1fs)\n",
              fit_ind$res$gamma, fit_ind$res$gamma_var, fit_ind$res$gamma_total_var, fit_ind$res$iter, t_ind))
  cat("max relative differences (SS vs IND):\n")
  print(signif(d, 3))
  cat("PASS (all < 1e-6):", all(abs(d) < 1e-6), "\n")
  invisible(list(ss = fit_ss, ind = fit_ind, d = d))
}

#Defaults of the SS function with the usual simulation settings
res1 <- compare("alpha-beta-gamma (default), 1 locus", G_x, G_y,
                tol = 1e-4, calc_cs_x = TRUE, calc_cs_y = TRUE, num_samples = 50)

res2 <- compare("beta-gamma-alpha, 1 locus", G_x, G_y,
                tol = 1e-4, beta_gamma_alpha = TRUE, num_samples = 50)

res3 <- compare("no pleiotropy, 1 locus", G_x, G_y,
                tol = 1e-4, include_pleiotropy = FALSE, num_samples = 50)

res4 <- compare("L_x = 5, L_y = 3, sigma2_gamma_prior = 1", G_x, G_y,
                L_x = 5, L_y = 3, tol = 1e-4, sigma2_gamma_prior = 1, num_samples = 50)

#Two loci: split the variants in half
idx1 <- 1:430; idx2 <- 431:ncol(G_x)
res5 <- compare("2 loci", list(G_x[, idx1], G_x[, idx2]), list(G_y[, idx1], G_y[, idx2]),
                tol = 1e-4, calc_cs_x = TRUE, num_samples = 50)

#Initialising with SuSiE fits
susie_x <- susie_suff_stat(XtX = crossprod(sweep(G_x, 2, colMeans(G_x))),
                           Xty = crossprod(sweep(G_x, 2, colMeans(G_x)), x - mean(x)),
                           yty = sum((x - mean(x))^2), n = length(x),
                           L = 10, estimate_residual_variance = FALSE)
susie_y <- susie_suff_stat(XtX = crossprod(sweep(G_y, 2, colMeans(G_y))),
                           Xty = crossprod(sweep(G_y, 2, colMeans(G_y)), y - mean(y)),
                           yty = sum((y - mean(y))^2), n = length(y),
                           L = 10, estimate_residual_variance = FALSE)
res6 <- compare("SuSiE initialization, 1 locus", G_x, G_y,
                tol = 1e-4, susie_init = list(susie_x), susie_init_y = list(susie_y),
                num_samples = 50)

#Credible sets agree?
cat("\nCredible sets (default fit), exposure (IND):\n"); print(res1$ind$cs_x)
cat("Credible sets (default fit), exposure (SS):\n"); print(res1$ss$cs_x)
cat("Credible sets (default fit), outcome (IND):\n"); print(res1$ind$cs_y)
cat("identical cs_x:", identical(res1$ind$cs_x, res1$ss$cs_x), " cs_y:", identical(res1$ind$cs_y, res1$ss$cs_y), "\n")

#Sanity: ELBO should be non-decreasing
cat("\nELBO non-decreasing (default fit, up to 1e-6):", all(diff(res1$ind$elbo) > -1e-6), "\n")
