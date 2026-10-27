# Compares finemappingMR_sampleOverlap() (Z-scores + correlation matrix, rho = 0)
# with run_bayes_method_ss() (sufficient statistics) and finemappingMR_ind_data()
# (individual-level data) on the data saved by data-raw/sim_ind_data_gamma0.2.R.
#
# Z-scores are the marginal-regression (with intercept) t-statistics from the
# individual-level data, R is the correlation of the exposure genotypes.
#
# Usage (from the package root):
#   Rscript data-raw/compare_sampleOverlap.R [data.rds]
# (the package does not need susieR)

lib_loc <- "/net/mario/cluster/lib/R/site-library-jammy-40"
.libPaths(c(lib_loc, .libPaths()))
library("stringr", lib.loc = lib_loc)
library("pkgload", lib.loc = lib_loc)
pkgload::load_all(".", quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
data_file <- if(length(args) > 0) args[1] else "data-raw/sim_ind_data_gamma0.2_N5000_run1.rds"
sim <- readRDS(data_file)

x <- sim$x; y <- sim$y
G_x <- matrix(as.numeric(sim$G_x), nrow(sim$G_x), ncol(sim$G_x))
G_y <- matrix(as.numeric(sim$G_y), nrow(sim$G_y), ncol(sim$G_y))
n_x <- length(x); n_y <- length(y); p <- ncol(G_x)

#Marginal regression Z-scores (identical to the OLS t-statistic with an intercept)
marg_z <- function(G, y){
  n <- length(y)
  r <- as.numeric(cor(G, y))
  r * sqrt((n - 2)/(1 - r^2))
}
Z_x <- marg_z(G_x, x)
Z_y <- marg_z(G_y, y)
R <- cor(G_x) #correlation from the exposure data set
cat("max |Z_x| =", round(max(abs(Z_x)), 2), " max |Z_y| =", round(max(abs(Z_y)), 2), "\n")

#Sufficient statistics
xc <- x - mean(x); yc <- y - mean(y)
Gxc <- sweep(G_x, 2, colMeans(G_x)); Gyc <- sweep(G_y, 2, colMeans(G_y))

num_samples <- 200
set.seed(1)
fit_ss <- run_bayes_method_ss(list(crossprod(Gxc)), list(crossprod(Gxc, xc)), sum(xc^2),
                              list(crossprod(Gyc)), list(crossprod(Gyc, yc)), sum(yc^2),
                              n_x, n_y, tol = 1e-4, calc_cs_x = TRUE, calc_cs_y = TRUE,
                              num_samples = num_samples)
set.seed(1)
fit_ind <- finemappingMR_ind_data(x, y, G_x, G_y, tol = 1e-4, calc_cs_x = TRUE,
                                  calc_cs_y = TRUE, num_samples = num_samples)
set.seed(1)
fit_z <- finemappingMR_sampleOverlap(Z_x, Z_y, R, rho = 0, n_x = n_x, n_y = n_y,
                                     tol = 1e-4, calc_cs_x = TRUE, calc_cs_y = TRUE,
                                     num_samples = num_samples)

pip <- function(alpha, V){
  a <- alpha[[1]][V[, 1] > 1e-9, , drop = FALSE]
  1 - apply(1 - a, 2, prod)
}
post_mean <- function(fit, which = c("b", "a")){
  which <- match.arg(which)
  if(which == "b") colSums(fit$alpha_b[[1]] * fit$mu_b[[1]]) else colSums(fit$alpha_a[[1]] * fit$mu_a[[1]])
}

fits <- list(SS = fit_ss, IND = fit_ind, Zscore = fit_z)
tab <- do.call(rbind, lapply(names(fits), function(nm){
  r <- fits[[nm]]$res
  data.frame(method = nm, gamma = r$gamma, gamma_var = r$gamma_var,
             gamma_total_var = r$gamma_total_var, iter = r$iter,
             n_eff_x = sum(fits[[nm]]$V_x[, 1] > 1e-9), n_eff_y = sum(fits[[nm]]$V_y[, 1] > 1e-9))
}))
cat("\n== Gamma (true = 0.2) ==\n"); print(tab, row.names = FALSE, digits = 4)

cat("\nX-causal:", sim$idx_X_causal, " Y-causal:", sim$idx_Y_causal, "\n")

cat("\n== Exposure PIPs: top variants ==\n")
pips <- sapply(fits, function(f) pip(f$alpha_b, f$V_x))
top <- order(-apply(pips, 1, max))[1:12]
print(round(cbind(snp = top, pips[top, ]), 3))
cat("cor of exposure PIPs: SS-IND", round(cor(pips[, "SS"], pips[, "IND"]), 4),
    " SS-Zscore", round(cor(pips[, "SS"], pips[, "Zscore"]), 4), "\n")

cat("\n== Outcome (pleiotropy) PIPs: top variants ==\n")
pips_y <- sapply(fits, function(f) if(any(f$V_y[, 1] > 1e-9)) pip(f$alpha_a, f$V_y) else rep(0, p))
top_y <- order(-apply(pips_y, 1, max))[1:8]
print(round(cbind(snp = top_y, pips_y[top_y, ]), 3))

cat("\n== Credible sets ==\n")
for(nm in names(fits)){
  cat(nm, "exposure:\n"); cs <- fits[[nm]]$cs_x[[1]]
  if(is.list(cs)) for(s in cs) cat("   size", length(s), ": ", head(s, 8), if(length(s) > 8) "...", "\n") else cat("   none\n")
  cat(nm, "outcome:\n"); cs <- fits[[nm]]$cs_y[[1]]
  if(is.list(cs)) for(s in cs) cat("   size", length(s), ": ", head(s, 8), if(length(s) > 8) "...", "\n") else cat("   none\n")
}

cat("\n== Posterior mean effects on exposure (per-allele vs standardised scale differ) ==\n")
cat("cor of posterior mean b: SS-IND", round(cor(post_mean(fit_ss), post_mean(fit_ind)), 4),
    " SS-Zscore", round(cor(post_mean(fit_ss), post_mean(fit_z)), 4), "\n")
