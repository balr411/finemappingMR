# Checks that the credible set helper functions copied from susieR (in_CS,
# get_purity, muffled_cov2cor; see R/credible_set_methods.R) give identical
# output to the susieR originals.
#
# Usage (from the package root; needs susieR 0.12.x to compare against):
#   SUSIE_LIB=<dir containing susieR 0.12.x> Rscript data-raw/check_credible_sets_vs_susieR.R

lib_loc <- "/net/mario/cluster/lib/R/site-library-jammy-40"; .libPaths(c(lib_loc, .libPaths()))
library("susieR", lib.loc = Sys.getenv("SUSIE_LIB")); library("pkgload", lib.loc = lib_loc)
pkgload::load_all(".", quiet = TRUE)
S <- asNamespace("susieR"); ok <- TRUE
chk <- function(label, a, b){ r <- isTRUE(all.equal(a, b, tolerance = 0)); cat(sprintf("%-40s identical: %s\n", label, r)); ok <<- ok && r }
set.seed(7)
for(cov in c(0.5, 0.9, 0.95, 0.99)){
  A <- matrix(rexp(10*60), 10, 60)^3; A <- A/rowSums(A); A[3,] <- 0; A[5, 1] <- 1; A[5, -1] <- 0   # include degenerate rows
  chk(paste("in_CS coverage", cov), in_CS(A, cov), S$in_CS(A, cov))
}
X <- matrix(rnorm(300*40), 300); X[,2] <- X[,1] + rnorm(300, sd = .1); Rm <- cor(X); Cm <- crossprod(scale(X, scale = FALSE))
chk("muffled_cov2cor", muffled_cov2cor(Cm), S$muffled_cov2cor(Cm))
Cz <- Cm; Cz[7, ] <- 0; Cz[, 7] <- 0
chk("muffled_cov2cor (zero-variance column)", suppressWarnings(muffled_cov2cor(Cz)), suppressWarnings(S$muffled_cov2cor(Cz)))
for(pos in list(1, c(1,2), c(1,2,5,9), 1:40)){
  for(sq in c(FALSE, TRUE)) for(rf in c(FALSE, requireNamespace("Rfast", quietly = TRUE))) {
    chk(sprintf("get_purity |cs|=%d squared=%s rfast=%s", length(pos), sq, rf), get_purity(pos, X = NULL, Xcorr = Rm, squared = sq, use_rfast = rf), S$get_purity(pos, X = NULL, Xcorr = Rm, squared = sq, use_rfast = rf))
  }
}
set.seed(3); a <- get_purity(1:40, NULL, Rm, n = 10, use_rfast = FALSE); set.seed(3); b <- S$get_purity(1:40, NULL, Rm, n = 10, use_rfast = FALSE); chk("get_purity subsampling (n=10)", a, b)
cat("\nALL IDENTICAL:", ok, "\n")
