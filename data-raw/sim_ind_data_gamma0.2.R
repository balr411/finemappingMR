# Simulates individual-level data for testing finemappingMR_ind_data().
#
# Modified from
# Scripts/18.GWASBrewer_cis_noGammaPrior/sim_combined_cis_noGammaPrior.R
# in that:
#  - only a direct causal effect gamma = 0.2 of X on Y is simulated
#  - the outcome heritability is 0.01 + gamma^2*0.05 (i.e. there is pleiotropy,
#    as in Scripts/38.sims_GWASBrewer_NealeLab/01_01.sim_randomCausal.R)
#  - N_exposure = N_outcome = 5000
#  - none of the traditional MR / LD-pruning / single-variant analyses are run
#    and no model is fit; the individual-level data is just generated and saved
#
# As in the original, we use LD block 26 of the HapMap chromosome 19 LD pattern
# that ships with GWASBrewer (863 SNPs). 5 variants are causal for X. As the
# outcome has direct (pleiotropic) heritability, num_pleio variants (not causal
# for X) are causal for Y.
#
# Usage (from the package root), all arguments optional:
#   Rscript data-raw/sim_ind_data_gamma0.2.R run=1 num_pleio=2 N_exposure=5000 \
#       N_outcome=5000 output=data-raw/sim_ind_data_gamma0.2_N5000_run1.rds

lib_loc <- "/net/mario/cluster/lib/R/site-library-jammy-40"
library("Matrix", lib.loc = lib_loc)
library("GWASBrewer", lib.loc = lib_loc)
#hapsim has to be loaded explicitly (it is not in the site library, so fall back to the default paths)
library("hapsim", lib.loc = c(lib_loc, .libPaths()))

#Read key=value arguments
args <- commandArgs(trailingOnly = TRUE)
opts <- list(run = "1", num_causal = "5", num_pleio = "2",
             N_exposure = "5000", N_outcome = "5000",
             output = "data-raw/sim_ind_data_gamma0.2_N5000_run1.rds")
for(a in args){
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if(length(kv) != 2 || !(kv[1] %in% names(opts))){
    stop("Unrecognised argument: ", a, ". Use key=value with keys: ", paste(names(opts), collapse = ", "))
  }
  opts[[kv[1]]] <- kv[2]
}

run <- as.numeric(opts$run)
num_causal <- as.numeric(opts$num_causal)
num_pleio <- as.numeric(opts$num_pleio)
N_exposure <- as.numeric(opts$N_exposure)
N_outcome <- as.numeric(opts$N_outcome)

gamma_true_curr <- 0.2
herit_x <- 0.05
herit_y <- 0.01 + gamma_true_curr^2*0.05

#Simulate GWAS for the cis-region of the 26th LD block
J <- 863
idx_start <- sum(sapply(ld_mat_list[1:25], nrow)) + 1
idx_end <- idx_start + nrow(ld_mat_list[[26]]) - 1
af_curr <- AF[idx_start:idx_end]
ld_curr <- ld_mat_list[[26]]
stopifnot(nrow(ld_curr) == J, length(af_curr) == J)

#Specify causal relationship
G <- matrix(c(0, gamma_true_curr,
              0, 0), 2, 2, byrow = TRUE)
colnames(G) <- row.names(G) <- c("X", "Y")

#Pick exactly which variants are causal
set.seed(run)
idx_X_causal <- sample(1:J, num_causal)
set.seed(run)
idx_Y_causal <- sample(c(1:J)[-idx_X_causal], num_pleio)

idx_X_causal <- idx_X_causal[order(idx_X_causal)]
idx_Y_causal <- idx_Y_causal[order(idx_Y_causal)]

pi_mat <- matrix(0, nrow = J, ncol = 2)
pi_mat[idx_X_causal, 1] <- 1
pi_mat[idx_Y_causal, 2] <- 1

#Simulate summary statistics (this also gives the true effects)
N_vec <- c(N_exposure, N_outcome)
set.seed(run)
sim_dat1_LD <- sim_mv(G = G,
                      J = J,
                      N = N_vec,
                      h2 = c(herit_x, herit_y),
                      pi = pi_mat,
                      R_LD = list(ld_curr),
                      af = af_curr,
                      est_s = TRUE)

#Now get the individual-level data; the first N_exposure rows are the exposure
#study and the last N_outcome rows are the outcome study
set.seed(run)
new_ind_dat1 <- resample_inddata(N = N_vec,
                                 dat = sim_dat1_LD,
                                 R_LD = list(ld_curr),
                                 af = af_curr,
                                 calc_sumstats = TRUE)

rows_x <- 1:N_exposure
rows_y <- (N_exposure + 1):(N_exposure + N_outcome)

x <- new_ind_dat1$Y[rows_x, 1]
y <- new_ind_dat1$Y[rows_y, 2]
G_x <- new_ind_dat1$X[rows_x, ]
G_y <- new_ind_dat1$X[rows_y, ]

#Genotypes are 0/1/2 so store them as integers to keep the file small
stopifnot(all(G_x == round(G_x)), all(G_y == round(G_y)))
G_x <- matrix(as.integer(G_x), nrow = nrow(G_x), ncol = ncol(G_x))
G_y <- matrix(as.integer(G_y), nrow = nrow(G_y), ncol = ncol(G_y))

sim_data <- list(
  # Individual-level data
  x = as.numeric(x),
  y = as.numeric(y),
  G_x = G_x,
  G_y = G_y,

  # Truth
  gamma_true = gamma_true_curr,
  idx_X_causal = idx_X_causal,
  idx_Y_causal = idx_Y_causal,
  beta_joint = new_ind_dat1$beta_joint, # true joint effects on X (col 1) and Y (col 2)
  h2 = c(x = herit_x, y = herit_y),
  sigma2_x = new_ind_dat1$Sigma_E[1, 1], # true residual variances
  sigma2_y = new_ind_dat1$Sigma_E[2, 2],

  # Single-variant summary statistics (from the individual-level data)
  beta_hat = new_ind_dat1$beta_hat,
  s_estimate = new_ind_dat1$s_estimate,

  # LD and allele frequencies used to simulate
  af = af_curr,
  ld_expected = ld_curr,

  # Settings
  settings = list(run = run, N_exposure = N_exposure, N_outcome = N_outcome,
                  num_causal = num_causal, num_pleio = num_pleio, J = J,
                  ld_block = 26,
                  gwasbrewer_version = as.character(packageVersion("GWASBrewer", lib.loc = lib_loc)))
)

saveRDS(sim_data, file = opts$output)
cat("Saved simulated data to", opts$output, "\n")
