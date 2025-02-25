#' @title Get credible sets for exposure or outcome
#'
#' @description Function that gets credible sets for the exposure or outcome.
#' Note that much of the code is directly copied from the susieR package, but
#' modified so that a SuSiE object is not needed
#'
#' @param V The L x M matrix of prior variance estimates
#'
#' @param alpha The list of length M of L x p matrices of pips for each region
#'
#' @param G_t_G The list of length M of covariance matrices for the genotypes
#'
#' @return A list of lists containing the credible sets for each region
#'
#' @export

get_cs <- function(V, alpha, G_t_G){
  M <- length(alpha)

  cs_x <- list()
  for(m in 1:M){
    include_idx <- V[,m] > 1e-9
    status <- susieR:::in_CS(alpha[[m]], coverage = 0.95)
    cs <- lapply(1:nrow(status),function(x) which(status[x,]!=0))
    claimed_coverage <- sapply(1:length(cs),
                               function (x) sum(alpha[[m]][x,][cs[[x]]]))
    include_idx <- include_idx * (lapply(cs,length) > 0)

    include_idx <- include_idx * (!duplicated(cs)) #Generally only done if dedup = TRUE in SuSiE but I am going to do it
    include_idx <- as.logical(include_idx)

    if (sum(include_idx) == 0){
      cs_x[[m]] <- NA
    }else{
      cs <- cs[include_idx]
      claimed_coverage <- claimed_coverage[include_idx]

      # Compute and filter by "purity".
      use_rfast <- requireNamespace("Rfast",quietly = TRUE)

      purity <- NULL
      for (j in 1:length(cs)) {
        purity <- rbind(purity, matrix(susieR:::get_purity(cs[[j]], X = NULL, Xcorr = susieR:::muffled_cov2cor(G_t_G[[m]]), squared = FALSE, n = 100, use_rfast), 1, 3)) #double check the muffled_cov2cor?
      }

      purity <- as.data.frame(purity)
      colnames(purity) <- c("min.abs.corr","mean.abs.corr","median.abs.corr")
      min_abs_corr <- 0.5
      squared <- FALSE
      threshold <- ifelse(squared,min_abs_corr^2,min_abs_corr)
      is_pure <- which(purity[,1] >= threshold)

      if (length(is_pure) > 0) {
        cs        <- cs[is_pure]
        purity    <- purity[is_pure,]

        cs_x[[m]] <- cs
      }else{
        cs_x[[m]] <- NA
      }
    }
  }

  return(cs_x)

}
