#' @title Relevant functions for running finemappingMR with Z-score derivation
#'
#' @keywords internal

negloglik_rss_alpha <- function(sl2, RRinvR, Z_star_l, n_y){
  omega_j <- n_y * RRinvR + 1/sl2

  return(-matrixStats::logSumExp(-0.5*log(sl2*omega_j) + Z_star_l^2/(2 * omega_j)))
}

negloglik_rss_b <- function(sl2, RRinvR, Z_star_l, n_x, n_y, Egamma2){

  omega_j <- (n_x + n_y*Egamma2) * RRinvR + 1/sl2

  return(-matrixStats::logSumExp(-0.5*log(sl2*omega_j) + Z_star_l^2/(2 * omega_j)))
}
