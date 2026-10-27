# Helper functions for computing credible sets without depending on susieR.
#
# NOTE: the functions in this file (in_CS, in_CS_x, n_in_CS_x, get_purity and
# muffled_cov2cor) were COPIED from the susieR package (version 0.12.35,
# https://github.com/stephenslab/susieR) so that credible sets can be computed
# without a SuSiE object and without calling susieR internals. They are used by
# get_cs() and get_cs_rss(). They are unchanged except that this note and the
# license text below were added.
#
# susieR is distributed under the BSD 3-clause license:
#
#   Copyright (c) 2017-2022, Gao Wang, Peter Carbonetto, Yuxin Zou,
#   Kaiqian Zhang, Matthew Stephens
#   All rights reserved.
#
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are met:
#
#   1. Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.
#   2. Redistributions in binary form must reproduce the above copyright notice,
#      this list of conditions and the following disclaimer in the documentation
#      and/or other materials provided with the distribution.
#   3. Neither the name of the copyright holder nor the names of its
#      contributors may be used to endorse or promote products derived from this
#      software without specific prior written permission.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#   POSSIBILITY OF SUCH DAMAGE.

# Copied from susieR: indicator matrix of which variants are in the credible
# set of each effect (res is the L x p matrix of PIPs)
in_CS <- function (res, coverage = 0.9)
{
  if (inherits(res, "susie"))
    res = res$alpha
  return(t(apply(res, 1, function(x) in_CS_x(x, coverage))))
}

# Copied from susieR
in_CS_x <- function (x, coverage = 0.9)
{
  n = n_in_CS_x(x, coverage)
  o = order(x, decreasing = TRUE)
  result = rep(0, length(x))
  result[o[1:n]] = 1
  return(result)
}

# Copied from susieR
n_in_CS_x <- function (x, coverage = 0.9)
  sum(cumsum(sort(x, decreasing = TRUE)) < coverage) + 1

# Copied from susieR: min, mean and median absolute correlation among the
# variants in a credible set (Xcorr is the correlation matrix; X is not needed
# when Xcorr is given)
get_purity <- function (pos, X, Xcorr, squared = FALSE, n = 100, use_rfast)
{
  if (missing(use_rfast))
    use_rfast = requireNamespace("Rfast", quietly = TRUE)
  if (use_rfast) {
    get_upper_tri = Rfast::upper_tri
    get_median = Rfast::med
  }
  else {
    get_upper_tri = function(R) R[upper.tri(R)]
    get_median = stats::median
  }
  if (length(pos) == 1)
    return(c(1, 1, 1))
  else {
    if (length(pos) > n)
      pos = sample(pos, n)
    if (is.null(Xcorr)) {
      X_sub = X[, pos]
      X_sub = as.matrix(X_sub)
      value = abs(get_upper_tri(muffled_corr(X_sub)))
    }
    else value = abs(get_upper_tri(Xcorr[pos, pos]))
    if (squared)
      value = value^2
    return(c(min(value), sum(value)/length(value), get_median(value)))
  }
}

# Copied from susieR
muffled_cov2cor <- function (x)
  withCallingHandlers(cov2cor(x), warning = function(w) {
    if (grepl("had 0 or NA entries; non-finite result is doubtful",
              w$message))
      invokeRestart("muffleWarning")
  })

# Copied from susieR (only needed when get_purity is called with Xcorr = NULL,
# which this package does not do)
muffled_corr <- function (x)
  withCallingHandlers(cor(x), warning = function(w) {
    if (grepl("the standard deviation is zero", w$message))
      invokeRestart("muffleWarning")
  })
