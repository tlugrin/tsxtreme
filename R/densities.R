## Copyright (C) 2017-2025 Thomas Lugrin
## Implementation of current density/distribution/etc. functions
## (Scope:) Any
## List of functions: - dlapl
##                    - plapl
##                    - qlapl
##                    - rlapl
##################################################

##################################################
## LAPLACE DISTRIBUTION
#' The Laplace distribution
#' 
#' Density, distribution function, quantile function and random generation for
#' the Laplace distribution.
#' 
#' @details
#' If `loc` or `scale` are not specified, they assume the default
#' values of 0 and 1 respectively.
#' 
#' The Laplace distribution has density
#' \deqn{f(x) = \exp(- |x-\mu|/\sigma)/(2\sigma)}
#' where \eqn{\mu} is the location parameter and \eqn{\sigma} is the scale
#' parameter.
#' 
#' @section Warning: 
#' Some checks are done previous to standard evaluation, but vector computations
#' have not yet been tested thoroughly! Typically vectors not having lengths
#' multiple of each other return an error.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of samples.
#' @param loc vector of location parameters.
#' @param scale vector of scale parameters; must be non-negative.
#' @param log,log_p logical; if TRUE, probabilities `p` are given as `log(p)`;
#'   defaults to FALSE.
#' @param lower_tail logical; if TRUE (default), probabilities are
#'   \eqn{P(X\le x)}, otherwise \eqn{P(X>x)}.
#' @param log.p,lower.tail `r lifecycle::badge('deprecated')` use `log_p` and
#'   `lower_tail` instead.
#' @returns `dlapl` gives the density, `plapl` gives the distribution function,
#'   `qlapl` gives the quantile function, and `rlapl` generates random deviates.
#'   
#'   The length of the result is determined by `n` in `rlapl`, and is the
#'   maximum of the lengths of the numerical arguments for the other functions.
#'   Standard `R` vector operations are to be assumed.
#'   
#'   If `sd==0`, the limit as `sd` decreases to 0 is returned, i.e., a point
#'   mass at `loc`. The case `sd<0` is an error and nothing is returned.
#' @seealso [stats::dexp()] for the exponential distribution which is the
#'   positive part of the Laplace distribution.
#' @examples
#' ## evaluate the density function on a grid of values
#' x  <- seq(from=-5, to=5, by=0.1)
#' fx <- dlapl(x, loc=1, scale=.5)
#'
#' ## generate random samples of a mixture of Laplace distributions
#' rnd <- rlapl(1000, loc=c(-5,-3,2), scale=0.5)
#' 
#' ## an alternative:
#' rnd <- runif (1000)
#' rnd <- qlapl(rnd, loc=c(-5,-3,2), scale=0.5)
#' 
#' ## integrate the Laplace density on [a,b]
#' a <- -1
#' b <- 7
#' integral <- plapl(b)-plapl(a)
#' 
#' @export
dlapl <- function(x, loc = 0, scale = 1, log = FALSE) {
  if (any(scale < 0)) stop("scale parameter must be non-negative")
  if (log) {
    ret <- -abs(x-loc)/scale-log(2)-log(scale)
    ret[is.nan(ret)] <- 1
    return(ret)
  } else {
    ret <- exp(-abs(x-loc)/scale)/2/scale
    ret[is.nan(ret)] <- 1
    return(ret)
  }
}

#' @rdname dlapl
#' @export
plapl <- function(q, loc = 0, scale = 1, lower_tail = TRUE, log_p = FALSE,
                  lower.tail = deprecated(), log.p = deprecated()) {
  if (lifecycle::is_present(lower.tail)) {
    lifecycle::deprecate_warn("0.4.0", "plapl(lower.tail)", "plapl(lower_tail")
    lower_tail <- lower.tail
  }
  if (lifecycle::is_present(log.p)) {
    lifecycle::deprecate_warn("0.4.0", "plapl(log.p)", "plapl(log_p)")
    log_p <- log.p
  }
  if (any(scale < 0)) stop("scale parameter must be non-negative")
  ret <- sign(q-loc)*(1-exp(-abs(q-loc)/scale))/2
  ret[is.nan(ret)] <- 0.5
  if (lower_tail) {
    if (log_p)
      return(log(0.5+ret))
    else
      return(0.5+ret)
  } else {
    if (log_p)
      return(log(0.5-ret))
    else
      return(0.5-ret)
  }
}

#' @rdname dlapl
#' @export
qlapl <- function(p, loc = 0, scale = 1, lower_tail = TRUE, log_p = FALSE,
                  lower.tail = deprecated(), log.p = deprecated()) {
  if (lifecycle::is_present(lower.tail)) {
    lifecycle::deprecate_warn("0.4.0", "plapl(lower.tail)", "plapl(lower_tail")
    lower_tail <- lower.tail
  }
  if (lifecycle::is_present(log.p)) {
    lifecycle::deprecate_warn("0.4.0", "plapl(log.p)", "plapl(log_p)")
    log_p <- log.p
  }
  if (any(p < 0) || any(p > 1)) stop("p must lie between 0 and 1")
  if (any(scale < 0)) stop("scale must be non-negative")
  if (log_p) p <- exp(p)
  if (!lower_tail) p <- 1-p
  ret <- loc-scale*sign(p-0.5)*log(1-2*abs(p-0.5))
  nas <- is.nan(ret)
  if (sum(nas)) {
    p <- rep(p, length.out = length(ret))
    p <- p[nas]
    ret[nas] <- ifelse(p == 1, Inf, -Inf)
  }
  return(ret)
}

#' @rdname dlapl
#' @export
rlapl <- function(n, loc = 0, scale = 1) {
  if (any(scale < 0)) stop("scale must be non-negative")
  u <- runif (n)
  return(qlapl(u, loc = loc, scale = scale))
}
