## Copyright (C) 2017-2025 Thomas Lugrin
## runs estimator for theta(x,m)
## (Scope:) empirical estimate to be compared to E+T typically
## List of functions: - runs
##                    - compute_runs
##                    - block_bootstrap_sample
##################################################

##################################################
#' Runs estimator
#' 
#' Compute the empirical estimator of the extremal index using the runs method
#' (Smith & Weissman, 1994, JRSSB).
#' 
#' @details
#' Consider a stationary time series \eqn{(X_t)}. A characterisation of the
#' extremal index is
#' \deqn{\theta(x,m) = Pr(X_1\le x,\ldots,X_m\le x \mid X_0\ge x).}{\theta(x,m) = Pr(X_1\le x,\dots,X_m\le x | X_0\ge x).}
#' In the limit when \eqn{x} and \eqn{m} tend to \eqn{\infty} appropriately,
#' \eqn{\theta} corresponds to the asymptotic inverse mean cluster size. It also
#' links the generalised extreme value distribution of the independent series
#' \eqn{(Y_t)}, with the same marginal distribution as \eqn{(X_t)},
#' \deqn{G_Y(z)=G_X^\theta(z),}
#' with \eqn{G_X} and \eqn{G_Y} the extreme value distributions of \eqn{(X_t)}
#' and \eqn{(Y_t)} respectively.
#' 
#' `nlag` corresponds to the _run-length_ \eqn{m} and `probs` is a set of values
#' for \eqn{x}.
#' The _runs_ estimator is computed, which consists of counting the proportion
#' of clusters to the number of exceedances of a threshold \eqn{x}; two
#' exceedances of the threshold belong to different clusters if there are at
#' least \eqn{m+1} non-exceedances inbetween.
#' 
#' @param ts a vector, the time series for which to estimate the threshold-based
#'   extremal index \eqn{\theta(x,m)}, with \eqn{x} a probability level and
#'   \eqn{m} a run-length (see details).
#' @param lapl logical; is `ts` on the Laplace scale already? The default
#'   (FALSE) assumes unknown marginal distribution.
#' @param nlag the run-length; an integer larger or equal to 1.
#' @param u.mar `r lifecycle::badge('deprecated')` use `u_mar` instead.
#' @param u_mar marginal threshold (probability); used when transforming the
#'   time series to Laplace scale if `lapl` is `FALSE`; if `lapl` is `TRUE`,
#'   it is nevertheless used when bootstrapping, since the bootstrapped series
#'   generally do not have Laplace marginal distributions.
#' @param probs vector of probabilities; the values of `x` for which to evaluate
#'   \eqn{\theta(x,m)}.
#' @param method.mar `r lifecycle::badge('deprecated')` use `method_mar` instead.
#' @param method_mar a character string defining the method used to estimate the
#'   marginal GPD; either `"mle"` for maximum likelihood or `"mom"` for method
#'   of moments or `"pwm"` for probability weighted moments methods. Defaults to
#'   `"mle"`.
#' @param block.length `r lifecycle::badge('deprecated')` use `block_length`
#'   instead.
#' @param block_length integer; the block length used for the block-bootstrapped
#'   confidence intervals.
#' @param R.boot `r lifecycle::badge('deprecated')` use `R_boot`instead.
#' @param R_boot integer; the number of samples used for the block bootstrap.
#' @param levels vector of probabilites; the quantiles of the posterior
#'   distribution of the extremal index \eqn{\theta(x,m)} to output.
#' @returns An object of class [depmeasure()] containing:
#'   \item{theta }{matrix; estimates of the extremal index \eqn{\theta(x,m)} with rows corresponding to the \code{probs} values of \eqn{x} and the columns to the runs estimate and the chosen \code{levels}-quantiles of the bootstrap distribution.}
#'   \item{nbr_exc }{numeric vector; number of exceedances for each threshold corresponding to the elements in \code{probs}.}
#'   \item{probs }{\code{probs}.}
#'   \item{levels }{numeric vector; \code{probs} converted to the original scale of \code{ts}.}
#'   \item{nlag }{\code{nlag}.}
#' @seealso [theta2fit()], [thetafit()]
#' @examples
#' ## generate data from an AR(1)
#' ## with Gaussian marginal distribution
#' n   <- 10000
#' dep <- 0.5
#' ar    <- numeric(n)
#' ar[1] <- rnorm(1)
#' for (i in 2:n)
#'   ar[i] <- rnorm(1, mean = dep*ar[i-1], sd = 1-dep^2)
#' ## transform to Laplace scale
#' ar <- qlapl(pnorm(ar))
#' ## compute empirical estimate
#' theta <- thetaruns(ts=ar, u_mar = 0.95, probs = c(0.95, 0.98, 0.99))
#' ## output
#' plot(theta, ylim=c(.2,1))
#' abline(h=1, lty="dotted")
#' @export
#' @aliases empfit
thetaruns <- function(ts,
                      lapl = FALSE,
                      nlag = 1,
                      u_mar = 0,
                      probs = seq(u_mar, 0.995, length.out = 30),
                      method_mar = c("mle", "mom", "pwm"),
                      R_boot = 0,
                      block_length = (nlag + 1)*5,
                      levels = c(.025, .975),
                      u.mar = deprecated(),
                      method.mar = deprecated(),
                      R.boot = deprecated(),
                      block.length = deprecated()) {
  if (lifecycle::is_present(u.mar)) {
    lifecycle::deprecate_warn("0.4.0", "thetaruns(u.mar)", "thetaruns(u_mar)")
    u_mar <- u.mar
  }
  if (lifecycle::is_present(method.mar)) {
    lifecycle::deprecate_warn("0.4.0", "thetaruns(method.mar)",
                              "thetaruns(method_mar)")
    method_mar <- method.mar
  }
  if (lifecycle::is_present(R.boot)) {
    lifecycle::deprecate_warn("0.4.0", "thetaruns(R.boot)", "thetaruns(R_boot)")
    R_boot <- R.boot
  }
  if (lifecycle::is_present(block.length)) {
    lifecycle::deprecate_warn("0.4.0", "thetaruns(block.length)",
                              "thetaruns(block_length)")
    block_length <- block.length
  }
  ## marginal re-scale to Laplace if needed
  if (!lapl)
    ts_L <- scale_ts(ts = ts, u = u_mar, method = method_mar)$ts_L
  data <- lags_matrix(ts=ts_L, nlag=nlag)
  mesh <- probs
  mesh_O <- scale_to_original(p = mesh, ts = ts, u = u_mar,
                              gpd_pars = gpd(ts, u=u_mar, method=method_mar)$pars)
  mesh_L <- qexp(1-(1-mesh)*2)
  ret <- depmeasure("runs")
  ret$probs <- probs
  ret$levels <- mesh_O
  ret$nlag <- nlag
  ## set containers
  nbr_vert <- length(mesh)
  nbr_quant<- 1+length(levels)
  theta    <- matrix(0, nrow = nbr_vert, ncol = nbr_quant)
  colnames(theta) <- c("estimate", paste(levels*100,"%", sep=""))
  ## compute and store
  th_est             <- compute_runs(data, mesh_L)
  theta[,"estimate"] <- th_est[[1]]
  nbr_exc            <- th_est[[2]]
  ret$nbr_exc <- nbr_exc
  ## compute bootstrapped quantiles
  if (R_boot > 0 & nbr_quant > 1) {
    ts_R <- block_bootstrap_sample(ts, block_length, R_boot)
    th_R <- matrix(0, nrow = nbr_vert, ncol=R_boot)
    for (r in 1:R_boot) {
      ts_boot  <- scale_ts(ts_R[,r], u = u_mar, method = method_mar)$ts_L
      ts_mat   <- lags_matrix(ts_boot, nlag = nlag)
      th_R[,r] <- compute_runs(ts_mat, mesh_L)[[1]]
    }
    theta[,-1] <- t(apply(th_R, 1, quantile, levels))
  }else{
    theta <- theta[,"estimate",drop=FALSE]
  }
  ret$theta <- theta
  return(ret)
}


#' Compute the extremal index using the runs method
#' 
#' Called by [thetaruns()].
#' 
#' @param data matrix, lags matrix of a time series
#' @param mesh vector, x in theta(x,m) on Laplace scale
#' @returns list, estimates of theta across `mesh` - nbr of exceedances
#' across `mesh`.
#' @keywords internal
compute_runs <- function(data, mesh) {
  nbr_vert <- length(mesh)
  theta <- nbr_exc <- numeric(nbr_vert)
  for (i in 1:nbr_vert) {
    exc        <- data[data[,1]>=mesh[i],,drop=FALSE]
    nbr_exc[i] <- dim(exc)[1]
    theta[i]   <- sum(apply(exc[,-1,drop=FALSE], 1, max)<mesh[i])/nbr_exc[i]
  }
  return(list(theta,nbr_exc))
}

#' Compute a block bootstrap sample from a time series
#' 
#' Called by [thetaruns()].
#' 
#' @param ts vector of reals, time series which to estimate theta(x,m) for
#' @param block_length: integer>0, used for the block-bootstrap CI
#' @param R integer>0, nbr of repetitions used for the block-bootstrap CI
#' @returns matrix, each column is a block-bootstrapped time series
#' @keywords internal
block_bootstrap_sample <- function(ts, block_length, R) {
  ## handle cases when block length does not divide the time series
  rest <- length(ts)%%block_length
  if (rest > 0)
    ts <- ts[-((length(ts)-rest+1):length(ts))]
  ts_mat <- matrix(ts, nrow=block_length, byrow=FALSE)
  nbr_bl <- floor(length(ts)/block_length)# == dim(ts_mat)[2]
  ts_R   <- matrix(0, nrow=length(ts), ncol=R)
  for (r in 1:R) {
    sample   <- sample.int(nbr_bl, size=nbr_bl, replace=TRUE)
    ts_R[,r] <- ts_mat[,sample]
  }
  return(ts_R)
}
