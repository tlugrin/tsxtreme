## Copyright (C) 2017-2025 Thomas Lugrin
## Marginal transforms
## (Scope:) H+T/E+T fits, marginal transformations, thresholding
## List of functions: - format_ts
##                    - scale_ts
##                    - scale_to_margin
##                    - threshold_margin
##                    - lags_matrix
##################################################


##################################################
## FORMAT TIME SERIES FOR H+T MODEL

#' Pre-process time series to feed Heffernan--Tawn model
#' 
#' Called by [etfit()] and [et2fit()], not exposed to the user.
#' 
#' @param ts vector of reals, (stationary!) time series on an arbitrary scale.
#' @param u_mar scalar, marginal (GPD) threshold given as a probability.
#' @param u_dep scalar, conditional threshold given as a probability.
#' @param method string, "mle" for max likelihood or "mom" for the method of
#'   moments or "pwm" for proba weighted moments
#' @param nlag integer, number of lags to be considered.
#' @param lapl logical, is [ts] on Laplace scale already?
#' @returns A matrix with [nlag]+1 columns on the Laplace scale, with first
#'   column above threshold.
#' @keywords internal
format_ts <- function(ts, u_mar, u_dep, method = c("mle","mom","pwm"), nlag,
                      lapl = FALSE) {
  if (!lapl)
    ts <- scale_ts(ts = ts, u = u_mar, method = method)$ts_L
  data <- lags_matrix(ts = ts, nlag = nlag)
  data <- threshold_margin(data = data, u = u_dep, quantile = FALSE)
  return(data)
}


##################################################
## MARGINAL TRANSFORMATIONS

#' Marginal transformation to Laplace scale
#' 
#' Called by [format_ts()], not exposed to the user.
#' 
#' @param ts vector of reals, (stationary!) time series.
#' @param u scalar, threshold given as a probability.
#' @param method string, "mle" for max likelihood or "mom" for the method of
#'   moments or "pwm" for proba weighted moments.
#' @returns A list, original time series with different marginal distributions +
#'   GPD parameters.
#' @keywords internal
scale_ts <- function(ts, u, method = c("mle","mom","pwm")) {
  u    <- quantile(ts, u)
  n    <- length(ts)
  ts_U <- numeric(n)
  ei   <- (ts>u)
  p_u  <- sum(ei)/n
  pars <- gpd(ts, u, method)$pars
  if (pars[2] == 0) {
    ts_U[ei] <- 1 - p_u*exp(-(ts[ei]-u)/pars[1])
  }
  else{
    ts_U[ei]  <- 1 - p_u*(1+pars[2]/pars[1]*(ts[ei]-u))^(-1/pars[2])
  }
  ts_U[!ei] <- rank(ts[!ei])/(n+1)
  ts_L <- numeric(n)
  ts_L[ts_U<0.5]  <- log(2*ts_U[ts_U<0.5])
  ts_L[ts_U>=0.5] <- -log(2*(1-ts_U[ts_U>=0.5]))
  return(list(ts_U = ts_U, ts_G = -log(-log(ts_U)), ts_F = -1/log(ts_U),
              ts_L = ts_L, pars = c(pars, p_u)))
}


#' Use fitted transform to scale data to different scales
#' 
#' Currently not used, not exposed to the user.
#' 
#' @param x (vector of) reals, (a range of) observations for which change of
#'   scale is needed.
#' @param ts vector of reals, time series in original scale.
#' @param u scalar, threshold given as a probability.
#' @param gpd_pars vector of 3 reals, scale-shape-probability of exceeding u.
#' @returns A list, original observation(s) scaled to different marginal
#'   distributions: Uniform, Gumbel, Fréchet and Laplace.
#' @keywords internal
scale_to_margin <- function(x, ts, u, gpd_pars) {
  u   <- quantile(ts,u)
  n   <- length(ts)
  x_U <- x
  if (gpd_pars[2] == 0) {
    x_U[x>u] <- 1 - gpd_pars[3]*exp(-(x[x>u]-u)/gpd_pars[1])
  }
  else{
    x_U[x>u] <- 1 - gpd_pars[3]*(1+gpd_pars[2]/gpd_pars[1]*(x[x>u]-u))^(-1/gpd_pars[2])
  }
  x_U[x<=u] <- approx(sort(ts), y = (1:n)/(n+1), xout = x_U[x<=u],
                      method = "linear", ties = "ordered", rule = 2)$y
  if (x_U<0.5) x_L <- log(2*x_U)
  else x_L <- -log(2*(1-x_U))
  return(list(x_U = x_U, x_G = -log(-log(x_U)), x_F = -1/log(x_U), x_L = x_L))
}


#' Use fitted transform to scale probabilities back to the original scale
#' 
#' Called by [thetafit()], [chifit()], [theta2fit()]
#' 
#' @param p (vector of) reals, (a range of) probabilities which to transform to
#'   original [ts] scale.
#' @param ts vector of reals, time series in original scale.
#' @param u scalar, threshold given as a probability.
#' @param gpd_pars vector of 2 reals, scale-shape parameters of the GPD fitted
#'   to [ts].
#' @returns A vector of the same length as [p] on [ts] scale.
#' @keywords internal
scale_to_original <- function(p, ts, u, gpd_pars) {
  u_O <- quantile(ts, u)
  exc <- p > u
  p[exc]  <- qgpd((p[exc]-u)/(1-u), loc = u_O, scale = gpd_pars[1],
                  shape = gpd_pars[2])
  p[!exc] <- quantile(ts, p[!exc])
  return(p)
}



##################################################
## PRE-PROCESSING FOR HT & ET

#' Pre-process matrix of lagged data for Heffernan--Tawn fit
#' 
#' Called by [format_ts()], not exposed to the user.
#' 
#' @param data matrix of reals, replicates in rows - variables in columns.
#' @param u scalar, threshold given as a probability or quantile (cf. quantile).
#' @param quantile boolean, `TRUE` means `u` is given as a quantile --- `FALSE`
#'   means as a probability.
#' @returns A matrix of reals, only rows for which the first column of
#'  observations exceed `u` are kept.
#' @keywords internal
threshold_margin <- function(data, u = 0.95, quantile = FALSE)
{
  if (!quantile) u <- quantile(data[,1], u)
  ei      <- (data[,1]>u)
  return(data[ei,])
}


#' Build matrix of lagged observations
#' 
#' Called by [format_ts()] and [thetaruns()], not exposed to the user.
#' 
#' @param ts vector of reals, (stationary!) time series.
#' @param nlag scalar, number of lags of the output lagged matrix
#'   (=> number of columns == `nlag`+1).
#' @returns A matrix of lagged versions of `ts` in its columns.
#' @keywords internal
lags_matrix <- function(ts, nlag) {
  n      <- length(ts)
  dim    <- nlag+1
  ts_mat <- matrix(0, nrow = n-dim+1, ncol = dim)
  for(i in 1:dim) {
    ts_mat[,i] <- ts[i:(n-dim+i)]
  }
  return(ts_mat)
}
