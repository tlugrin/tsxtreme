## Copyright (C) 2017-2025 Thomas Lugrin
## Marginal transforms
## (Scope:) H+T/E+T fits, marginal transformations, thresholding
## List of functions: - format.ts
##                    - scale.ts
##                    - scale.to.margin
##                    - threshold.margin
##                    - lags.matrix
##################################################


##################################################
## FORMAT TIME SERIES FOR H+T MODEL

#' Pre-process time series to feed Heffernan--Tawn model
#' 
#' Called by [etfit()] and [et2fit()], not exposed to the user.
#' 
#' @param ts vector of reals, (stationary!) time series on an arbitrary scale.
#' @param u.mar scalar, marginal (GPD) threshold given as a probability.
#' @param u.dep scalar, conditional threshold given as a probability.
#' @param method string, "mle" for max likelihood or "mom" for the method of
#'   moments or "pwm" for proba weighted moments
#' @param nlag integer, number of lags to be considered.
#' @param lapl logical, is [ts] on Laplace scale already?
#' @returns A matrix with [nlag]+1 columns on the Laplace scale, with first
#'   column above threshold.
#' @keywords internal
format.ts <- function(ts, u.mar, u.dep, method=c("mle","mom","pwm"), nlag, lapl=FALSE){
  if(!lapl)
    ts <- scale.ts(ts=ts, u=u.mar, method=method)$ts.L
  data <- lags.matrix(ts=ts, nlag=nlag)
  data <- threshold.margin(data=data, u=u.dep, quantile=FALSE)
  return(data)
}


##################################################
## MARGINAL TRANSFORMATIONS

#' Marginal transformation to Laplace scale
#' 
#' Called by [format.ts()], not exposed to the user.
#' 
#' @param ts vector of reals, (stationary!) time series.
#' @param u scalar, threshold given as a probability.
#' @param method string, "mle" for max likelihood or "mom" for the method of
#'   moments or "pwm" for proba weighted moments.
#' @returns A list, original time series with different marginal distributions +
#'   GPD parameters.
#' @keywords internal
scale.ts <- function(ts, u, method=c("mle","mom","pwm")){
  u    <- quantile(ts, u)
  n    <- length(ts)
  ts.U <- numeric(n)
  ei   <- (ts>u)
  p.u  <- sum(ei)/n
  pars <- gpd(ts, u, method)$pars
  if(pars[2] == 0){
    ts.U[ei] <- 1 - p.u*exp(-(ts[ei]-u)/pars[1])
  }
  else{
    ts.U[ei]  <- 1 - p.u*(1+pars[2]/pars[1]*(ts[ei]-u))^(-1/pars[2])
  }
  ts.U[!ei] <- rank(ts[!ei])/(n+1)
  ts.L <- numeric(n)
  ts.L[ts.U<0.5]  <- log(2*ts.U[ts.U<0.5])
  ts.L[ts.U>=0.5] <- -log(2*(1-ts.U[ts.U>=0.5]))
  return(list(ts.U=ts.U, ts.G=-log(-log(ts.U)), ts.F=-1/log(ts.U), ts.L=ts.L, pars=c(pars,p.u)))
}


#' Use fitted transform to scale data to different scales
#' 
#' Currently not used, not exposed to the user.
#' 
#' @param x (vector of) reals, (a range of) observations for which change of
#'   scale is needed.
#' @param ts vector of reals, time series in original scale.
#' @param u scalar, threshold given as a probability.
#' @param gpd.pars vector of 3 reals, scale-shape-probability of exceeding u.
#' @returns A list, original observation(s) scaled to different marginal
#'   distributions: Uniform, Gumbel, Fréchet and Laplace.
#' @keywords internal
scale.to.margin <- function(x, ts, u, gpd.pars){
  u   <- quantile(ts,u)
  n   <- length(ts)
  x.U <- x
  if(gpd.pars[2] == 0){
    x.U[x>u] <- 1 - gpd.pars[3]*exp(-(x[x>u]-u)/gpd.pars[1])
  }
  else{
    x.U[x>u] <- 1 - gpd.pars[3]*(1+gpd.pars[2]/gpd.pars[1]*(x[x>u]-u))^(-1/gpd.pars[2])
  }
  x.U[x<=u] <- approx(sort(ts), y=(1:n)/(n+1), xout=x.U[x<=u], method="linear", ties="ordered", rule=2)$y
  if(x.U<0.5) x.L <- log(2*x.U)
  else x.L <- -log(2*(1-x.U))
  return(list(x.U=x.U, x.G=-log(-log(x.U)), x.F=-1/log(x.U), x.L=x.L))
}


#' Use fitted transform to scale probabilities back to the original scale
#' 
#' Called by [thetafit()], [chifit()], [theta2fit()]
#' 
#' @param p (vector of) reals, (a range of) probabilities which to transform to
#'   original [ts] scale.
#' @param ts vector of reals, time series in original scale.
#' @param u scalar, threshold given as a probability.
#' @param gpd.pars vector of 2 reals, scale-shape parameters of the GPD fitted
#'   to [ts].
#' @returns A vector of the same length as [p] on [ts] scale.
#' @keywords internal
scale.to.original <- function(p, ts, u, gpd.pars){
  u.O <- quantile(ts, u)
  exc <- p > u
  p[exc]  <- qgpd((p[exc]-u)/(1-u), loc=u.O, scale=gpd.pars[1], shape=gpd.pars[2])
  p[!exc] <- quantile(ts, p[!exc])
  return(p)
}



##################################################
## PRE-PROCESSING FOR HT & ET

#' Pre-process matrix of lagged data for Heffernan--Tawn fit
#' 
#' Called by [format.ts()], not exposed to the user.
#' 
#' @param data matrix of reals, replicates in rows - variables in columns.
#' @param u scalar, threshold given as a probability or quantile (cf. quantile).
#' @param quantile boolean, `TRUE` means `u` is given as a quantile --- `FALSE`
#'   means as a probability.
#' @returns A matrix of reals, only rows for which the first column of
#'  observations exceed `u` are kept.
#' @keywords internal
threshold.margin <- function(data, u=0.95, quantile=FALSE)
{
  if(!quantile) u <- quantile(data[,1], u)
  ei      <- (data[,1]>u)
  return(data[ei,])
}


#' Build matrix of lagged observations
#' 
#' Called by [format.ts()] and [thetaruns()], not exposed to the user.
#' 
#' @param ts vector of reals, (stationary!) time series.
#' @param nlag scalar, number of lags of the output lagged matrix
#'   (=> number of columns == `nlag`+1).
#' @returns A matrix of lagged versions of `ts` in its columns.
#' @keywords internal
lags.matrix <- function(ts, nlag){
  n      <- length(ts)
  dim    <- nlag+1
  ts.mat <- matrix(0, nrow=n-dim+1, ncol=dim)
  for(i in 1:dim){
    ts.mat[,i] <- ts[i:(n-dim+i)]
  }
  return(ts.mat)
}
