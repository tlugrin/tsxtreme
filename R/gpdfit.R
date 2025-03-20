## Copyright (C) 2017-2025 Thomas Lugrin
## Generalised Pareto distribution fit
## (Scope:) GPD fit, marginal transform
## List of functions: - nllh_gpd
##                    - grad_gpd
##                    - nllh_gpd_0
##                    - gpd.mle
##                    - gpd.mom
##                    - gpd.pwm
##################################################

##################################################
## LIKELIHOOD-BASED GPD FIT
#' Negative log-likelihood and its gradient for the GPD
#' 
#' Called by [gpd.mle()], not exposed to the user.
#' 
#' @param par vector of 2 reals, scale-shape.
#' @param u scalar, threshold given as a quantile.
#' @param x vector of reals, data on which to fit the GPD.
#' @returns A scalar, value of the negative log-likelihood of the GPD.
#' @keywords internal
nllh_gpd <- function(par, u, x){
  exced <- x[x>u]-u
  n     <- length(exced)
  if(par[1] <= 0) return(Inf)
  if(abs(par[2]) < 0.00001){
    nllh     <- n*log(par[1])+sum(exced)/par[1]
  }
  else{
    pos_part <- pmax(1+par[2]/par[1]*exced, 0)
    if(min(pos_part) <= 0) return(Inf)
    nllh     <- n*log(par[1])+(1/par[2]+1)*sum(log(pos_part))
  }
  return(nllh)
}

#' @rdname nllh_gpd
#' @keywords internal
grad_gpd <- function(par, u, x){
  exced <- x[x>u]-u
  n     <- length(exced)
  if(abs(par[2]) < 0.00001){
    d_s <- n/par[1] - sum(exced)/par[1]^2
    d_x <- 0
  }
  else{
    d_s <- ( n-(1/par[2]+1)*sum(1/(1+par[1]/(par[2]*exced))) )/par[1]
    d_x <- -1/par[2]^2*sum(log(1+par[2]*exced/par[1]) + (1/par[2]+1)*sum(1/(par[1]/exced+par[2])))
  }
  return(c(d_s,d_x))
}

#' Negtive log-likelihood for the GPD with null shape
#' 
#' Called by [gpd.mle()] when `shape_null_hyp==TRUE`, not exposed to the user.
#' 
#' @param par scalar, scale (shape is assumed 0).
#' @param u scalar, threshold given as a quantile.
#' @param x vector of reals, data on which to fit the GPD.
#' @returns A scalar, value of the negative log-likelihood of the GPD
#' @keywords internal
nllh_gpd_0 <- function(par, u, x){
  nllh_gpd(c(par,0), u, x)
}


#' Maximum likelihood estimation for the GPD
#' 
#' Called by [scale_ts()], not exposed to the user.
#' 
#' @param ts vector of reals, stationary time series.
#' @param u scalar, threshold given as a quantile.
#' @param hessian boolean, TRUE means computing the hessian matrix of the
#'   likelihood.
#' @returns A list, threshold - MLEs - log-likelihood - standard deviations -
#'   covariance matrix.
#' @keywords internal
gpd.mle <- function(ts, u, hessian=FALSE){
  opt <- optim(c(sd(ts),0.1), fn=nllh_gpd, gr=grad_gpd, u=u, x=ts,
               hessian=hessian)
  ret <- list(u=u, pars=opt$par, llh=-opt$value)
  if(hessian){
    ret$var     <- solve(opt$hessian)
    ret$pars_sd <- sqrt(diag(var))
  }
  return(ret)
}

##################################################
## MOMENT-BASED GPD FIT
#' Method of moments and probability weighted moments to fit the GPD
#' 
#' Called by [scale_ts()], not exposed to the user.
#' 
#' @param ts vector of reals, stationary time series.
#' @param  u scalar, threshold given as a quantile.
#' @returns A list, threshold - MOM estimates.
#' @keywords internal
gpd.mom <- function(ts, u){
  ts <- ts-u
  m <- mean(ts)
  s <- var(ts)
  scale <- m*(m^2/s+1)/2
  shape <- (m^2/s-1)/2
  return(list(u = u, pars = c(scale, shape)))
}

#' @rdname gpd.mom
#' @keywords internal
gpd.pwm <- function(ts, u){
  ts <- ts-u
  ts <- sort(ts, decreasing = FALSE)
  n  <- length(ts)
  a0 <- mean(ts)
  a1 <- sum(ts[-n]*((n-1):1))/(n-1)/n
  scale <- 2*a0*a1/(a0-2*a1)
  shape <- a0/(a0-2*a1)-2
  return(list(u = u, pars=c(scale, shape)))
}


##################################################
## GENERIC FUNCTION
#' GPD fit wrapper
#' 
#' @inheritParams gpd.mom
#' @keywords internal
gpd <- function(ts, u, method){
  fct_name <- paste("gpd.", method[1], sep="")
  eval(call(fct_name, ts = ts, u = u))
}
