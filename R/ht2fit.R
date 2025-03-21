## Copyright (C) 2017-2025 Thomas Lugrin
## Gathering functions from old files into one reference file
## (Scope:) Heffernan-Tawn model fit in 2 stages
## List of functions: - ht2fit
##                    - nllh_ht
##                    - grad_ht
##                    - ht2step_2d
##                    - ht2step
##################################################

#' Dependence model fit (stepwise)
#' 
#' The conditional Heffernan--Tawn model is used to fit the dependence in time
#' of a stationary series. A standard 2-stage procedure is used.
#' 
#' @details
#' Consider a stationary time series \eqn{(X_t)} with Laplace marginal
#' distribution; the fitting procedure consists of fitting
#' \deqn{X_t = \alpha_t\times x_0 + x_0^{\beta_t}\times Z_t,\quad t=1,\ldots,m,}{X_t = \alpha_t * x_0 + x_0^{\beta_t} * Z_t, t=1,\ldots,m,}
#' with \eqn{m} the number of lags considered. A likelihood is maximised
#' assuming \eqn{Z_t\sim N(\mu_t, \sigma^2_t)}{Z_t ~ N(\mu_t, \sigma^2_t)},
#' then an empirical distribution for the \eqn{Z_t} is derived using the
#' estimates of \eqn{\alpha_t} and \eqn{\beta_t} and the relation
#' \deqn{\hat Z_t = \frac{X_t - \hat\alpha_t\times x_0}{x_0^{\hat\beta_t}}.}{Z_t = (X_t - \alpha_t * x_0) / x_0^{\beta_t}.}
#' 
#' [conditions] implements additional conditions suggested by
#' Keef, Papastathopoulos and Tawn (2013) on the ordering of conditional
#' quantiles. These conditions help with getting a consistent fit by shrinking
#' the domain in which \eqn{(\alpha,\beta)} live.
#' 
#' @param ts numeric vector; time series to be fitted.
#' @param u.mar `r lifecycle::badge('deprecated')` use `u_mar` instead.
#' @param u_mar marginal threshold; used when transforming the time series to
#'   Laplace scale.
#' @param u.dep `r lifecycle::badge('deprecated')` use `u_dep` instead.
#' @param u_dep dependence threshold; level above which the dependence is
#'   modelled. `u_dep` can be lower than `u_mar`.
#' @param lapl logical; is \code{ts} on the Laplace scale already? The default
#'   (FALSE) assumes unknown marginal distribution.
#' @param method.mar `r lifecycle::badge('deprecated')` use `method_mar`
#'   instead.
#' @param method_mar a character string defining the method used to estimate the
#'   marginal GPD; either \code{"mle"} for maximum likelihood of \code{"mom"}
#'   for method of moments. Defaults to \code{"mle"}.
#' @param nlag integer; number of lags to be considered when modelling the
#'   dependence in time.
#' @param conditions logical; should conditions on \eqn{\alpha} and \eqn{\beta}
#'   be set? (see Details) Defaults to `TRUE`.
#' @returns An object of class [stepfit()].
#' @seealso [depfit()], [theta2fit()]
#' @examples
#' ## generate data from an AR(1)
#' ## with Gaussian marginal distribution
#' n   <- 10000
#' dep <- 0.5
#' ar    <- numeric(n)
#' ar[1] <- rnorm(1)
#' for(i in 2:n)
#'   ar[i] <- rnorm(1, mean=dep*ar[i-1], sd=1-dep^2)
#' plot(ar, type="l")
#' plot(density(ar))
#' grid <- seq(-3,3,0.01)
#' lines(grid, dnorm(grid), col="blue")
#' 
#' ## rescale margin
#' ar <- qlapl(pnorm(ar))
#' 
#' ## fit model without constraints...
#' fit1 <- dep2fit(ts=ar, u_mar = 0.95, u_dep=0.98, conditions=FALSE)
#' fit1$a; fit1$b
#' 
#' ## ...and compare with a fit with constraints
#' fit2 <- dep2fit(ts=ar, u_mar = 0.95, u_dep=0.98, conditions=TRUE)
#' fit2$a; fit2$b# should be similar, as true parameters lie well within the constraints
#' 
#' @export
dep2fit <- function(ts,
                    u_mar = 0,
                    u_dep,
                    lapl = FALSE,
                    method_mar = c("mle","mom","pwm"),
                    nlag = 1,
                    conditions = TRUE,
                    u.mar = deprecated(),
                    u.dep = deprecated(),
                    method.mar = deprecated()) {
  if (lifecycle::is_present(u.mar)) {
    lifecycle::deprecate_warn("0.4.0", "dep2fit(u.mar)", "dep2fit(u_mar)")
    u_mar <- u.mar
  }
  if (lifecycle::is_present(u.dep)) {
    lifecycle::deprecate_warn("0.4.0", "dep2fit(u.dep)", "dep2fit(u_dep)")
    u_dep <- u.dep
  }
  if (lifecycle::is_present(method.mar)) {
    lifecycle::deprecate_warn("0.4.0", "dep2fit(method.mar)",
                              "dep2fit(method_mar)")
    method_mar <- method.mar
  }
  data_up <- format_ts(ts = ts, u_mar = u_mar, u_dep = u_dep,
                       method = method_mar, lapl = lapl, nlag = nlag)
  ht2step(data = data_up, conditions = conditions)
}

##################################################
## LIKELIHOODS

#' Negative log-likelihood for the Heffernan--Tawn model
#' 
#' Called by [ht2step_2d()]
#' 
#' @param par vector of reals, alpha-beta-mu-sigma^2
#' @param data 2-column matrix, first column is X>u - second column is Y|X>u
#' @param conditions boolean, TRUE means alpha and beta must satisfy AD/AI constraints
#' @returns scalar, value of the negative log-likelihood of H+T
#' @keywords internal
nllh_ht <- function(par,data,conditions) {
  if (par[4] <= 0) { return(Inf) }
  if (conditions) {
    if (!conditions_verify(par[1], par[2], 0, data=data) ||
         !conditions_verify(par[1], par[2], 1, data=data)) { return(Inf) }
  } else {
    if (par[1]< -1 || par[1]>1 || par[2]<0 || par[2]>1) { return(Inf) }
  }
  sigma <- par[4]*data[,1]^(2*par[2])
  mu    <- par[1]*data[,1] + par[3]*data[,1]^par[2]
  return(sum(log(sigma) + (data[,2]-mu)^2/sigma))
}


#' @rdname nllh_ht
#' @keywords internal
grad_ht <- function(par,data,conditions) {
  sig <- par[4]*data[,1]^(2*par[2])
  mu  <- par[1]*data[,1] + par[3]*data[,1]^par[2]
  cen <- data[,2]-mu
  
  d_a <- -2*sum(cen*data[,1]/sig)
  d_b <- 2*sum(log(data[,1]) * (1 - cen*(par[3]*data[,1]^par[2]+cen)/sig))
  d_m <- -2*sum(cen/(data[,1]^(par[2])*par[4]))
  d_s <- sum((1-cen^2/sig)/par[4])
  
  return(c(d_a,d_b,d_m,d_s))
}


##################################################
## FITS

#' Heffernan--Tawn model stepwise fit in two dimensions
#'
#' Called by [ht2step()]
#' 
#' @param data 2-column matrix, first column is \eqn{X>u} (Laplace) - second
#'   column is \eqn{Y|X>u}
#' @param conditions boolean, TRUE means alpha and beta must satisfy AD/AI constraints
#' @returns list, MLE for alpha - beta - z hat - standard deviations of alpha and beta
#' @keywords internal
ht2step_2d <- function(data, conditions = TRUE) {
  ret <- stepfit()
  if (conditions) {
    bds <- conditions_bounds(0,FALSE,data)
    a   <- runif (1,bds[1],bds[2])
    bds <- conditions_bounds(a,TRUE,data)
    b   <- runif (1,bds[1],bds[2])
  } else {
    a <- runif (1, -1, 1)
    b <- runif (1, 0, 1)
  }
  res <- optim(par = c(a,b,0,1), fn = nllh_ht, gr = grad_ht, data = data,
               conditions = conditions, hessian = TRUE, method = "BFGS")
  ret$a <- res$par[1]
  ret$b <- res$par[2]
  ret$res <- (data[,2] - data[,1]*ret$a)/data[,1]^ret$b
  ret$pars_se <- sqrt(diag(solve(res$hessian))[1:2])
  return(ret)
}


#' Heffernan--Tawn model stepwise fit in the general case
#' 
#' Called by [dep2fit()]
#' 
#' @param data matrix of reals, first column is X>u (Laplace) - other columns are |X>u
#' @param conditions boolean, TRUE means alpha and beta must satisfy AD/AI constraints
#' @returns An object of class [stepfit()].
#' @keywords internal
ht2step <- function(data, conditions = TRUE) {
  n       <- dim(data)[1]
  d       <- dim(data)[2]
  ret   <- stepfit()
  ret$a <- ret$b <- numeric(d-1)
  ret$res     <- matrix(0, nrow = n, ncol = d-1,
                        dimnames = list(NULL, paste("Z",1:(d-1), sep = "")))
  ret$pars_se <- matrix(NA, nrow = 2, ncol = d-1)
  ret$nlag    <- d-1
  for(i in 2:d) {
    fit           <- ht2step_2d(data[,c(1,i)], conditions)
    ret$a[i-1]    <- fit$a
    ret$b[i-1]    <- fit$b
    ret$res[,i-1]     <- fit$res
    ret$pars_se[,i-1] <- fit$pars_se
  }
  return(ret)
}
