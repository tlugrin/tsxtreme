## Copyright (C) 2017-2025 Thomas Lugrin
## user inteface for htfit class
## (Scope:) Fit H+T model using Bayesian semiparametrics
## List of functions: - depfit
##                    - htfit
##################################################

##################################################
## H+T R interface
#' Conditional tail model Bayesian fit
#' 
#' Bayesian semiparametrics are used to fit the Eastoe--Tawn model to a time
#' series. Options are available to impose a structure in time on the model.
#' 
#' @details
#' \code{submodel} can be \code{"fom"} to impose a first order Markov structure
#' on the model parameters \eqn{\alpha_j} and \eqn{\beta_j}
#' (see [thetafit()] for more details); it can take the value \code{"none"} to
#' impose no particular structure in time; it can also be \code{"ugm"} which
#' can be applied to density estimation, as it corresponds to setting
#' \eqn{\alpha=\beta=0} (see examples).#' 
#' 
#' @param ts numeric vector; time series to be fitted.
#' @param u.mar marginal threshold; used when transforming the time series to
#'   the Laplace scale.
#' @param u.dep dependence threshold; level above which the dependence is
#'   modelled. \code{u.dep} can be lower than \code{u.mar}.
#' @param lapl logical; is \code{ts} on the Laplace scale already? The default
#'   (FALSE) assumes unknown marginal distribution.
#' @param method.mar a character string defining the method used to estimate the
#'   marginal GPD; either \code{"mle"} for maximum likelihood or \code{"mom"}
#'   for method of moments or \code{"pwm"} for probability weighted moments.
#'   Defaults to \code{"mle"}.
#' @param nlag integer; number of lags to be considered when modelling the
#'   dependence in time.
#' @param par an object of class "bayesparams".
#' @param submodel a character string; "fom" for \emph{first order Markov},
#'   "none" for \emph{no particular time structure}, or "ugm" for
#'   _univariate Gaussian mixture_ (see Details).
#' @returns An object of class [bayesfit()].
#' @seealso [thetafit()], [chifit()]
#' @examples
#' ## generate data from an AR(1)
#' ## with Gaussian marginal distribution
#' n   <- 10000
#' dep <- 0.5
#' ar    <- numeric(n)
#' ar[1] <- rnorm(1)
#' for(i in 2:n)
#'   ar[i] <- rnorm(1, mean=dep*ar[i-1], sd=1-dep^2)
#'
#' ## rescale the margin
#' ar <- qlapl(pnorm(ar))
#' 
#' ## fit the data
#' params <- bayesparams()
#' params$maxit <- 100# bigger numbers would be
#' params$burn  <- 10 # more sensible...
#' params$thin  <- 4
#' fit <- depfit(ts=ar, u.mar=0.95, u.dep=0.98, par=params)
#'
#' ########
#' ## density estimation with submodel=="ugm"
#' data <- MASS::galaxies/1e3
#' dens <- depfit(ts=data, par=params, submodel="ugm")
#' @export
depfit <- function(ts, u.mar=0, u.dep=u.mar, lapl=FALSE, method.mar=c("mle","mom","pwm"), nlag=1,
                   par=bayesparams(), submodel=c("fom","none","ugm")){
  if(submodel[1]=="ugm"){
    data.up <- cbind(1,ts)
  }else if(submodel[1] %in% c("fom","none")){
    data.up <- format.ts(ts=ts, u.mar=u.mar, u.dep=u.dep, method=method.mar,
                         nlag=nlag, lapl=lapl)
  }else{
    stop(paste("submodel",submodel,"unknown in depfit()"))
  }
  htfit(data=data.up, prop.a=par$prop.a, prop.b=par$prop.b,
        prior.mu=par$prior.mu, prior.nu=par$prior.nu, prior.eta=par$prior.eta,
        trunc=par$trunc, comp.saved=par$comp.saved,
        maxit=par$maxit, burn=par$burn, thin=par$thin,
        adapt=par$adapt, batch.size=par$batch.size, start.ab=par$start.ab,
        mode=par$mode, submodel=submodel)
}

#' Raw method for the semiparametric Bayesian fit
#' 
#' Wrapper function for `C_et_interface()`. Called by [depfit()], not exposed
#' to the user.
#' 
#' @param data matrix, first column is > threshold -- other columns are lagged
#'   versions of the first one.
#' @inheritParams bayesparams
#' @inheritParams depfit
#' @returns An object of class [bayesfit()].
#' @keywords internal
htfit <- function(data,
                  prop.a, prop.b,
                  prior.mu=c(0,10), prior.nu=c(2,1/2), prior.eta=c(4,1),
                  trunc=100, comp.saved=10,
                  maxit=10000, burn=2000, thin=4,
                  adapt=2000, batch.size=125,
                  start.ab = c("guesstimate", "prior"),
                  mode=1, submodel=c("fom","none","ugm")){# "ugm": univariate Gaussian mixture
  n    <- dim(data)[1]
  nlag <- dim(data)[2]-1
  if(length(prop.a) < 5){ prop.a <- rep(prop.a[1],5) }
  if(length(prop.b) < 3){ prop.b <- rep(prop.b[1],3) }
  ## build trace containers
  tr.len <- floor((maxit-burn)/thin)
  t.a   <- t.b <- matrix(0, nrow=tr.len, ncol=nlag)
  t.sig <- t.mu <- array(0, dim = c(tr.len,comp.saved,nlag))
  t.w   <- matrix(0, nrow=tr.len, ncol=comp.saved)
  t.g   <- numeric(tr.len)
  t.ci  <- matrix(0, nrow=tr.len, ncol=n)
  t.noo <- matrix(0, nrow=tr.len, ncol=comp.saved)
  t.noc <- numeric(tr.len)
  t.sd  <- array(0, dim=c(tr.len,8,nlag))#for RAMA
  if(start.ab[1] == "guesstimate"){
    fit2step <- ht2step(data) # lags > 1 will be properly initialised in C routine
    start.a <- fit2step$a     # depending on [submodel]
    start.b <- fit2step$b
  }else if(start.ab[1] == "prior"){
    start.a <- start.b <- rep(Inf,nlag) # will be properly initialised in C routine
  }else{
    stop("In htfit(): invalid starting value strategy for start.ab")
  }
  if(submodel[1] == "ugm"){ submodel <- 0 }
  else if(submodel[1]=="fom"){ submodel <- 1 }
  else if(submodel[1]=="none"){ submodel <- 2 }
  else{ stop("In htfit(): invalid submodel.") }
  ## call C(++) wrapper
  fit = .C(C_et_interface, as.double(data), as.integer(n), as.integer(nlag),
           as.integer(trunc), as.integer(comp.saved), as.integer(maxit),
           as.integer(burn), as.integer(thin), as.integer(adapt), as.integer(batch.size),
           as.double(prop.a), as.double(prop.b),
           as.double(prior.mu), as.double(prior.nu), as.double(prior.eta),
           as.integer(mode), as.integer(submodel),
           as.double(t.a), as.double(t.b), as.double(t.sig),
           as.double(t.mu), as.double(t.w), as.double(t.g),
           as.integer(t.ci), as.integer(t.noo), as.integer(t.noc), as.double(t.sd),
           as.double(start.a), as.double(start.b))
  ## format C output
  ret <- bayesfit()
  ret$a    <- matrix(fit[[18]], nrow=tr.len, ncol=nlag)
  ret$b    <- matrix(fit[[19]], nrow=tr.len, ncol=nlag)
  ret$sd   <- array(fit[[20]], dim = c(tr.len,comp.saved,nlag))
  ret$mean <- array(fit[[21]], dim = c(tr.len,comp.saved,nlag))
  ret$w       <- matrix(fit[[22]], nrow=tr.len, ncol=comp.saved)
  ret$prec    <- fit[[23]]
  ret$ci      <- matrix(fit[[24]], nrow=tr.len, ncol=n)
  ret$noo     <- matrix(fit[[25]], nrow=tr.len, ncol=comp.saved)
  ret$noc     <- fit[[26]]
  ret$prop.sd <- array(fit[[27]], dim=c(tr.len, 8, nlag))
  ret$len <- tr.len
  ret$nlag <- nlag
  ## add names
  ind <- 1:nlag
  colnames(t.a) <- paste("alpha(",ind,")", sep="")
  colnames(t.b) <- paste("beta(",ind,")", sep="")
  #return in R format
  return(ret)
}
