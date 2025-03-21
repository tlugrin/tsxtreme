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
#' @param u.mar `r lifecycle::badge('deprecated')` use `u_mar` instead.
#' @param u_mar marginal threshold; used when transforming the time series to
#'   the Laplace scale.
#' @param u.dep `r lifecycle::badge('deprecated')` use `u_dep` instead.
#' @param u_dep dependence threshold; level above which the dependence is
#'   modelled. `u_dep` can be lower than `u_mar`.
#' @param lapl logical; is \code{ts} on the Laplace scale already? The default
#'   (FALSE) assumes unknown marginal distribution.
#' @param method.mar `r lifecycle::badge('deprecated')` use `method_mar`
#'   instead.
#' @param method_mar a character string defining the method used to estimate the
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
#' fit <- depfit(ts=ar, u_mar=0.95, u_dep=0.98, par=params)
#'
#' ########
#' ## density estimation with submodel=="ugm"
#' data <- MASS::galaxies/1e3
#' dens <- depfit(ts=data, par=params, submodel="ugm")
#' @export
depfit <- function(ts,
                   u_mar = 0,
                   u_dep = u_mar,
                   lapl = FALSE,
                   method_mar = c("mle","mom","pwm"),
                   nlag = 1,
                   par = bayesparams(),
                   submodel = c("fom","none","ugm"),
                   u.mar = deprecated(),
                   u.dep = deprecated(),
                   method.mar = deprecated()) {
  if (lifecycle::is_present(u.mar)) {
    lifecycle::deprecate_warn("0.4.0", "depfit(u.mar)", "depfit(u_mar)")
    u_mar <- u.mar
  }
  if (lifecycle::is_present(u.dep)) {
    lifecycle::deprecate_warn("0.4.0", "depfit(u.dep)", "depfit(u_dep)")
    u_dep <- u.dep
  }
  if (lifecycle::is_present(method.mar)) {
    lifecycle::deprecate_warn("0.4.0", "depfit(method.mar)",
                              "depfit(method_mar)")
    method_mar <- method.mar
  }
  if (submodel[1]=="ugm") {
    data_up <- cbind(1,ts)
  } else if (submodel[1] %in% c("fom","none")) {
    data_up <- format_ts(ts = ts, u_mar = u_mar, u_dep = u_dep,
                         method = method_mar, nlag = nlag, lapl = lapl)
  } else {
    stop(paste("submodel", submodel, "unknown in depfit()"))
  }
  htfit(data = data_up,
        prop_a = par$prop_a,
        prop_b = par$prop_b,
        prior_mu = par$prior_mu,
        prior_nu = par$prior_nu,
        prior_eta = par$prior_eta,
        trunc = par$trunc,
        comp_saved = par$comp_saved,
        maxit = par$maxit,
        burn = par$burn,
        thin = par$thin,
        adapt = par$adapt,
        batch_size = par$batch_size,
        start_ab = par$start_ab,
        conditions = par$conditions,
        mode = par$mode,
        submodel = submodel)
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
                  prop_a,
                  prop_b,
                  prior_mu = c(0,10),
                  prior_nu = c(2,1/2),
                  prior_eta = c(4,1),
                  trunc = 100,
                  comp_saved = 10,
                  maxit = 10000,
                  burn = 2000,
                  thin = 4,
                  adapt = 2000,
                  batch_size = 125,
                  start_ab = c("guesstimate", "prior"),
                  conditions = TRUE,
                  mode = 1,
                  submodel = c("fom","none","ugm")) {
  n    <- dim(data)[1]
  nlag <- dim(data)[2]-1
  if (length(prop_a) < 5) { prop_a <- rep(prop_a[1], 5) }
  if (length(prop_b) < 3) { prop_b <- rep(prop_b[1], 3) }
  ## build trace containers
  tr_len <- floor((maxit-burn)/thin)
  t_a   <- t_b <- matrix(0, nrow = tr_len, ncol = nlag)
  t_sig <- t_mu <- array(0, dim = c(tr_len, comp_saved, nlag))
  t_w   <- matrix(0, nrow = tr_len, ncol = comp_saved)
  t_g   <- numeric(tr_len)
  t_ci  <- matrix(0, nrow = tr_len, ncol = n)
  t_noo <- matrix(0, nrow = tr_len, ncol = comp_saved)
  t_noc <- numeric(tr_len)
  t_sd  <- array(0, dim = c(tr_len, 8, nlag))#for RAMA
  if (submodel[1] == "ugm") {
    start_a <- start_b <- 0
  } else {
    if (start_ab[1] == "guesstimate") {
      fit2step <- ht2step(data, conditions) # lags > 1 will be properly
      start_a <- fit2step$a                 # initialised in C routine
      start_b <- fit2step$b                 #  depending on [submodel]
      if (mode == 0) {
        message("DEBUG: Guesstimates for alpha and beta: ", start_a,
                " and ", start_b)
      }
    } else if (start_ab[1] == "prior") {
      start_a <- start_b <- rep(Inf, nlag)# will be properly initialised later
    } else {
      stop("In htfit(): invalid starting value strategy for start_ab")
    }
  }
  submodel <- switch(submodel[1], "ugm" = 0, "fom" = 1, "none" = 2,
                     stop("In htfit(): invalid submodel."))
  ## call C(++) wrapper
  fit = .C(C_et_interface,
           as.double(data),
           as.integer(n),
           as.integer(nlag),
           as.integer(trunc),
           as.integer(comp_saved),
           as.integer(maxit),
           as.integer(burn),
           as.integer(thin),
           as.integer(adapt),
           as.integer(batch_size),
           as.double(prop_a),
           as.double(prop_b),
           as.double(prior_mu),
           as.double(prior_nu),
           as.double(prior_eta),
           as.integer(conditions),
           as.integer(mode),
           as.integer(submodel),
           as.double(t_a),
           as.double(t_b),
           as.double(t_sig),
           as.double(t_mu),
           as.double(t_w),
           as.double(t_g),
           as.integer(t_ci),
           as.integer(t_noo),
           as.integer(t_noc),
           as.double(t_sd),
           as.double(start_a),
           as.double(start_b))
  ## format C output
  ret <- bayesfit()
  ret$a    <- matrix(fit[[19]], nrow = tr_len, ncol = nlag)
  ret$b    <- matrix(fit[[20]], nrow = tr_len, ncol = nlag)
  ret$sd   <- array(fit[[21]], dim = c(tr_len, comp_saved, nlag))
  ret$mean <- array(fit[[22]], dim = c(tr_len, comp_saved, nlag))
  ret$w       <- matrix(fit[[23]], nrow = tr_len, ncol = comp_saved)
  ret$prec    <- fit[[24]]
  ret$ci      <- matrix(fit[[25]], nrow = tr_len, ncol = n)
  ret$noo     <- matrix(fit[[26]], nrow = tr_len, ncol = comp_saved)
  ret$noc     <- fit[[27]]
  ret$prop_sd <- array(fit[[28]], dim=c(tr_len, 8, nlag))
  ret$len <- tr_len
  ret$nlag <- nlag
  ## add names
  ind <- 1:nlag
  colnames(t_a) <- paste("alpha(",ind,")", sep = "")
  colnames(t_b) <- paste("beta(",ind,")", sep = "")
  #return in R format
  return(ret)
}
