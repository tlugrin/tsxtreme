## Copyright (C) 2017-2025 Thomas Lugrin
## user interface for etfit class
## (Scope:) Fit Eastawn model using Bayesian semiparametrics
## List of functions: - thetafit
##                    - chifit
##################################################

##################################################
## E+T R interface

#' Semiparametric Bayesian estimation of extremal dependence measures
#' 
#' `thetafit` gives posterior samples, mean, median and other chosen quantiles
#' for the extremal index \eqn{\theta(x,m)} and \code{chifit} does the same for
#' the coefficient of extremal dependence \eqn{\chi_m(x)}.
#' Appropriate marginal transforms can be automatically carried out using
#' standard procedures, or performed by the user prior to the extremal
#' dependence model fit. Existing Bayesian fits can be recycled to only estimate
#' a dependence measure \eqn{\theta(x,m)} or \eqn{chi_m(x)}.
#' 
#' @details
#' The sub-asymptotic extremal index is defined as
#' \deqn{\theta(x,m) = Pr(X_1 < x,\ldots,X_m < x | X_0 > x),}
#' whose limit as \eqn{x} and \eqn{m} go to \eqn{\infty}appropriately is the
#' extremal index \eqn{\theta}. The extremal index can be interpreted as the
#' inverse of the asymptotic mean cluster size (see \code{\link{thetaruns})}.
#' 
#' The sub-asymptotic coefficient of extremal dependence is
#' \deqn{\chi_m(x) = Pr(X_m > x | X_0 > x),}
#' whose limit \eqn{\chi} defines asymptotic dependence (\eqn{\chi > 0}) or
#' asymptotic independence (\eqn{\chi = 0}).
#' 
#' Both types of extremal dependence measures can be estimated either using a
#' * **proportion method** (\code{method == "prop"}), sampling from the
#'   conditional probability given \eqn{X_0 > x} and counting the proportion of
#'   sampled points falling in the region of interest, or
#' * **Monte Carlo integration** (\code{method == "MCi"}), sampling replicates
#'   from the marginal exponential tail distribution and evaluating the
#'   conditional tail distribution in these replicates, then taking their mean
#'   as an approximation of the integral.
#'   
#' `submodel == "fom"` imposes a first order Markov structure to the model,
#' namely a geometrical decrease in \eqn{\alpha} and a constant \eqn{\beta}
#' across lags, i.e. \eqn{\alpha_j = \alpha^j} and \eqn{\beta_j = \beta},
#' \eqn{j=1,\ldots,m}.
#' 
#' @param ts a vector, the time series for which to estimate the extremal index
#'   \eqn{\theta(x,m)} or the coefficient of extremal dependence
#'   \eqn{\chi_m(x)}, with \eqn{x} a probability level and \eqn{m} a run-length
#'   (see Details).
#' @param lapl logical; \code{TRUE} indicates that \code{ts} has a marginal
#'   Laplace distribution. If \code{FALSE} (default), \code{method_mar} is used
#'   to transform the marginal distribution of \code{ts} to the Laplace scale.
#' @param nlag the run-length; an integer larger or equal to 1.
#' @param R the number of samples per MCMC iteration drawn from the sampled `S`
#'   posterior distributions; used for the estimation of the dependence measure.
#' @param S the number of posterior distributions sampled from the MCMC trace
#'   to be used for the estimation of the dependence measure.
#' @param u.mar `r lifecycle::badge('deprecated')` use `u_mar` instead.
#' @param u_mar probability; threshold used for marginal transformation if
#'   \code{lapl} is \code{FALSE}. Ignored otherwise.
#' @param u.dep `r lifecycle::badge('deprecated')` use `u_dep` instead.
#' @param u_dep probability; threshold used for the extremal dependence model.
#' @param probs vector of probabilities; the values of \eqn{x} for which to
#'   evaluate \eqn{\theta(x,m)} or \eqn{\chi_m(x)}.
#' @param method.mar `r lifecycle::badge('deprecated')` use `method_mar`
#'   instead.
#' @param method_mar a character string defining the method used to estimate the
#'   marginal GPD; either \code{"mle"} for maximum likelihood of \code{"mom"}
#'   for method of moments or \code{"pwm"} for probability weighted moments
#'   methods. Defaults to \code{"mle"}.
#' @param method a character string defining the method used to estimate the
#'   dependence measure; either \code{"prop"} for proportions or \code{"MCi"}
#'   for Monte Carlo integration (see details).
#' @param fit logical; \code{TRUE} means that the dependence model must be
#'   fitted and the values in \code{par} are used to calibrate the MCMC.
#'   Otherwise `prev_fit` is required and provides the necessary dependence fit.
#' @param prev.fit `r lifecycle::badge('deprecated')` use `prev_fit` instead.
#' @param prev_fit an object of class [bayesfit()], e.g., the result from a
#'   previous call to [depfit()]. Required if \code{fit} is `FALSE`.
#' @param par an object of class [bayesparams()] to be used for the fit of
#'   the dependence model.
#' @param submodel a character string, either \code{"fom"} for
#'   _first order Markov_ or \code{"none"} for no specification.
#' @param levels a vector of probabilites, coverage levels; the quantiles of the
#'   posterior distribution of the extremal measure to be computed.
#' @returns An object of class [depmeasure()] of `type` either "theta" or "chi".
#' @seealso [depfit()], [theta2fit()], [thetaruns()]
#' @examples
#' ## generate data from an AR(1)
#' ## with Gaussian marginal distribution
#' n   <- 10000
#' dep <- 0.5
#' ar    <- numeric(n)
#' ar[1] <- rnorm(1)
#' for (i in 2:n)
#'   ar[i] <- rnorm(1, mean = dep*ar[i-1], sd = 1-dep^2)
#' plot(ar, type = "l")
#' plot(density(ar))
#' grid <- seq(-3,3,0.01)
#' lines(grid, dnorm(grid), col = "blue")
#' 
#' ## rescale the margin (focus on dependence)
#' ar <- qlapl(pnorm(ar))
#' 
#' ## fit the data
#' params <- bayesparams()
#' params$maxit <- 100 # bigger numbers would be
#' params$burn  <- 10  # more sensible...
#' params$thin  <- 4
#' theta <- thetafit(ts = ar, R = 500, S = 100, u_mar = 0.95, u_dep = 0.98,
#'                   probs = c(0.98, 0.999), par = params)
#' ## or, same thing in two steps to control fit output before computing theta:
#' fit <- depfit(ts = ar, u_mar = 0.95, u_dep = 0.98, par = params)
#' plot(fit)
#' theta <- thetafit(ts = ar, R = 500, S = 100, u_mar = 0.95, u_dep = 0.98,
#'                   probs = c(0.98, 0.999), fit = FALSE, prev_fit = fit)
#' @export
thetafit <- function(ts,
                     lapl = FALSE,
                     nlag = 1,
                     R = 1000,
                     S = 500,
                     u_mar = 0,
                     u_dep,
                     probs = seq(u_dep, 0.9999, length.out = 30),
                     method_mar = c("mle","mom","pwm"),
                     method = c("prop","MCi"),
                     fit = TRUE,
                     prev_fit = bayesfit(),
                     par = bayesparams(),
                     submodel = c("fom","none"),
                     levels = c(.025,.975),
                     u.mar = deprecated(),
                     u.dep = deprecated(),
                     method.mar = deprecated(),
                     prev.fit = deprecated()) {
  if (lifecycle::is_present(u.mar)) {
    lifecycle::deprecate_warn("0.4.0", "thetafit(u.mar)", "thetafit(u_mar)")
    u_mar <- u.mar
  }
  if (lifecycle::is_present(u.dep)) {
    lifecycle::deprecate_warn("0.4.0", "thetafit(u.dep)", "thetafit(u_dep)")
    u_dep <- u.dep
  }
  if (lifecycle::is_present(method.mar)) {
    lifecycle::deprecate_warn("0.4.0", "thetafit(method.mar)",
                              "thetafit(method_mar)")
    method_mar <- method.mar
  }
  if (lifecycle::is_present(prev.fit)) {
    lifecycle::deprecate_warn("0.4.0", "thetafit(prev.fit)",
                              "thetafit(prev_fit)")
    prev_fit <- prev.fit
  }
  data_up <- format_ts(ts=ts, u_mar=u_mar, u_dep=u_dep, method=method_mar,
                       nlag=nlag, lapl=lapl)
  ret <- etfit(data=data_up, R=R, S=S, probs=probs, method=method,
              fit=fit, prev_fit=prev_fit, par=par, submodel=submodel, levels=levels)
  mesh_O     <- scale_to_original(p=probs, ts=ts, u=u_mar, gpd_pars=gpd(ts, u=u_mar, method=method_mar)$pars)
  ret$levels <- mesh_O
  return(ret)
}

#' @rdname thetafit
#' @keywords internal
etfit <- function(data, R, S, probs, method,
                  fit, prev_fit, par, submodel, levels) {
  data_up <- data
  n       <- dim(data_up)[1]
  mesh    <- probs
  
  if (fit) {
    if (!is.bayesparams(par)) stop('par must be of class "bayesparams".')
    fit <- htfit(data_up,
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
                 mode = par$mode,
                 submodel = submodel)
    if (par$mode < 2)
      print(paste("Fit of DDP done. Mean value for alpha & beta:",
                  signif (apply(fit$a, 2, mean),3),"and",
                  signif (apply(fit$b, 2, mean),3)))
  } else {
    if (!is.bayesfit(prev_fit)) stop('prev_fit must be of class "bayesfit".')
    if (prev_fit$len == 0)
      stop('prev_fit must contain a proper instance of class "bayesfit".')
    fit <- prev_fit
  }
  # set up parameters for theta
  tictoc::tic()
  nit      <- dim(fit$a)[1]
  nlag     <- dim(fit$a)[2]
  mesh_L   <- qexp(1-(1-mesh)*2)
  nbr_vert <- length(mesh_L)
  sim_U    <- runif (R)
  it <- sample.int(nit, S, replace=TRUE)
  a  <- as.double(fit$a[it,])
  b  <- as.double(fit$b[it,])
  m  <- fit$mean[it,,, drop=FALSE]
  s  <- fit$sd[it,,, drop=FALSE]
  w  <- fit$w[it,, drop=FALSE]
  # prepare returned value
  ret       <- depmeasure("theta")
  ret$fit   <- fit
  ret$probs <- mesh
  ret$nlag  <- nlag
  # method of proportions: generate residuals
  if (grepl("prop",method[1])) {
    sim_Z <- r_res(R=R, S=S, nlag=nlag, w=w, mu=m, sig=s)
  }
  nbr_quant <- length(levels)+2
  th    <- matrix(0, nrow=nbr_vert, ncol=nbr_quant,
                 dimnames=list(NULL,c("mean","median",paste(levels*100,"%",sep=""))))
  distr <- matrix(0, nrow=nbr_vert, ncol=S) # compute coverage/RMSE when true theta is known
  for (i in 1:nbr_vert) {
    sim_L <- -log(2*(1-mesh[i]-(1-mesh[i])*sim_U))
    ## method of proportions
    if (grepl("prop", method[1])) {
      sim_Y <- sim_L%*%t(a) + exp(log(sim_L)%*%t(b))*sim_Z# [RxS*(m-1)]
      sim_Y <- matrix(sim_Y, nrow=R*S, ncol=nlag)
      max_Y  <- apply(sim_Y, 1, max)
      max_Y  <- matrix(max_Y, nrow=R, ncol=S)
      th_tmp <- colSums(max_Y < mesh_L[i])/R
      distr[i,] <- th_tmp
      th[i,]    <- c(mean(th_tmp), median(th_tmp), quantile(th_tmp, levels))
    } else if (grepl("MCi",method[1])) {## Monte Carlo integration
      mesh_Z <- (mesh_L[i]-sim_L%*%t(a))/exp(log(sim_L)%*%t(b))# [RxS(m-1)]
      mesh_Z <- array(mesh_Z, dim=c(R,S,nlag))
      HZ     <- vapply(1:S, p_res, z=mesh_Z, w=w, mu=m, sig=s,
                        FUN.VALUE=numeric(R))# [RxS]
      th_samp <- colMeans(HZ)# [S]
      distr[i,] <- th_samp
      th[i,]    <- c(mean(th_samp), median(th_samp), quantile(th_samp, levels))
    }
  }
  # return and summaries
  ret$theta <- th
  ret$distr <- distr
  time <- toc(silent=TRUE)
  if (par$mode < 2) {
    print(paste("Time elapsed on estimating theta: ",
                time%/%60," min ",round(time-60*(time%/%60), 1)," sec.", sep=""))
    print(paste("Time per x-value: ",round(time/nbr_vert, 1)," sec.", sep=""))
  }
  return(ret)
}



##################################################
## CHI(X) ESTIMATION

#' @rdname thetafit
#' @export
chifit <- function(ts,
                   lapl = FALSE,
                   nlag = 1,
                   R = 1000,
                   S = 500,
                   u_mar = 0,
                   u_dep,
                   probs = seq(u_dep, 0.9999, length.out = 30),
                   method_mar = c("mle","mom","pwm"),
                   method = c("prop","MCi"),
                   fit = TRUE,
                   prev_fit = bayesfit(),
                   par = bayesparams(),
                   submodel = c("fom","none"),
                   levels = c(.025,.975),
                   u.mar = deprecated(),
                   u.dep = deprecated(),
                   method.mar = deprecated(),
                   prev.fit = deprecated()) {
  if (lifecycle::is_present(u.mar)) {
    lifecycle::deprecate_warn("0.4.0", "chifit(u.mar)", "chifit(u_mar)")
    u_mar <- u.mar
  }
  if (lifecycle::is_present(u.dep)) {
    lifecycle::deprecate_warn("0.4.0", "chifit(u.dep)", "chifit(u_dep)")
    u_dep <- u.dep
  }
  if (lifecycle::is_present(method.mar)) {
    lifecycle::deprecate_warn("0.4.0", "chifit(method.mar)",
                              "chifit(method_mar)")
    method_mar <- method.mar
  }
  if (lifecycle::is_present(prev.fit)) {
    lifecycle::deprecate_warn("0.4.0", "chifit(prev.fit)", "chifit(prev_fit)")
    prev_fit <- prev.fit
  }
  # pre-process data
  data_up <- format_ts(ts = ts, u_mar = u_mar, u_dep = u_dep,
                       method = method_mar, nlag = nlag, lapl = lapl)
  n       <- dim(data_up)[1]
  mesh    <- probs
  mesh_O  <- scale_to_original(p = mesh, ts = ts, u = u_mar,
                               gpd_pars = gpd(ts, u = u_mar,
                                              method = method_mar)$pars)
  # H+T fit
  if (fit) {
    if (!is.bayesparams(par)) stop('par must be of class "bayesparams".')
    fit <- htfit(data_up,
                 prop_a = par$prop_a,
                 prop_b = par$prop_b,
                 prior_mu = par$prior_mu,
                 prior_nu = par$prior_nu,
                 prior_eta = par$prior_eta,
                 trunc = par$trunc,
                 comp_saved = par$comp_saved,
                 maxit=par$maxit,
                 burn = par$burn,
                 thin = par$thin,
                 adapt = par$adapt,
                 batch_size = par$batch_size,
                 mode = par$mode,
                 submodel = submodel)
    if (par$mode < 2)
      print(paste("Fit of DDP done. Mean value for alpha & beta:",
                  signif (apply(fit$a, 2, mean),3),"and",
                  signif (apply(fit$b, 2, mean),3)))
  } else {
    if (!is.bayesfit(prev_fit)) stop('prev_fit must be of class "bayesfit".')
    if (prev_fit$len == 0)
      stop('prev_fit must contain a proper instance of class "bayesfit".')
    fit <- prev_fit
  }
  tictoc::tic()
  ret <- depmeasure("chi")
  ret$bayesfit <- fit
  ret$probs    <- mesh
  ret$levels   <- mesh_O
  ret$nlag     <- nlag
  # call etfit on (X_1,X_j) to get distr and/or distr_MC
  len       <- dim(fit$a)[1]
  n_comp    <- dim(fit$mu)[2]
  nbr_vert  <- length(mesh)
  nbr_quant <- length(levels)+2
  chi <- matrix(0, nrow=nbr_vert, ncol=nbr_quant)
  dimnames(chi) <- list(NULL,c("mean","median",paste(levels*100,"%",sep="")))
  subfit <- bayesfit()
  subfit$a <- fit$a[,nlag, drop=FALSE]
  subfit$b <- fit$b[,nlag, drop=FALSE]
  subfit$sd   <- fit$sd[,,nlag, drop=FALSE]
  subfit$mean <- fit$mean[,,nlag, drop=FALSE]
  subfit$w   <- fit$w
  subfit$len <- fit$len
  fit_sample <- etfit(data = data_up[,c(1,nlag+1)],
                      R = R,
                      S = S,
                      probs = mesh,
                      method=method,
                      fit = FALSE,
                      prev_fit = subfit,
                      submodel = submodel,
                      levels = levels)
  # compute chi_j(x) = Pr(X_j>x | X_1>x)
  chi[,"mean"]   <- 1-fit_sample$theta[,"mean"]
  chi[,"median"] <- 1-fit_sample$theta[,"median"]
  if (nbr_quant>2) {
    chi[,3:nbr_quant] <- t(vapply(seq_along(probs), function(i) {
      quantile(1-fit_sample$distr[i,], probs=levels)
    }, levels))
  }
  distr <- 1-fit_sample$distr
  ret$chi   <- chi
  ret$distr <- distr
  ## print time
  time <- toc(silent=TRUE)
  if (par$mode < 2) {
    print(paste("Time elapsed on estimating chi: ",
                time%/%60," min ",round(time-60*(time%/%60), 1)," sec.", sep=""))
    print(paste("Time per x-value: ",round(time/nbr_vert, 1)," sec.", sep=""))
  }
  return(ret)
}
