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
#'   Laplace distribution. If \code{FALSE} (default), \code{method.mar} is used
#'   to transform the marginal distribution of \code{ts} to the Laplace scale.
#' @param nlag the run-length; an integer larger or equal to 1.
#' @param R the number of samples per MCMC iteration drawn from the sampled `S`
#'   posterior distributions; used for the estimation of the dependence measure.
#' @param S the number of posterior distributions sampled from the MCMC trace
#'   to be used for the estimation of the dependence measure.
#' @param u.mar probability; threshold used for marginal transformation if
#'   \code{lapl} is \code{FALSE}. Ignored otherwise.
#' @param u.dep probability; threshold used for the extremal dependence model.
#' @param probs vector of probabilities; the values of \eqn{x} for which to
#'   evaluate \eqn{\theta(x,m)} or \eqn{\chi_m(x)}.
#' @param method.mar a character string defining the method used to estimate the
#'   marginal GPD; either \code{"mle"} for maximum likelihood of \code{"mom"}
#'   for method of moments or \code{"pwm"} for probability weighted moments
#'   methods. Defaults to \code{"mle"}.
#' @param method a character string defining the method used to estimate the
#'   dependence measure; either \code{"prop"} for proportions or \code{"MCi"}
#'   for Monte Carlo integration (see details).
#' @param silent logical (\code{FALSE}); verbosity.
#' @param fit logical; \code{TRUE} means that the dependence model must be
#'   fitted and the values in \code{par} are used to calibrate the MCMC.
#'   Otherwise `prev.fit` is required and provides the necessary dependence fit.
#' @param prev.fit an object of class [bayesfit()], e.g., the result from a
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
#' for(i in 2:n)
#'   ar[i] <- rnorm(1, mean=dep*ar[i-1], sd=1-dep^2)
#' plot(ar, type="l")
#' plot(density(ar))
#' grid <- seq(-3,3,0.01)
#' lines(grid, dnorm(grid), col="blue")
#' 
#' ## rescale the margin (focus on dependence)
#' ar <- qlapl(pnorm(ar))
#' 
#' ## fit the data
#' params <- bayesparams()
#' params$maxit <- 100 # bigger numbers would be
#' params$burn  <- 10  # more sensible...
#' params$thin  <- 4
#' theta <- thetafit(ts=ar, R=500, S=100, u.mar=0.95, u.dep=0.98,
#'                   probs = c(0.98, 0.999), par=params)
#' ## or, same thing in two steps to control fit output before computing theta:
#' fit <- depfit(ts=ar, u.mar=0.95, u.dep=0.98, par=params)
#' plot(fit)
#' theta <- thetafit(ts=ar, R=500, S=100, u.mar=0.95, u.dep=0.98,
#'                   probs = c(0.98, 0.999), fit=FALSE, prev.fit=fit)
#' @export
thetafit <- function(ts, lapl=FALSE, nlag=1, R=1000, S=500,
                  u.mar=0, u.dep, probs=seq(u.dep,0.9999,length.out=30),
                  method.mar=c("mle","mom","pwm"), method=c("prop","MCi"),
                  silent=FALSE,
                  fit=TRUE, prev.fit=bayesfit(), par=bayesparams(),
                  submodel=c("fom","none"), levels=c(.025,.975)){
  data.up <- format.ts(ts=ts, u.mar=u.mar, u.dep=u.dep, method=method.mar,
                       nlag=nlag, lapl=lapl)
  ret <- etfit(data=data.up, R=R, S=S, probs=probs, method=method, silent=silent,
              fit=fit, prev.fit=prev.fit, par=par, submodel=submodel, levels=levels)
  mesh.O     <- scale.to.original(p=probs, ts=ts, u=u.mar, gpd.pars=gpd(ts, u=u.mar, method=method.mar)$pars)
  ret$levels <- mesh.O
  return(ret)
}

#' @rdname thetafit
#' @keywords internal
etfit <- function(data, R, S, probs, method,
                  silent,
                  fit, prev.fit, par, submodel, levels){
  data.up <- data
  n       <- dim(data.up)[1]
  mesh    <- probs
  
  if(fit){
    if(!is.bayesparams(par)) stop("par must be of class 'bayesparams'.")
    fit <- htfit(data.up, prop.a=par$prop.a, prop.b=par$prop.b,
                prior.mu=par$prior.mu, prior.nu=par$prior.nu, prior.eta=par$prior.eta,
                trunc=par$trunc, comp.saved=par$comp.saved, maxit=par$maxit,
                burn=par$burn, thin=par$thin,
                adapt=par$adapt, batch.size=par$batch.size, mode=par$mode, submodel=submodel)
    if(!silent)
      print(paste("Fit of DDP done. Mean value for alpha & beta:",
                  signif(apply(fit$a, 2, mean),3),"and",
                  signif(apply(fit$b, 2, mean),3)))
  }else{
    if(!is.bayesfit(prev.fit)) stop("prev.fit must be of class 'bayesfit'.")
    if(prev.fit$len == 0) stop("prev.fit must contain a proper instance of class 'bayesfit'.")
    fit <- prev.fit
  }
  # set up parameters for theta
  tictoc::tic()
  nit      <- dim(fit$a)[1]
  nlag     <- dim(fit$a)[2]
  mesh.L   <- qexp(1-(1-mesh)*2)
  nbr.vert <- length(mesh.L)
  sim.U    <- runif(R)
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
  if(grepl("prop",method[1])){
    sim.Z <- r.res(R=R, S=S, nlag=nlag, w=w, mu=m, sig=s)
  }
  nbr.quant <- length(levels)+2
  th    <- matrix(0, nrow=nbr.vert, ncol=nbr.quant,
                 dimnames=list(NULL,c("mean","median",paste(levels*100,"%",sep=""))))
  distr <- matrix(0, nrow=nbr.vert, ncol=S) # compute coverage/RMSE when true theta is known
  for(i in 1:nbr.vert){
    sim.L <- -log(2*(1-mesh[i]-(1-mesh[i])*sim.U))
    ## method of proportions
    if(grepl("prop", method[1])){
      sim.Y <- sim.L%*%t(a) + exp(log(sim.L)%*%t(b))*sim.Z# [RxS*(m-1)]
      sim.Y <- matrix(sim.Y, nrow=R*S, ncol=nlag)
      max.Y  <- apply(sim.Y, 1, max)
      max.Y  <- matrix(max.Y, nrow=R, ncol=S)
      th.tmp <- colSums(max.Y < mesh.L[i])/R
      distr[i,] <- th.tmp
      th[i,]    <- c(mean(th.tmp), median(th.tmp), quantile(th.tmp, levels))
    }## Monte Carlo integration
    else if(grepl("MCi",method[1])){
      mesh.Z <- (mesh.L[i]-sim.L%*%t(a))/exp(log(sim.L)%*%t(b))# [RxS(m-1)]
      mesh.Z <- array(mesh.Z, dim=c(R,S,nlag))
      HZ      <- vapply(1:S, p.res, z=mesh.Z, w=w, mu=m, sig=s, FUN.VALUE=numeric(R))# [RxS]
      th.samp <- colMeans(HZ)# [S]
      distr[i,] <- th.samp
      th[i,]    <- c(mean(th.samp), median(th.samp), quantile(th.samp, levels))
    }
  }
  # return and summaries
  ret$theta <- th
  ret$distr <- distr
  time <- toc(silent=TRUE)
  if(!silent){
    print(paste("Time elapsed on estimating theta: ",
                time%/%60," min ",round(time-60*(time%/%60), 1)," sec.", sep=""))
    print(paste("Time per x-value: ",round(time/nbr.vert, 1)," sec.", sep=""))
  }
  return(ret)
}



##################################################
## CHI(X) ESTIMATION

#' @rdname thetafit
#' @export
chifit <- function(ts, lapl=FALSE, nlag=1, R=1000, S=500,
                   u.mar=0, u.dep, probs=seq(u.dep,0.9999,length.out=30),
                   method.mar=c("mle","mom","pwm"), method=c("prop","MCi"),
                   silent=FALSE,
                   fit=TRUE, prev.fit=bayesfit(), par=bayesparams(),
                   submodel=c("fom","none"), levels=c(.025,.975)){
  # pre-process data
  data.up <- format.ts(ts=ts, u.mar=u.mar, u.dep=u.dep, method=method.mar,
                       nlag=nlag, lapl=lapl)
  n       <- dim(data.up)[1]
  mesh    <- probs
  mesh.O  <- scale.to.original(p=mesh, ts=ts, u=u.mar, gpd.pars=gpd(ts, u=u.mar, method=method.mar)$pars)
  # H+T fit
  if(fit){
    if(!is.bayesparams(par)) stop("par must be of class 'bayesparams'.")
    fit <- htfit(data.up, prop.a=par$prop.a, prop.b=par$prop.b,
                prior.mu=par$prior.mu, prior.nu=par$prior.nu, prior.eta=par$prior.eta,
                trunc=par$trunc, comp.saved=par$comp.saved, maxit=par$maxit,
                burn=par$burn, thin=par$thin,
                adapt=par$adapt, batch.size=par$batch.size, mode=par$mode, submodel=submodel)
    if(!silent)
      print(paste("Fit of DDP done. Mean value for alpha & beta:",
                  signif(apply(fit$a, 2, mean),3),"and",
                  signif(apply(fit$b, 2, mean),3)))
  }else{
    if(!is.bayesfit(prev.fit)) stop("prev.fit must be of class 'bayesfit'.")
    if(prev.fit$len == 0) stop("prev.fit must contain a proper instance of class 'bayesfit'.")
    fit <- prev.fit
  }
  tictoc::tic()
  ret <- depmeasure("chi")
  ret$bayesfit <- fit
  ret$probs    <- mesh
  ret$levels   <- mesh.O
  ret$nlag     <- nlag
  # call etfit on (X_1,X_j) to get distr and/or distr.MC
  len       <- dim(fit$a)[1]
  n.comp    <- dim(fit$mu)[2]
  nbr.vert  <- length(mesh)
  nbr.quant <- length(levels)+2
  chi <- matrix(0, nrow=nbr.vert, ncol=nbr.quant)
  dimnames(chi) <- list(NULL,c("mean","median",paste(levels*100,"%",sep="")))
  subfit <- bayesfit()
  subfit$a <- fit$a[,nlag, drop=FALSE]
  subfit$b <- fit$b[,nlag, drop=FALSE]
  subfit$sd   <- fit$sd[,,nlag, drop=FALSE]
  subfit$mean <- fit$mean[,,nlag, drop=FALSE]
  subfit$w   <- fit$w
  subfit$len <- fit$len
  fit.sample <- etfit(data=data.up[,c(1,nlag+1)], R=R, S=S, probs=mesh, method=method,
                      silent=silent,
                      fit=FALSE, prev.fit=subfit, submodel=submodel, levels=levels)
  # compute chi_j(x) = Pr(X_j>x | X_1>x)
  chi[,"mean"]   <- 1-fit.sample$theta[,"mean"]
  chi[,"median"] <- 1-fit.sample$theta[,"median"]
  if(nbr.quant>2){
    chi[,3:nbr.quant] <- t(vapply(seq_along(probs), function(i){
      quantile(1-fit.sample$distr[i,], probs=levels)
    }, levels))
  }
  distr <- 1-fit.sample$distr
  ret$chi   <- chi
  ret$distr <- distr
  ## print time
  time <- toc(silent=TRUE)
  if(!silent){
    print(paste("Time elapsed on estimating chi: ",
                time%/%60," min ",round(time-60*(time%/%60), 1)," sec.", sep=""))
    print(paste("Time per x-value: ",round(time/nbr.vert, 1)," sec.", sep=""))
  }
  return(ret)
}