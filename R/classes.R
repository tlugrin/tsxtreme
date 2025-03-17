## Copyright (C) 2017-2025 Thomas Lugrin
## Definition of methods and related functions
## (Scope:) single framework for 1- and 2-stage functions
## (Note:) summary is tricky especially when multi-dimensional!
## List of functions: - stepfit
##                    - is.stepfit
##                    - print/summary.stepfit
##                    - plot.stepfit
##                    - bayesfit
##                    - is.bayesfit
##                    - print/summary.bayesfit
##                    - plot.bayesfit
##                    - bayesparams
##                    - is.bayesparams
##                    - print/summary.bayesparams
##                    - depmeasure
##                    - is.depmeasure
##                    - print/summary.depmeasure
##                    - plot.depmeasure
##                    - is_int
##################################################

##################################################
## CLASS "STEPFIT"
#' Estimates from the stepwise fit
#' 
#' Create, test or show objects of class "stepfit".
#' 
#' @param x an arbitrary \R object.
#' @returns An object of class "stepfit" with the following elements:
#'   \item{a, b }{Heffernan--Tawn model parameters of length `nlag`; \eqn{\alpha} controls the conditional extremal expectation, while \eqn{\beta} controls the conditional extremal expectation and variance.}
#'   \item{res }{matrix of fitted residuals with `nlag` columns.}
#'   \item{pars.se} `r lifecycle::badge('deprecated')` use `pars_se` instead.
#'   \item{pars_se }{2-column matrix of estimated standard errors for \code{a} (first column) and \code{b} (second column), given by the hessian matrix of the likelihood function used in the first step of the inference procedure.}
#'   \item{nlag }{number of lags.}
#' @seealso [bayesfit()], [depmeasure()]
#' @export
stepfit <- function() {
  x <- list(a=numeric(0), b=numeric(0), res=matrix(0,nrow=1,ncol=1), pars_se=matrix(0,nrow=1,ncol=2), nlag=0)
  class(x) <- "stepfit"
  return(x)
}

# validator
#' @export
is.stepfit <- function(x) {
  if (length(names(x))) {
    if ("pars.se" %in% names(x)) {
      lifecycle::deprecate_warn(when = "0.4.0",
                                what = I('The "stepfit" element "pars.se"'),
                                with = I('"pars_se"'))
      x$pars_se <- x$pars.se
    }
    conds <- names(stepfit()) %in% names(x)# minimal requirement
    conds <- sum(conds) == length(names(stepfit()))
    if (conds) {
      conds <- is.numeric(x$a) && is.numeric(x$b) && is.matrix(x$res) && is.matrix(x$pars_se)
      if (!conds) return(FALSE)
      conds <- dim(x$pars_se)[1]==2 && length(x$a)==length(x$b) && dim(x$pars_se)[2]==length(x$a) && length(x$nlag) == 1
      return(inherits(x, "stepfit") && conds)
    }
  }
  return(FALSE)
}

## viewer
#' @export
print.stepfit <- function(x, ...) {
  summary.stepfit(x, ...)
}

#' @export
summary.stepfit <- function(object, ...) {
  if (!is.stepfit(object)) {
    warning('"object" does not comply with the "stepfit" class')
    summary.default(object, ...)
  }
  cat(" Parameters:\n")
  for (ind in 1:object$nlag) {
    cat(paste("alpha(",ind,")  ",signif (object$a,3)," (",signif (object$pars_se[1,ind],3),")\n", sep=""), sep="")
    cat(paste("beta(",ind,")   ",signif (object$b,3)," (",signif (object$pars_se[2,ind],3),")\n", sep=""), sep="")
  }
  cat(" Residuals:\n")
  print(summary(object$res))
}

## plotter
#' @rdname stepfit
#' @export
plot.stepfit <- function(x, ...) {
  stopifnot(is.stepfit(x))
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  par(ask=TRUE)
  def <- list(xlab="z", ylab=bquote(h(z)), main="")
  if (length(list(...)) > 0) {
    pm <- !pmatch(names(def), names(list(...)), 0)# all non-matches
    def <- def[pm]
  }
  for (j in 1:x$nlag) {
    if ("ylab" %in% names(def)) def$ylab <- bquote(bquote(h(z[.(j)])))
    if ("xlab" %in% names(def)) def$xlab <- bquote(bquote(z[.(j)]))
    den <- density(x$res[,j])
    def$x <- x$res[,j]
    do.call(hist, c(def, prob=TRUE, list(...)))
    lines(den, lty="dashed")
  }
}


##################################################
## CLASS "BAYESFIT"
#' Traces from the MCMC output
#' 
#' Test or show objects of class "bayesfit".
#' 
#' @details
#' Default plot shows samples of residual densities (\code{which==1}), residual
#' distribution with credible interval (\eqn{5\%} and \eqn{95\%} posterior
#' quantile; \code{which==2}), and joint posterior distribution of \eqn{\alpha}
#' and \eqn{\beta} (\code{which==3}) for each lag successively. \code{which}
#' can be any composition of 1,2 and 3.
#' 
#' @param x an arbitrary \R object.
#' @param which a vector with values in {1,2,3} where 1 is to plot residual
#'   density functions, 2 is for residual distribution functions and 3 is for
#'   a contour plot of the posterior distribution function of (alpha, beta).
#' 
#' @returns An object of class "bayesfit" with MCMC traces for:
#'   \item{a, b }{Heffernan--Tawn model parameters (matrices).}
#'   \item{sd, mean }{means and standard deviations of the residual mixture components (3-dimensional arrays).}
#'   \item{w }{weights of the mixture components (matrix).}
#'   \item{prec }{precision parameter of the Dirichlet process (vector).}
#'   \item{ci }{auxiliary variable; components' indices for each observation (matrix).}
#'   \item{noo }{number of observation in each mixture component (matrix).
#'   \item{noc }{number of non-empty components in the mixture, i.e., for which at least one data points has a component index pointing to it.}
#'   \item{prop.sd } `r lifecycle::badge('deprecated')` use `prop_sd` instead.
#'   \item{prop_sd }{standard deviations for the proposal distributions of \code{a} and \code{b} (3-dimensional array).}
#'   And \code{len}, the length of the traces, i.e., the number of iterations
#'   saved after burning and thinning; \code{nlag} stores the number of lags
#'   present in the data.
#' @seealso [bayesparams()], [stepfit()]
#' @export
bayesfit <- function() {
  z_mat <- matrix(0, nrow=1, ncol=1)
  z_ar  <- array(0, dim=c(1,1,1))
  x <- list(a=z_mat, b=z_mat, sd=z_ar, mean=z_ar, w=z_mat, prec=0,
              ci=z_mat, noo=z_mat, noc=0, prop_sd=z_ar,
              len=0, nlag=0)
  class(x) <- "bayesfit"
  return(x)
}

## validator
#' @export
is.bayesfit <- function(x) {
  if (length(names(x))) {
    if ("prop.sd" %in% names(x)) {
      lifecycle::deprecate_warn(when = "0.4.0",
                                what = I('The "bayesfit" element "prop.sd"'),
                                with = I('"prop_sd"'))
      x$prop_sd <- x$prop.sd
    }
    conds <- names(bayesfit()) %in% names(x)
    conds <- sum(conds) == length(names(bayesfit()))
    if (all(conds)) {
      conds <- is.matrix(x$a) && is.matrix(x$b) &&
        is.array(x$sd) && is.array(x$mean) &&
        is.matrix(x$w) && is.numeric(x$prec) && is.array(x$prop_sd) &&
        is.matrix(x$ci) && is.matrix(x$noo) && is.numeric(x$noc) && is.numeric(x$len) && is.numeric(x$nlag) &&
        is_int(x$len) && length(x$len)==1 && is_int(x$nlag) && length(x$nlag)==1 &&
        is_int(x$noc) && is_int(x$noo) && is_int(x$ci)
      return(inherits(x, "bayesfit") && conds)
    }
  }
  return(FALSE)
}

## viewer
#' @export
print.bayesfit <- function(x, ...) {
  summary.bayesfit(x, ...)
}

#' @export
summary.bayesfit <- function(object, ...) {
  if (!is.bayesfit(object)) {
    warning('"object" does not comply with the "bayesfit" class')
    summary.default(object, ...)
  }
  cat(" Posterior medians:\n")
  for (j in 1:object$nlag) {
    cat("alpha(",j,")  ",median(object$a[,j]),"\n", sep="")
    cat("beta(",j,")   ",median(object$b[,j]),"\n", sep="")
    cat("mean(",j,",",1,") ",median(object$mean[,1,j]),"\n", sep="")
    cat("sd(",j,",",1,")   ",median(object$sd[,1,j]),"\n", sep="")
  }
}

#' @rdname bayesfit
#' @export
plot.bayesfit <- function(x, which=1:3, ...) {
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  stopifnot(is.bayesfit(x))
  stopifnot(is_int(which))
  stopifnot(min(which) >= 1 && max(which) <= 3)
  which <- sort(unique(which))
  ugm <- all(x$a==0) && all(x$b==0)
  if (ugm) which <- intersect(which,1:2)
  if (x$nlag > 1) par(ask=TRUE)
  if (length(which) == 2) par(mfrow=c(1,2))
  else if (length(which) > 2) par(mfrow=c(2,2))
  def <- list(xlab="x", ylab="y", main="")
  if (length(list(...)) > 0) {
    pm <- !pmatch(names(def), names(list(...)), 0)# all non-matches
    def <- def[pm]
  }
  for (j in 1:x$nlag) {
    if (1 %in% which) {
      if ("xlab" %in% names(def)) def$xlab <- bquote(bquote(z[.(j)]))
      if ("ylab" %in% names(def)) def$ylab <- bquote(bquote(h(z[.(j)])))
      if ("main" %in% names(def)) def$main <- "Sample of residual densities"
        grid <- seq(-100,100,2)
        de <- residual_densities(x,j,grid)
        grid <- seq(min(de$x),max(de$x),length.out=501)
        de <- residual_densities(x,j,grid)
        xy <- list(x=de$x, y=de$y[,1])
        if (!("ylim" %in% names(list(...)))) xy$ylim <- c(0,max(de$y))
        do.call(plot, c(xy, type="l",def,list(...)))
        for (it in 1:dim(de$y)[2])
          lines(de$x, de$y[,it])
    }
    if (2 %in% which) {
      if ("xlab" %in% names(def)) def$xlab <- bquote(bquote(z[.(j)]))
      if ("ylab" %in% names(def)) def$ylab <- bquote(bquote(H(z[.(j)])))
      if ("main" %in% names(def)) def$main <- "Residual distribution"
      grid <- seq(-100,100,2)
      ds <- residual_distributions(x,j,grid)
      grid <- seq(min(ds$x),max(ds$x),length.out=501)
      ds <- residual_distributions(x,j,grid)
      xy <- list(x=ds$x, y=ds$y[,1])
      do.call(plot, c(xy, type="l",def,list(...)))
      for (q in 2:3)
        lines(ds$x, ds$y[,q], lty="dashed")
    }
    if (3 %in% which) {
      if ("xlab" %in% names(def)) def$xlab <- bquote(bquote(alpha[.(j)]))
      if ("ylab" %in% names(def)) def$ylab <- bquote(bquote(beta[.(j)]))
      if ("main" %in% names(def)) def$main <- "Joint posterior of H-T par."
      lims <- c(max(min(x$a[,j]),-1),min(max(x$a[,j]),1),max(min(x$b[,j]),0),min(max(x$b[,j]),1))
      kd <- kde2d(x$a[,j], x$b[,j], lims=lims)
      do.call(contour, c(kd,nlevels=5,def,list(...)))
    }
  }
}

#' Compute the residual density/distribution function on a grid of x-values
#' 
#' Helper function for plotting the mixture density or distribution function
#' of the residuals as expressed in the conditional tail model. Given a broad
#' range of x-values, trims the left and right tails to cover the centre of the
#' distribution.
#' 
#' @details
#' `r lifecycle::badge('experimental')` Initially intended for internal use
#' only, but exposed to the user as it may prove useful to adjust plots to
#' different situations. Only minimal checks are performed on the inputs.
#' You have been warned! 
#' 
#' @param x a bayesfit object
#' @param lag the dimension of the residual density/distribution function to
#'   consider; starts at 1
#' @param grid vector of x-values to evaluate the function at
#' @returns A list with two elements, x and y, where x is a subset of `grid`.
#' @export
residual_densities <- function(x, lag, grid) {
  stopifnot(is.bayesfit(x))
  nsamp<- min(25,x$len)# nbr of samples
  noc  <- dim(x$mean)[2]
  nobs <- sum(x$noo[1,])
  ind <- sample.int(x$len, nsamp, replace=FALSE)
  dens <- vapply(ind, function(it,j) {
    de <- vapply(1:noc, function(c,j,it) {
        x$noo[it,c]*dnorm(grid, x$mean[it,c,j], x$sd[it,c,j])
      }, grid, j=j, it=it)
      rowSums(de)# sample density
    }, grid, j=lag)
  dens <- dens/nobs
  maxs <- vapply(1:nsamp, function(j) max(dens[,j]), 0)
  wmin <- which.min(maxs)# expected to be the flattest
  minimax <- min(maxs)
  fact <- 1e3# trim tails which are at least [fact] times smaller than minimax
  lo <- 1; up <- length(grid)
  i <- 1
  while(i < up && dens[i,wmin] < minimax/fact)
    i <- i+1
  lo <- max(lo,i-3); i <- up
  while(i > lo && dens[i,wmin] < minimax/fact)
    i <- i-1
  up <- min(up,i+3)
  return(list(x=grid[lo:up], y=dens[lo:up,]))
}

#' @rdname residual_densities
#' @export
residual_distributions <- function(x, lag, grid) {
  stopifnot(is.bayesfit(x))
  nsamp<- min(500,x$len)# nbr of samples
  noc  <- dim(x$mean)[2]
  nobs <- sum(x$noo[1,])
  ind <- sample.int(x$len, nsamp, replace=FALSE)
  dstr <- vapply(ind, function(it, j) {
    ds <- vapply(1:noc, function(c,j,it) {
      x$noo[it,c]*pnorm(grid, x$mean[it,c,j], x$sd[it,c,j])
    }, grid, j=j, it=it)
    rowSums(ds)# sample distribution
  }, grid, j=lag)
  dstr <- dstr/nobs
  dstr <- vapply(seq_along(grid), function(gr) {
    quantile(dstr[gr,], probs = c(.5,.05,.95))
  }, numeric(3))
  dstr <- t(dstr)
  # trim tails of distribution
  lo <- 1; up <- length(grid)
  if (grid[2] - grid[1] > 1) {
    eps_prob <- 1e-3
    i <- 1
    while(i < up && dstr[i,1] < eps_prob)
      i <- i+1
    lo <- max(lo,i-1); i <- up
    while(i > lo && dstr[i,1] > 1 - eps_prob)
      i <- i - 1
    up <- min(up, i + 1)
  }
  return(list(x = grid[lo:up], y = dstr[lo:up,]))
}


##################################################
## CLASS "BAYESPARAMS"
#' Parameters for the Bayesian semi-parametric approach
#' 
#' Create, test or show objects of class "bayesparams". Objects of this class
#' are used as a meta-parameter for the methods fitting the Bayesian approach.
#' 
#' @details 
#' \code{prop_a} is a vector of length 5 with the standard deviations for each
#' region of the RAMA for the (Gaussian) proposal for \eqn{\alpha}. If a scalar
#' is given, 5 identical values are assumed.
#' 
#' \code{prop_b} is a vector of length 3 with the standard deviations for each
#' region of the RAMA for the (Gaussian) proposal for \eqn{\beta}. If a scalar
#' is provided, 3 identical values are assumed.
#' 
#' \code{comp_saved} has no impact on the calculations: its only purpose is to
#' prevent from storing huge amounts of empty components.
#' 
#' The regional adaption scheme targets a \eqn{0.44} acceptance probability. It
#' splits \eqn{[-1;1]} in \eqn{5} regions for \eqn{\alpha} and \eqn{[0;1]} in
#' \eqn{3} regions for \eqn{\beta}. The decision to increase/decrease the
#' proposal standard deviation is based on the past \code{batch_size} MCMC
#' iterations, so too low values yield inefficient adaption, while too large
#' values yield slow adaption.
#' 
#' Default values for the hyperparameters are chosen to get reasonably
#' uninformative priors.
#' 
#' `r lifecycle::badge('deprecated')` The **dot-notation** for elements of
#' the "bayesparams" object is deprecated in favour of the **snake-notation**.
#' For example, `prop.a` is now `prop_a`, `prior.mu` is now `prior_mu`,
#' and so on. The dot-notation elements will be removed in a subsequent version
#' of the package.
#' 
#' @param x an arbitrary \R object.
#' @param prop_a,prop_b standard deviation for the Gaussian proposal
#'   distributions of the Heffernan--Tawn model parameters; see also Details.
#' @param prior_mu mean and standard deviation parameters for the Gaussian prior
#'   distribution of the components' means.
#' @param prior_nu shape and rate parameters for the inverse gamma prior
#'   distribution of the components' variances.
#' @param prior_eta shape and scale parameters for the gamma prior distribution
#'   of the precision parameter used in the Dirichlet process.
#' @param trunc integer; truncation parameter for approximating the infinite sum
#'   in the stick-breaking process.
#' @param comp_saved number of components to be stored and returned.
#' @param maxit total number of iterations to perform (fewer are saved based on
#'   \code{burn} and \code{thin}).
#' @param burn number of first iterations to discard.
#' @param thin spacing between iterations to be saved. Defaults to 1, i.e., all
#'   iterations are saved.
#' @param adapt number of iterations during which an adaptive scheme (RAMA:
#'   Regional Adaptive Monte-Carlo Algorithm) is applied to the proposal
#'   variances of \eqn{\alpha} and \eqn{\beta}; 0 means no adaption.
#' @param batch_size size of batches used in the adaption algorithm. It has no effect if \code{adapt==0}.
#' @param start_ab either "guesstimate" to initialise start values for alpha and
#'   beta from a stepwise fit or "prior" to use random draws from their
#'   respective prior distributions.
#' @param mode verbosity; 0 for debug mode, 1 (default) for standard output,
#'   and 2 for silent.
#' @returns An object of class "bayesparams".
#' @seealso [bayesfit()], [depmeasure()]
#' @examples
#' is.bayesparams(bayesparams()) # TRUE
#' ## use defaults, change max number of iteration of MCMC
#' par <- bayesparams(maxit=1e5)
#' @export
bayesparams <- function(prop_a = 0.02,
                        prop_b = 0.02,
                        prior_mu = c(0,10),
                        prior_nu = c(2,1/2),
                        prior_eta = c(2,2),
                        trunc = 100,
                        comp_saved = 15,
                        maxit = 30000,
                        burn = 5000,
                        thin = 1,
                        adapt = 5000,
                        batch_size = 125,
                        start_ab = c("guesstimate", "prior"),
                        mode = 1,
                        prop.a = deprecated(),
                        prop.b = deprecated(),
                        prior.mu = deprecated(),
                        prior.nu = deprecated(),
                        prior.eta = deprecated(),
                        comp.saved = deprecated(),
                        batch.size = deprecated()) {
  if (lifecycle::is_present(prop.a)) {
    lifecycle::deprecate_warn("0.4.0", "thetaruns(prop.a)", "thetaruns(prop_a)")
    prop_a <- prop.a
  }
  if (lifecycle::is_present(prop.b)) {
    lifecycle::deprecate_warn("0.4.0", "thetaruns(prop.b)", "thetaruns(prop_b)")
    prop_b <- prop.b
  }
  if (lifecycle::is_present(prior.mu)) {
    lifecycle::deprecate_warn("0.4.0", "thetaruns(prior.mu)",
                              "thetaruns(prior_mu)")
    prior_mu <- prior.mu
  }
  if (lifecycle::is_present(prior.nu)) {
    lifecycle::deprecate_warn("0.4.0", "thetaruns(prior.nu)",
                              "thetaruns(prior_nu)")
    prior_nu <- prior.nu
  }
  if (lifecycle::is_present(prior.eta)) {
    lifecycle::deprecate_warn("0.4.0", "thetaruns(prior.eta)",
                              "thetaruns(prior_eta)")
    prior_eta <- prior.eta
  }
  if (lifecycle::is_present(comp.saved)) {
    lifecycle::deprecate_warn("0.4.0", "thetaruns(comp.saved)",
                              "thetaruns(comp_saved)")
    comp_saved <- comp.saved
  }
  if (lifecycle::is_present(batch.size)) {
    lifecycle::deprecate_warn("0.4.0", "thetaruns(batch.size)",
                              "thetaruns(batch_size)")
    batch_size <- batch.size
  }
  x <- list(prop_a = prop_a,
            prop_b = prop_b,
            prior_mu = prior_mu,
            prior_nu = prior_nu,
            prior_eta = prior_eta,
            trunc = trunc,
            comp_saved = comp_saved,
            maxit = maxit,
            burn = burn,
            thin = thin,
            adapt = adapt,
            batch_size = batch_size,
            start_ab = start_ab[1],
            mode = mode)
  class(x) <- "bayesparams"
  if (is.bayesparams(x)) return(x)
  else stop("Wrong type of argument. See help for details.")
}

## validator
#' @export
is.bayesparams <- function(x) {
  if (length(names(x))) {
    deprecated_names <- c("prop.a", "prop.b", "prior.mu", "prior.eta",
                          "comp.saved", "batch.size")
    if (any(deprecated_names %in% names(x))) {
      lifecycle::deprecate_warn(
        when = "0.4.0",
        what = I('Dot-notation elements such as "prop.a"'),
        with = I('snake-notation such as "prop_a"'))
      dn <- deprecated_names[deprecated_names %in% names(x)]
      dn_new <- sub(".", "_", dn)
      x[dn_new] <- x[dn]
    }
    names_bp <- c("prop_a","prop_b","prior_mu","prior_nu","prior_eta",
                  "trunc","comp_saved","maxit","burn","thin",
                  "adapt","batch_size","start_ab","mode")
    st_ab <- c("guesstimate", "prior")
    conds <- names_bp %in% names(x)
    if (all(conds)) {
      conds <- is.numeric(x$prop_a) && is.numeric(x$prop_b) &&
        length(x$prop_a)==1 && length(x$prop_b)==1 && x$prop_a>0 && x$prop_b>0 &&
        is.numeric(x$prior_mu) && is.numeric(x$prior_nu) && is.numeric(x$prior_eta) &&
        length(x$prior_mu)==2 && length(x$prior_nu)==2 && length(x$prior_eta)==2 &&
        x$prior_mu[2] > 0 && all(x$prior_nu>0) && all(x$prior_eta>0) &&
        length(x$trunc)==1 && is_int(x$trunc) && length(x$comp_saved)==1 && is_int(x$comp_saved) &&
        length(x$maxit)==1 && is_int(x$maxit) && length(x$burn)==1 && is_int(x$burn) &&
        length(x$thin)==1 && is_int(x$burn) && length(x$adapt)==1 && is_int(x$adapt) &&
        length(x$batch_size)==1 && is_int(x$batch_size) &&
        is.character(x$start_ab) && x$start_ab %in% st_ab &&
        length(mode)==1 && x$mode %in% c(0,1,2)
      return(inherits(x, "bayesparams") && conds)
    }     
  }
  return(FALSE)
}

## viewer
#' @export
print.bayesparams <- function(x, ...) {
  summary.bayesparams(x, ...)
}

#' @export
summary.bayesparams <- function(object, ...) {
  if (!is.bayesparams(object)) {
    warning('"object" does not comply with the "bayesparams" class')
    summary.default(object, ...)
  }else{
    cat("proposal for alpha              ",object$prop_a,"\n", sep="")
    cat("proposal for beta               ",object$prop_b,"\n", sep="")
    cat("start values for alpha and beta ",ifelse(object$start_ab[1]=="prior", "drawn from prior", "guesstimated from 2-step fit"),"\n", sep="")
    cat("prior par. for means            (",object$prior_mu[1],", ",object$prior_mu[2],")\n", sep="")
    cat("prior par. for sd               (",object$prior_nu[1],", ",object$prior_nu[2],")\n", sep="")
    cat("prior par. for precision par.   (",object$prior_eta[1],", ",object$prior_eta[2],")\n", sep="")
    cat("DP truncation                   ",object$trunc,"\n", sep="")
    cat("nbr of components saved         ",object$comp_saved,"\n", sep="")
    cat("max. nbr of MCMC iterations     ",object$maxit,"\n", sep="")
    cat("length of burn-in               ",object$burn,"\n", sep="")
    cat("thinning par.                   ",object$thin,"\n", sep="")
    cat("length of adaption              ",object$adapt,"\n", sep="")
    cat("size of adaption batches        ",object$batch_size,"\n", sep="")
    cat("debug mode                      ",object$mode,"\n", sep="")
  }
  invisible()
}


##################################################
## CLASS "DEPMEASURE"
#' Dependence measures estimates
#' 
#' Test or show objects of class "depmeasure", a structure for all dependence
#' measures covered in this package.
#' 
#' @param type one of "theta", "chi", "steptheta" or "runs".
#' @param x an arbitrary \R object.
#' @returns An object of class "depmeasure" whose exact structure depends on
#'   its type. All types contain the following elements,
#'   \item{probs }{probability levels at which the dependence measure was estimated.}
#'   \item{levels }{corresponding quantile levels on the data's original scale.}
#'   \item{nlag }{number of lags.}
#'   
#'   Furthermore, if \code{type} is "theta" or "chi",
#'   \item{fit }{an object of class "bayesfit".}
#'   \item{distr }{a matrix with samples from the posterior distribution, with \code{length(levels)} rows and \code{S} columns (see [thetafit()]).}
#'   \item{theta,chi }{a matrix with the posterior estimates of the respective dependence measure at different thresholds \code{levels} (rows) and coverage intervals (named columns, first two for the mean and median).}
#'  
#'  If \code{type} is "steptheta",
#'  \item{fit }{an object of class "stepfit".}
#'  \item{theta }{a matrix with the estimate for theta at different thresholds \code{levels} (rows) and confidence levels (columns).}
#'    
#'  If \code{type} is "runs",
#'  \item{theta }{a matrix with the (empirical) runs estimator estimate and confidence intervals (columns) at the given \code{levels} (rows).}
#'  \item{nbr_exc }{the number of exceedances corresponding to the \code{levels}.}
#' @seealso [depmeasures()]
#' @export
depmeasure <- function(type=c("theta","chi","steptheta","runs")) {
  z_mat <- matrix(0, nrow=1, ncol=1)
  x <- list(probs=0, levels=0, nlag=1)
  if (type[1]=="theta" || type[1]=="chi") {
    x$fit   <- bayesfit()
    x$distr <- z_mat
    if (type[1]=="theta")
      x$theta <- z_mat
    else if (type[1]=="chi")
      x$chi <- z_mat
  } else if (type[1]=="steptheta") {
    x$fit <- stepfit()
    x$theta <- z_mat
  } else if (type[1]=="runs") {
    x$theta <- z_mat
    x$nbr_exc <- 0
  } else {
    stop("type not recognised. Must be either 'theta' or 'chi'.")
  }
  class(x) <- "depmeasure"
  return(x)
}

## validator
#' @export
is.depmeasure <- function(x) {
  if (length(names(x))) {
    conds <- names(depmeasure("theta")) %in% names(x)
    if (all(conds)) {
      conds <- is.bayesfit(x$fit) && is.matrix(x$theta)
      return(inherits(x, "depmeasure") && conds)
    }
    conds <- names(depmeasure("chi")) %in% names(x)
    if (all(conds)) {
      conds <- is.bayesfit(x$fit) && is.matrix(x$chi)
      return(inherits(x, "depmeasure") && conds) 
    }
    conds <- names(depmeasure("steptheta")) %in% names(x)
    if (all(conds)) {
      conds <- is.stepfit(x$fit) && is.matrix(x$theta)
      return(inherits(x, "depmeasure") && conds) 
    }
    conds <- names(depmeasure("runs")) %in% names(x)
    if (all(conds)) {
      conds <- is.matrix(x$theta)
      return(inherits(x, "depmeasure") && conds)
    }
  }
  return(FALSE)
}

## viewer
#' @export
print.depmeasure <- function(x, ...) {
  summary.depmeasure(x, ...)
}

#' @export
summary.depmeasure <- function(object, ...) {
  if (!is.depmeasure(object)) {
    warning('"object" does not comply with the "depmeasure" class')
    summary.default(object, ...)
  }
  cat("Bayes fit\n")
  summary(object$fit)
  cat(paste(names(object)[2],"posterior estimates\n"))
  print(summary(object[[2]]))
  invisible()
}

## plotter
#' @rdname depmeasure
#' @export
plot.depmeasure <- function(x, ...) {
  stopifnot(is.depmeasure(x))
  if ("nbr_exc" %in% names(x)) {
    ord <- order(colnames(x$theta))
    def <- list(xlab="x", ylab=bquote(bquote(theta(x,.(x$nlag)))), ylim=range(x$theta))
    if (length(list(...)) > 0) {
      pm <- !pmatch(names(def), names(list(...)), 0)# all non-matches
      def <- def[pm]
    }
    def$x <- x$levels; def$y <- x$theta[,"estimate"]
    do.call(plot, c(type="l", def, list(...)))
    ord <- ord[-length(ord)]
    # add credible/confidence lines
    if (length(ord) > 1) {
      for (i in 0:(floor(length(ord)/2)-1)) {
        lines(x$levels, x$theta[,ord[i*2+1]], lty=i+2)
        lines(x$levels, x$theta[,ord[i*2+2]], lty=i+2)
      }
    }
  } else if ("chi" %in% names(x)) {
    ord <- order(colnames(x$chi))
    def <- list(xlab="x", ylab=bquote(bquote(chi[.(x$nlag)](x))), ylim=range(x$chi))
    if (length(list(...)) > 0) {
      pm <- !pmatch(names(def), names(list(...)), 0)
      def <- def[pm]
    }
    def$x <- x$levels; def$y <- x$chi[,"median"]
    do.call(plot, c(def,type="l",list(...)))
    ord <- ord[-c(length(ord)-1,length(ord))]
    # add credible/confidence lines
    if (length(ord) > 1) {
      for (i in 0:(floor(length(ord)/2)-1)) {
        lines(x$levels, x$chi[,ord[i*2+1]], lty=i+2)
        lines(x$levels, x$chi[,ord[i*2+2]], lty=i+2)
      }
    }
  } else {
    ord <- order(colnames(x$theta))
    def <- list(xlab="x", ylab=bquote(bquote(theta(x,.(x$nlag)))), ylim=range(x$theta))
    if (length(list(...)) > 0) {
      pm <- !pmatch(names(def), names(list(...)), 0)
      def <- def[pm]
    }
    def$x <- x$levels; ifelse(is.stepfit(x$fit), def$y <- x$theta[,"estimate"], def$y <- x$theta[,"median"])
    do.call(plot, c(def,type="l",list(...)))
    ord <- ord[-length(ord)]
    if (is.bayesfit(x$fit)) ord <- ord[-length(ord)]
    # add credible/confidence lines
    if (length(ord) > 1) {
      for (i in 0:(floor(length(ord)/2)-1)) {
        lines(x$levels, x$theta[,ord[i*2+1]], lty=i+2)
        lines(x$levels, x$theta[,ord[i*2+2]], lty=i+2)
      }
    }
  }
}


##################################################
## ADDITIONAL VALIDATORS
#' A loose integer validation function
#' 
#' @param x an arbitrary \R object.
#' @keywords internal
is_int <- function(x) {
  if (all(x%/%1==x)) return(TRUE)
  else return(FALSE)
}
