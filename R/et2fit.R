## Copyright (C) 2017-2025 Thomas Lugrin
## user interface for E+T fit
## (Scope:) Model fit in 2 stages
## List of functions: - theta2fit
##                    - th2est
##                    - thboot
##                    - blockboot
##################################################

##################################################
## FIT

#' Estimate extremal index using a stepwise approach
#' 
#' The measure of dependence \eqn{\theta(x,m)} is derived using a method
#' described in Eastoe and Tawn (2012). Appropriate marginal transforms can be
#' performed automatically or separately.
#' 
#' @details
#' The standard procedure (\code{method=="prop"}) to estimating probabilities
#' from a Heffernan-Tawn fit best illustrated in the bivariate context
#' (\eqn{Y\mid X>u}{Y | X > u}):
#' 
#' 1. sample \eqn{X} from an exponential distribution above \eqn{v \ge u},
#' 2. sample \eqn{Z} (residuals) from their empirical distribution,
#' 3. compute \eqn{Y} using the relation
#'   \eqn{Y = \alpha\times X + X^\beta\times Z}{Y = \alpha * X + X^\beta * Z},
#' 4. estimate \eqn{Pr(X > v_x, Y > v_y)} by calculating the proportion \eqn{p}
#'   of \eqn{Y} samples above \eqn{v_y} and multiply \eqn{p} with the marginal
#'   survival distribution evaluated at \eqn{v_x}.
#'   
#' With \code{method="MCi"} a Monte Carlo integration approach is used, where
#' the survivor distribution of \eqn{Z} is evaluated at pseudo-residuals of the
#' form
#' \deqn{\frac{v_y - \alpha\times X}{X^\beta},}{(v_y - \alpha * X)/X^\beta,}
#' where \eqn{X} is sampled from an exponential distribution above \eqn{v_x}.
#' Taking the mean of these survival probabilities, we get the Monte Carlo
#' equivalent of \eqn{p} in the proportion approach.
#' 
#' @param ts numeric vector; time series to be fitted.
#' @param lapl logical; is \code{ts} on the Laplace scale already? The default
#'   (FALSE) assumes unknown marginal distribution.
#' @param nlag integer; number of lags to be considered when modelling the
#'   dependence in time.
#' @param R integer; the number of samples used for estimating \eqn{\theta(x,m)}.
#' @param u.mar marginal threshold; used when transforming the time series to
#'   Laplace scale if \code{lapl} is FALSE; not used otherwise.
#' @param u.dep dependence threshold; level above which the dependence is
#'   modelled. \code{u.dep} can be lower than \code{u.mar}.
#' @param probs vector of probabilities; the values of \eqn{x} for which to
#'   evaluate \eqn{\theta(x,m)}.
#' @param method.mar a character string defining the method used to estimate the
#'   marginal GPD; either \code{"mle"} for maximum likelihood of \code{"mom"}
#'   for method of moments or \code{"pwm"} for probability weighted moments
#'   methods. Defaults to \code{"mle"}.
#' @param method a character string defining the method used to estimate the
#'   dependence measure; either \code{"prop"} for proportions or \code{"MCi"}
#'   for Monte Carlo integration (see Details).
#' @param silent logical (\code{FALSE}); verbosity.
#' @param R.boot integer; the number of samples used for the block bootstrap for
#'   the confidence intervals.
#' @param block.length integer; the block length used for the block-bootstrapped
#'   confidence intervals.
#' @param levels vector of probabilities; the quantiles of the bootstrap
#'   distribution of the extremal measure to be computed.
#' @returns An object of class [depmeasure()] of `type` "steptheta".
#' @seealso [dep2fit()], [thetafit()], [thetaruns()]
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
#' fit <- theta2fit(ts=ar, u.mar=0.95, u.dep=0.98)
#' 
#' ## plot theta(x,1)
#' plot(fit)
#' abline(h=1, lty="dotted")
#' @export
theta2fit <- function(ts, lapl=FALSE, nlag=1, R=1000,
                         u.mar=0, u.dep, probs=seq(u.dep,0.9999,length.out=30),
                         method.mar=c("mle","mom","pwm"), method=c("prop","MCi"),
                         silent=FALSE,
                         R.boot=0, block.length=m*5, levels=c(.025,.975)){
  data.up <- format.ts(ts=ts, u.mar=u.mar, u.dep=u.dep, method=method.mar,
                       nlag=nlag, lapl=lapl)
  m       <- nlag+1
  n       <- dim(data.up)[1]
  mesh    <- probs
  mesh.O  <- scale.to.original(p=mesh, ts=ts, u=u.mar, gpd.pars=gpd(ts, u=u.mar, method=method.mar)$pars)
  
  fit <- ht2step(data.up)
  if(!silent)
    print(paste("Likelihood fit complete. Estimates for alpha, beta resp.:",
                signif(fit$a,3),"and",signif(fit$b,3)))
  n.vert <- length(mesh)
  sim.U  <- runif(R)
  ind    <- sample.int(n, size=R, replace=TRUE)
  sim.Z  <- fit$res[ind,, drop=FALSE]
  #compute theta for each m at each probability in mesh
  nbr.quant <- length(levels)+1
  th        <- matrix(NaN, nrow=n.vert, ncol=nbr.quant,
                    dimnames=list(NULL,c("estimate",paste(levels*100,"%",sep=""))))
  th[,"estimate"] <- th2est(fit, sim.U=sim.U, sim.Z=sim.Z, mesh=mesh, nlag=nlag, method=method)
  ret <- depmeasure("steptheta")
  ret$fit <- fit; ret$probs <- probs; ret$levels <- mesh.O; ret$nlag <- nlag
  ## compute bootstrap intervals (recursive)
  if(R.boot > 0 & nbr.quant > 1){
    boot <- thboot(par=list(ts=ts,u.mar=u.mar,u.dep=u.dep,nlag=nlag,method.mar=method.mar,R=R,mesh=mesh),
                      block.length=block.length, R=R.boot, method=method, levels=levels)
    n.vert <- dim(boot)[2]
    th[1:n.vert,2:nbr.quant] <- boot
  }else{
    th <- th[,"estimate",drop=FALSE]
  }
  ret$theta <- th
  return(ret)
}

#' Estimate the extremal index from a Heffernan--Tawn model fit
#' 
#' Simulate residuals from the fitted Heffernan--Tawn model and build (X,Y)
#' pairs with \eqn{X>u} to derive an estimate of \eqn{\theta(x,m)}.
#' Called by [theta2fit()], not exposed to the user.
#' 
#' @param fit output from [ht2step()], an object of class [stepfit()].
#' @param sim.U vector of reals, `R` uniform samples (see [ht2step()]).
#' @param sim.Z matrix of reals, `R` multivariate samples from the
#'   Heffernan--Tawn residual distribution.
#' @param mesh vector of reals, quantile values at which to evaluate
#'   \eqn{\theta} (uniform scale).
#' @param nlag number of lags \eqn{m} to consider
#' @param method (vector of) string(s), either "prop" or "MCi" for the method
#'   of proportions (empirical approach) or Monte-Carlo integration.
#' @returns Estimate of \eqn{\theta(x,m)} at each value of \eqn{x} given by
#'   [mesh].
#' @keywords internal
th2est <- function(fit, sim.U, sim.Z, mesh, nlag, method){
  n.vert <- length(mesh)
  R      <- length(sim.U)
  mesh.L <- qexp(1-(1-mesh)*2)
  nbr.quant <- length(levels)+1
  theta <- numeric(n.vert)
  # prepare MC integration...
  HZ.u   <- numeric(R)
  sort.Z <- apply(fit$res, 2, sort)
  for(i in 1:n.vert){
    sim.L <- -log(2*(1-mesh[i]-(1-mesh[i])*sim.U))
    # method of proportion
    if(sum(grepl("prop",method))){
      sim.Y <- as.matrix(sim.L%*%t(fit$a) + exp(log(sim.L)%*%t(fit$b))*sim.Z)# [Rx(m-1)]
      theta[i] <- sum(apply(as.matrix(sim.Y), 1, max)<mesh.L[i])/R
    }# Monte Carlo integration
    else if(sum(grepl("MCi",method))){
      sim.Z.MC <- (mesh.L[i]-sim.L%*%t(fit$a))/exp(log(sim.L)%*%t(fit$b))# [Rx(m-1)]
      HZ <- 1
      for(lag in 1:nlag){
        HZ <- HZ * p.res2(sim.Z.MC[,lag], sort.Z[,lag])
        if(i == 1){ HZ.u[lag,] <- HZ }
      }
      theta[i] <- mean(HZ)
    }else{
      stop("Unknown method. Should be either 'prop' or 'MCi'")
    }
  }
  return(theta)
}

#' Bootstrap confidence interval for the extremal index
#' 
#' Compute bootstrap confidence intervals for the extremal index estimated from
#' the stepwise fit of the conditional tail model. Called by [th2est()], not
#' exposed to the user.
#' 
#' @param par list, arguments to feed [theta2fit()].
#' @param block.length integer>0, length of blocks in block-bootstrap.
#' @param R integer>0, number of bootstrap samples to compute the confidence
#'   intervals.
#' @param method string, either "prop" or "MCi".
#' @param levels vector of probabilities, confidence levels, quantiles of
#'   bootstrap distribution to be output.
#' @returns A matrix where columns are quantiles and rows are x-levels of
#'   \eqn{\theta(x,m)}.
#' @keywords internal
thboot <- function(par, block.length, R, method, levels){
  nbr.vert <- length(par$mesh)
  ts.boot  <- blockboot(par$ts, block.length, R)
  n        <- block.length * floor(length(par$ts)/block.length)
  theta.R  <- apply(ts.boot, 2, function(ts,u.mar,u.dep,nlag,R,mesh,method.mar,method,shape.null.hyp,silent){
                                fit=theta2fit(ts=ts,u.mar=u.mar,u.dep=u.dep,
                                              lapl=TRUE,nlag=nlag,R=R,probs=mesh,
                                              method.mar=method.mar,method=method,silent=silent)
                                return(fit$theta[,"estimate"])}, 
                  nlag=par$nlag, R=par$R, u.mar=par$u.mar, u.dep=par$u.dep, mesh=par$mesh,
                  method.mar=par$method.mar, method=method, silent=TRUE)
  theta.R <- matrix(theta.R, nrow=nbr.vert, ncol=R)
  # remove NA values
  v <- 1
  while(v <= nbr.vert && sum(is.na(theta.R[v,])) == 0){ v <- v+1 }
  if(v == 1){
    warning("unable to compute uncertainty of theta")
    v <- 2
  }
  nbr.vert <- v-1
  theta.R  <- theta.R[1:nbr.vert,, drop=FALSE]
  # compute CI
  ci <- t(apply(theta.R, 1, quantile, levels))
  return(ci)
}

#' Basis method to block bootstrap a time series
#' 
#' Called by [thboot()], not exposed to the user.
#' 
#' @inheritParams thboot
#' @returns A matrix, [R] columns containing the block bootstrap samples of [ts]
#' @keywords internal
blockboot <- function(ts, block.length, R){
  rest <- length(ts)%%block.length
  if(rest > 0){ ts <- ts[-((length(ts)-rest+1):length(ts))] }
  ts.mat <- matrix(ts, nrow=block.length)
  nbr.bl <- floor(length(ts)/block.length)# == dim(ts.mat)[2]
  ts.R   <- matrix(0, nrow=nbr.bl*block.length, ncol=R)
  for(i in 1:R){
    sample   <- sample.int(nbr.bl, size=nbr.bl, replace=TRUE)
    ts.R[,i] <- c(ts.mat[,sample])
  }
  return(ts.R)
}
