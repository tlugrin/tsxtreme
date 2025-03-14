## Copyright (C) 2017-2025 Thomas Lugrin
## Functions related to H+T model
## (Scope:) Heffernan-Tawn model fit in 1 stage
## List of functions: - p.res
##                    - r.res
##################################################

#' The conditional tail model residual distribution (Bayesian)
#' 
#' Distribution function and random generation for the residual distribution of
#' the Heffernan--Tawn conditional tail model. Intended for the Bayesian
#' semiparametric approach. The `p.res` method is called by
#' [etfit()] for Monte-Carlo integration, the `r.res` method is called by
#' [etfit()] for the proportion method. Not exposed to the user.
#' 
#' @param s integer, ranges across sweeps.
#' @param z array of reals, trace of residuals.
#' @param w matrix of reals, trace of the weights of components --- rows sum
#'   to 1.
#' @param mu,sig arrays of reals, means and standard deviations of the stored
#'   mixture components output from [depfit()] with dimensions
#'   _nbr of sweeps_ x _number of components_ x _nbr of lags_.
#' @param R integer>0, number of samples from a mixture distribution sampled at
#'   a given iteration.
#' @param S integer>0, number of posterior iterations used for the estimation.
#' @param nlag integer>0, number of lags in the original data.
#' @returns For `p.res`: a matrix of reals, evaluation of the Heffernan--Tawn
#'   residual distribution function. For `r.res`: a matrix of reals, samples
#'   of the residual distribution across a sample of iterations.
#' @keywords internal
p.res <- function(s, z, w, mu, sig){
  R   <- dim(z)[1]; nlag = dim(z)[3]
  ret <- numeric(R)
  for(r in 1:R){
    prev <- 1
    for(lag in 1:nlag){
      prev <- prev*sum(w[s,]*pnorm(z[r,s,lag], mu[s,,lag], sig[s,,lag]))
    }
    ret[r] <- prev
  }
  return(ret)
}

#' @rdname p.res
r.res <- function(R, S, nlag, w, mu, sig){
  sim.Z <- array(0, dim=c(R,S,nlag))
  
  for(sw in 1:S){# loop on selected sweeps
    hwmany <- rmultinom(1, R, w[sw,])
    r      <- 1
    c      <- 1
    
    while(r<R){# sample from mixture
      if(hwmany[c]>0){
        if(dim(s)[3]==1){
          sim.Z[r:(r+hwmany[c]-1),sw,] <- rnorm(hwmany[c], m[sw,c,], s[sw,c,])
        }else{
          sim.Z[r:(r+hwmany[c]-1),sw,] <- rmvnorm(hwmany[c], m[sw,c,], diag(s[sw,c,]^2))
        }
      }
      r <- r+hwmany[c]
      c <- c+1
    }
  }
  sim.Z <- matrix(sim.Z, nrow=R, ncol=S*nlag)
  return(sim.Z)
}
