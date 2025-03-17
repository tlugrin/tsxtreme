## Copyright (C) 2017-2025 Thomas Lugrin
## Functions related to H+T model
## (Scope:) Heffernan-Tawn model fit in 2 stages
## List of functions: - p_res2
##                    - q_res2
##################################################

#' The conditional tail model residual distribution (stepwise)
#' 
#' The `p_res2` method is called by [th2est()], the `q_res2` method is called
#' by [verifies_conditions()] in order to get quantiles of the empirical
#' distribution function of the residual distribution given
#' \eqn{(\alpha,\beta)}. Not exposed to the user.
#' 
#' @param res vector of reals, quantiles at which to compute the distribution
#'   function of the residuals.
#' @param sorted_res vector, empirical distribution function.
#' @param p scalar, probability (to compute quantiles of residual distribution),
#'   in [0,1].
#' @param a scalar, alpha parameter, in [-1,1].
#' @param b scalar, beta parameter, in [0,1].
#' @param data bivariate vector, \eqn{(X,Y)} with \eqn{Y | X>u}.
#' @returns For `p_res2`: a vector of the same length as `res`, distribution of
#'   the residuals evaluated in `res`. For `q_res2`: a quantile of the residual
#'   distribution.
#' @keywords internal
p_res2 <- function(res, sorted_res){
  n <- length(sorted_res)
  a <- approx(sorted_res, y = (1:n)/(n+1), xout = res, method = "linear",
              ties = "ordered", rule = 2)$y
  return(a)
}


#' @rdname p_res2
#' @keywords internal
q_res2 <- function(p, a, b, data){
  if(dim(data)[2] != 2)
    stop("data are in the wrong format. Provide a 2-column matrix for (X,y), with Y|X>u")
  Z <- (data[,2] - a*data[,1])/data[,1]^b
  return(quantile(Z, p))
}
