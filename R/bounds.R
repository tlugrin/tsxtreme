## Copyright (C) 2017-2025 Thomas Lugrin
## Compute Keef, Papastathopoulos & Tawn bounds on (a,b)
## (Scope:) Heffernan-Tawn model fit in 2 stages
## List of functions: - conditions_verify
##                    - conditions_bounds
##                    - dicho_b_bounds
##                    - dicho_a_bounds
##################################################

##################################################
#' Bounds on (alpha,beta)
#' 
#' Mainly called by [conditions_bounds()] to check whether a given (alpha, beta)
#' pair lies within the Keef et al. boundaries.
#' 
#' @param a scalar, alpha parameter, in [-1,1]
#' @param b scalar, beta parameter, in [0,1]
#' @param p scalar, probability (to compute quantiles of residual distribution), in [0,1]
#' @param v scalar, threshold, recommended as beyond the maximum observation
#' @param data bivariate vector, (X,Y) with Y | X>u
#' @param boolean, does the couple (a,b) satisfy the conditions?
#' @returns A boolean, TRUE if the (alpha, beta) pair is valid, FALSE otherwise.
#' @keywords internal
conditions_verify <- function(a, b, p, v = -log(2*(1-0.99999999)), data) {
  z_q_pos <- q_res2(p, 1, 0, data)
  z_q     <- q_res2(p, a, b, data)
  z_q_neg <- q_res2(p, -1, 0, data)
  if (a < -1 || a > 1 || b >= 1 || b<0) { return(FALSE) }
  
  # Case I
  case1 <- FALSE
  if (a <= min(1,  1 - b*z_q*v^{b-1},  1 - v^{b-1}*z_q + z_q_pos/v)) { case_1 <- TRUE }
  else {
    cond1 <- (1 - b*z_q*v^{b-1} < a)
    if (b*z_q>0) cond2 <- ( (1-1/b) * exp( log(b*z_q)*(1/(1-b)) + log(1-a)*(-b/(1-b)) ) + z_q_pos > 0 )
    else        cond2 <- ( (1-1/b) * (b*z_q)^{1/(1-b)} * (1-a)^{-b/(1-b)} + z_q_pos > 0 )
    if (cond1 && cond2) { case1 <- TRUE }
    else { return(FALSE) }
  }
  # Case II
  case2 <- FALSE
  if (-a <= min(1 + b*z_q*v^{b-1}, 1 + v^{b-1}*z_q - z_q_neg/v)) {
    case2 <- TRUE
    return(TRUE)
  } else {
    cond1 <- ( 1 + b*v^{b-1}*z_q < -a )
    if (-b*z_q>0) cond2 <- ( (1-1/b) * exp( log(-b*z_q)*(1/(1-b)) + log(1+a)*(-b/(1-b)) ) - z_q_neg > 0 )
    else         cond2 <- ( (1-1/b) * (-b*z_q)^{1/(1-b)} * (1+a)^{-b/(1-b)} - z_q_neg > 0 )
    if (cond1 && cond2) { case2 <- TRUE; return(TRUE) }
    else { return(FALSE) }
  }
}

#' Compute bounds for alpha or beta
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' Called by [ht2fit()] to get the support of the (alpha, beta) pair during
#' initialisation. Exposed to the user for their own curiosity, at their own
#' risk!
#' 
#' @details
#' The function computes bounds either for \eqn{\alpha} and requires a (fixed)
#' value for \eqn{\beta}, or the other way around.
#' A twin C routine supports the MCMC algorithm when computing
#' the normalising constants for the truncated normal proposals of \eqn{\alpha}
#' and \eqn{\beta}.
#'
#' @param fix_par scalar, either the value of \eqn{\alpha} or that of \eqn{\beta}.
#' @param fix_alpha boolean, `TRUE` if [fix_par] is for \eqn{\alpha}, `FALSE`
#'   if [fix_par] is for \eqn{\beta}.
#' @param data 2-column matrix representing the bivariate vector
#'   \eqn{(X,Y)} with \eqn{Y | X>u}.
#' @param eps tolerance parameter for finding the boundary value,
#'   defaults to 0.001.
#' @param v value of a very high quantile on the Laplace scale.
#' @returns A numeric vector of length 2 with the lower and upper bounds for the
#'   given parameter (either alpha or beta).
#' @export 
conditions_bounds <- function(fix_par, fix_alpha=TRUE, data, eps=0.001, v=-log(2*(1-0.99999999))) {
  i <- o <- 2
  bounds <- numeric(2)
  if (fix_alpha) {
    if (conditions_verify(fix_par, 0, 1, v, data) &&
         conditions_verify(fix_par, 0, 0, v, data)) { i <- 0 }
    else { o <- 0 }
    if (conditions_verify(fix_par, 1, 1, v, data) &&
         conditions_verify(fix_par, 1, 0, v, data)) { i <- 1 }
    else { o <- 1 }
    if (i == 2) {
      grid <- seq(0, 1, 0.05)
      pt   <- 1
      while(i == 2 && pt <= length(grid)) {
        if (conditions_verify(fix_par, grid[pt], 1, v, data) &&
             conditions_verify(fix_par, grid[pt], 0, v, data)) { i <- grid[pt] }
        pt <- pt+1
      }
      if (i == 2) { stop("Refine grid on beta to get proper results... in [conditions_bounds].") }
    } else if (o == 2) {
      bounds[1] <- -1
      bounds[2] <- 1
      return(bounds)
    }
    bounds[1] <- i
    bounds[2] <- dicho_b_bounds(i, o, fix_par, data, eps, v)
  } else {
    if (conditions_verify(-1, fix_par, 1, v, data) &&
         conditions_verify(-1, fix_par, 0, v, data)) { i <- -1 }
    else { o <- -1 }
    if (conditions_verify(1, fix_par, 1, v, data) &&
         conditions_verify(1, fix_par, 0, v, data)) { i <- 1 }
    else { o <- 1 }
    if (i == 2) {
      grid <- seq(-1, 1, 0.05)
      pt   <- 1
      while(i == 2 && pt <= length(grid)) {
        if (conditions_verify(grid[pt], fix_par, 1, v, data) &&
             conditions_verify(grid[pt], fix_par, 0, v, data)) { i <- grid[pt] }
        pt <- pt+1
      }
      if (i == 2) { stop("Refine grid on alpha to get proper results... in [conditions_bounds].") }
    } else if (o == 2) {
      bounds[1] <- -1
      bounds[2] <- 1
      return(bounds)
    }
    bounds[1] <- i
    bounds[2] <- dicho_a_bounds(i, o, fix_par, data, eps, v)
  }
  if (bounds[1] > bounds[2]) {
    v         <- bounds[2]
    bounds[2] <- bounds[1]
    bounds[1] <- v
  }
  return(bounds)
}

#' Determine exact bound for alpha or beta
#' 
#' Pair of functions called by [conditions_bounds()] to get more accurate
#' bounds by dichotomy.
#' 
#' @param i current inner bound
#' @param o current outer bound
#' @param fix_par scalar, either value of alpha or of beta
#' @param data bivariate vector, (X,Y) with Y | X>u
#' @param eps scalar, tolerance, typically 0.001
#' @returns Scalar, the bound found.
#' @keywords internal
dicho_b_bounds <- function(i, o, fix_par, data, eps, v=-log(2*(1-0.99999999))) {
  mid <- mean(c(i,o))
  if (conditions_verify(fix_par, mid, 1, v, data) &&
       conditions_verify(fix_par, mid, 0, v, data)) { i<- mid }
  else { o <- mid }
  if (abs(i-o) < eps) { return(i) }
  else { return(dicho_b_bounds(i, o, fix_par, data, eps, v)) }
}

#' @rdname dicho_b_bounds
#' @keywords internal
dicho_a_bounds <- function(i, o, fix_par, data, eps, v=-log(2*(1-0.99999999))) {
  mid <- mean(c(i,o))
  
  if (conditions_verify(mid, fix_par, 1, v, data) &&
       conditions_verify(mid, fix_par, 0, v, data)) { i <- mid }
  else { o <- mid }
  
  if (abs(i-o) < eps) { return(i) }
  else return(dicho_a_bounds(i, o, fix_par, data, eps, v))
}