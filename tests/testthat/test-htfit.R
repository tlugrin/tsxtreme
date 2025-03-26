test_that("comp saved > stick-breaking truncation generates a warning", {
  ts <- stats::arima.sim(list(ar = 0.5), n = 1000)
  bpar <- bayesparams(maxit = 100, burn = 1, thin = 5, adapt = 0,
                      comp_saved = 10, trunc = 9, mode = 2)
  expect_warning(depfit(ts = ts, u_mar = 0.8, par = bpar),
                 "comp_saved.*trunc")
})

test_that("not meaningful maxit/burn/thin triad generates an error", {
  ts <- stats::arima.sim(list(ar = 0.5), n = 1000)
  bpar <- bayesparams(maxit = 100, burn = 101, thin = 1, adapt = 0, mode = 2)
  expect_error(depfit(ts = ts, u_mar = 0.8, par = bpar),
                 "maxit.*burn.*thin")
  bpar <- bayesparams(maxit = 100, burn = 100, thin = 1, adapt = 0, mode = 2)
  expect_error(depfit(ts = ts, u_mar = 0.8, par = bpar),
                 "maxit.*burn.*thin")
  bpar <- bayesparams(maxit = 100, burn = 96, thin = 5, adapt = 0, mode = 2)
  expect_error(depfit(ts = ts, u_mar = 0.8, par = bpar),
                 "maxit.*burn.*thin")
  bpar <- bayesparams(maxit = 100, burn = 10, thin = 0, adapt = 0, mode = 2)
  expect_error(depfit(ts = ts, u_mar = 0.8, par = bpar),
                 "thin")
})

# not meant to be exhaustive but should provide good coverage
test_that("depfit generates no errors across non-default parameter combinations", {
  # simulate time series
  set.seed(1)
  n <- 5000
  dep <- 0.7
  ts <- numeric(n)
  ts[1] <- rnorm(1, sd = dep)
  for (i in 2:n) ts[i] <- rnorm(1, mean=dep*ts[i-1], sd=1-dep^2)
  ts_L <- qlapl(pnorm(ts, sd = dep))
  # set up bayesparams object
  bpar <- bayesparams(prop_a = 0.1, #!=defaults
                      prop_b = 1,
                      prior_mu = c(2, 5),
                      prior_nu = c(1, 1),
                      prior_eta = c(1/2, 2),
                      batch_size = 1,
                      trunc = 9,
                      adapt = 0,
                      comp_saved = 3,
                      maxit = 99,
                      burn = 42,
                      thin = 7,
                      mode = 2)
  bpar$start_ab <- "guesstimate"
  bpar$conditions <- TRUE
  expect_no_error(depfit(ts = ts, u_mar = 0.9, u_dep = 0.81, lapl = FALSE,
                      method_mar = "mom", nlag = 2, par = bpar, submodel = "fom"))
  expect_no_error(depfit(ts = ts, u_mar = 0.9, u_dep = 0.82, lapl = FALSE,
                      method_mar = "mom", nlag = 2, par = bpar, submodel = "none"))
  bpar$conditions <- FALSE
  expect_error(depfit(ts = ts, u_mar = 0.9, u_dep = 0.83, lapl = FALSE,
                      method_mar = "pwm", nlag = 2, par = bpar, submodel = "fom"),
               "NaNs produced when rescaling")
  bpar$start_ab <- "prior"
  expect_error(depfit(ts = ts, u_mar = 0.9, u_dep = 0.84, lapl = FALSE,
                      method_mar = "pwm", nlag = 2, par = bpar, submodel = "none"),
               "NaNs produced when rescaling")
  expect_no_error(depfit(ts = ts_L, u_mar = 0.9, u_dep = 0.85, lapl = TRUE,
                         method_mar = "pwm", nlag = 2, par = bpar, submodel = "fom"))
 
})

test_that("depfit generates no errors on univariate Gaussian model", {
  data <- MASS::galaxies/1000
  bpar <- bayesparams(maxit = 100, burn = 1, thin = 5, adapt = 0, mode = 2)
  expect_no_error(depfit(ts = data, par = bpar, submodel = "ugm"))
})

test_that("depfit on silent mode is silent", {
  ts <- rlapl(500)
  bpar <- bayesparams(maxit = 50, burn = 0, thin = 1, adapt = 10, mode = 2,
                      conditions = FALSE, trunc = 11, comp_saved = 9)
  expect_no_message(depfit(ts = ts, par = bpar, lapl = TRUE, submodel = "ugm"))
  dep <- 0.8
  ar <- rnorm(1000, sd = dep)
  for (i in seq_along(ar)[-1]) ar[i] <- rnorm(1, mean = ar[i-1]*dep, sd = 1-dep^2)
  expect_no_message(depfit(ts = ar, u_mar = 0.8, par = bpar, lapl = FALSE,
                           submodel = "none"))
})
