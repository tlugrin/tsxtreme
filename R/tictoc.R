## Copyright (C) 2017-2025 Thomas Lugrin
## tic and toc functions are provided by tictoc package,
## but need a wrapper around the latter
##################################################

#####################################
## CLOCK FUNCTIONS
#' Timing function
#' 
#' Wrapper around [tictoc::tic()] and [tictoc::toc()]
#' 
#' @param silent boolean, should the function be verbose (`FALSE`, default), or
#'   only return a time value invisibly?
#' @returns A time difference (invisibly).#' 
#' @keywords internal
toc <- function(silent = FALSE)
{
  tt <- tictoc::toc(quiet = TRUE)#list of 3
  ttic <- tt$tic
  ttoc <- tt$toc
  if(!silent) print(ttoc - ttic)
  invisible(ttoc - ttic)
}
