#' Inconsistency
#'
#' @description Calculates measures of preference inconsistency for pairwise preferences
#'
#' @usage
#' inconsistency(x, numThreads=NULL, includeTies=TRUE)
#'
#' @param x an object of type \code{pairwise}.
#' @param numThreads The number of threads to use. If \code{numThreads<1}, the method uses that proportion of maximum available threads (rounded up).
#' @param includeTies
#' @details
#'
#' @export
inconsistency <- function(x, numThreads=NULL, includeTies=TRUE,m=3)
{

  # Convert NA values into magic number for C++
  x$comparisons[is.na(x$comparisons)] <- -99999

  if(is.null(numThreads))
  {
    numThreads <- Inf
  }

  numThreads <- min(numThreads,RcppParallel::defaultNumThreads())
  if(numThreads < 1)
  {
    numThreads <- ceiling(numThreads * RcppParallel::defaultNumThreads())
  }
  RcppParallel::setThreadOptions(numThreads = numThreads)

  out <- gpct:::countInconsistency(x$comparisons,x$weight,includeTies=includeTies, m=m)

  return(out)
}
