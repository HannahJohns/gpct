#' Permutation Tests
#'
#' @description Performs permutation testing
#'
#' @usage
#' winprop(x, splitTies=FALSE)
#'
#' @param x an object of type \code{pairwise}.
#' @param test an character vector giving the names of the tests to run \code{pairwise}.
#' @param nSims a positive number giving the number of permutations to run
#' @param splitTies A logical vector specifying how ties should be treated.
#'                  Should be \code{TRUE} for WMW Odds,
#'                  or \code{FALSE} for Agresti's GenOR.
#'
#' @details
#'
#'
#'
#' @examples
#'
#' @references
#'
#'
#'
#' @export
#'
permutationTest <- function(x, explanatory_var,
                            nSims=5000,
                             splitTies=TRUE, numThreads=NULL
                             )
{

  if(is.null(numThreads)){
    numThreads <- Inf
  }
  numThreads <- min(numThreads,RcppParallel::defaultNumThreads())
  if(numThreads < 1){
    numThreads <- ceiling(numThreads * RcppParallel::defaultNumThreads())
  }
  RcppParallel::setThreadOptions(numThreads = numThreads)

  if(prod(x$weight==1)!=1){stop("Permutation testing requires all observations to have weight of 1")}

  n <- nrow(x$comparisons)

  weights <- x$weight
  groupVector <- as.numeric(x$attributes[,explanatory_var])

  tol <- min(diff(sort(unique(groupVector))))/2

  shuffles <- do.call("cbind",lapply(1:n,function(i){floor(runif(nSims,min=i,max=n))}))
  results <- gpct:::get_permutationTest(x$comparisons,groupVector,weights,tol=tol,splitTies = splitTies,nSims = nSims,shuffles = shuffles)

  pVal <- mean(abs(log(results$result))<abs(log(results$result_null)))


  return(c(result=results$result,pVal=pVal))
}
