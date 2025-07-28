#' Pairwise Comparison (composite survival)
#'
#' @aliases pairwise_surv
#'
#' @description Calculates pairwise comparison for ordered survival data
#'
#' @usage
#' pairwise_surv(times, statuses, weight=NULL, attributes=NULL, data, simplify=FALSE, numThreads=1)
#'
#' @param times A character vector giving column names that contain event times.
#' @param statuses A character vector giving column names that contain event statuses (1 for event occurred, 0 for censored)
#' @param weight An optional character string giving the column name that contains the weight associated with each observation
#' @param attributes A character vector giving column names that contain attributes of each patient
#' @param data A data frame, containing the specified columns
#' @param numThreads The number of threads to use. If \code{numThreads<1}, the method uses that proportion of maximum available threads (rounded up).
#' @param simplify A boolean that currently does nothing but will eventually optionally simplify data down for faster computation
#'
#' @return An object of type \code{pairwise} containing the following:
#' \describe{
#'     \item{comparisons}{A square matrix containing all-to-all comparisons of outcomes}
#'     \item{group}{A vector giving the groups associated with each row}
#'     \item{weight}{The weight of each observation}
#'     \item{attributes}{A data frame containing attributes for each observation}
#'  }
#'
#' @details
#'
#' This function appears to perform best when \code{numThreads=1}.
#' This could be overhead outweighing benefits, or it could be something to do with linking to other functions? idk
#'
#'
#' @export
pairwise_surv <- function(times, statuses, weight=NULL, attributes=NULL, data,
                          direction=NULL,
                          censored_as_tie = TRUE, minDiff = 0,
                          numThreads=NULL)
{

  # tibbles break things
  data <- as.data.frame(data)

  if( censored_as_tie != TRUE |  minDiff != 0) {
    stop("censored.at.tie and minDiff aren't working yet, do not use them!")
  }

  outcome <- as.matrix(data[,times])

  outcome_exists <- matrix(0L,nrow(data),length(statuses))
  outcome_exists[data[,statuses]==1] <- 1L

  if(is.null(numThreads)){
    numThreads <- 1
  }
  numThreads <- min(numThreads,RcppParallel::defaultNumThreads())
  if(numThreads < 1){
    numThreads <- ceiling(numThreads * RcppParallel::defaultNumThreads())
  }
  RcppParallel::setThreadOptions(numThreads = numThreads)

  if(is.null(direction)){
    direction <- rep(1,length(times))
  }

  out <- {}
  out$comparisons <- gpct:::get_pw_surv(outcome_exists = outcome_exists,outcome = outcome,
                                           direction = direction, censored_as_tie=censored_as_tie,
                                           minDiff=minDiff)

  if(is.null(weight)){
    out$weight <- rep(1,nrow(data))
  } else {
    out$weight <- data[,weight]
  }

  tmpdf <- {}
  for(cov in attributes){
    if(is.null(tmpdf)) {
      tmpdf <- as.data.frame(data[,cov])
    } else {
      tmpdf <- cbind(tmpdf,as.data.frame(data[,cov]))
    }
    colnames(tmpdf)[ncol(tmpdf)] <- cov
  }

  out$attributes <- tmpdf


  class(out) <- "Pairwise"
  return(out)
}

#' Pairwise Comparison (ordered numeric)
#'
#' @aliases pairwise_numeric
#'
#' @description Calculates pairwise comparison for ordered numeric data
#'
#' @usage
#' pairwise_numeric(outcome, discriminant=NULL, weight=NULL, attributes=NULL, data, simplify=FALSE, numThreads=1)
#'
#' @param times A character vector giving column names that contain outcomes.
#' @param discriminant A numeric vector giving the minimal difference for preference to be established for each outcome.
#' @param weight An optional character string giving the column name that contains the weight associated with each observation
#' @param attributes A character vector giving column names that contain attributes
#' @param data A data frame, containing the specified columns
#' @param numThreads The number of threads to use. If \code{numThreads<1}, the method uses that proportion of maximum available threads (rounded up).
#' @param simplify A boolean that currently does nothing but will eventually optionally simplify data down for faster computation
#'
#' @return An object of type \code{pairwise} containing the following:
#' \describe{
#'     \item{comparisons}{A square matrix containing all-to-all comparisons of outcomes}
#'     \item{group}{A vector giving the groups associated with each row}
#'     \item{weight}{The weight of each observation}
#'     \item{attributes}{A data frame containing attributes for each observation}
#'  }
#'
#' @details
#'
#' This function appears to perform best when \code{numThreads=1}.
#' This could be overhead outweighing benefits, or it could be something to do with linking to other functions? idk
#'
#'
#' @examples
#'
#'
#'
#' @export
pairwise_numeric <- function(outcome, discriminant=NULL, weight=NULL,  attributes=NULL, data, numThreads=NULL, simplify=FALSE)
{
  # tibbles break things
  data <- as.data.frame(data)


  outcome <- as.matrix(data[,outcome])

  if(is.null(discriminant))
  {
    # If no discriminant value was given, use half the smallest
    # nonzero difference on each outcome
    discriminant <- apply(outcome,2,function(x){
        differences <- diff(sort(x))
        differences <- differences[differences>0]
        if(length(differences)>0)
        {
          return(min(differences)/2)
        }
        else
        {
          # All values are the same, pick any number
          return(0.5)
        }
      })
  }

  if(is.null(numThreads))
  {
    numThreads <- RcppParallel::defaultNumThreads()
  }
  numThreads <- min(numThreads,RcppParallel::defaultNumThreads())
  if(numThreads < 1)
  {
    numThreads <- ceiling(numThreads * RcppParallel::defaultNumThreads())
  }
  RcppParallel::setThreadOptions(numThreads = numThreads)

  out <- {}
  out$comparisons <- gpct:::get_pw_numeric(outcome = outcome,discriminant = discriminant)
  out$comparisons

  if(is.null(weight))
  {
    out$weight <- rep(1,nrow(data))
  }
  else
  {
    out$weight <- data[,weight]
  }

  tmpdf <- {}
  for(cov in attributes){
    if(is.null(tmpdf)) {
        tmpdf <- as.data.frame(data[,cov])
      } else {
        tmpdf <- cbind(tmpdf,as.data.frame(data[,cov]))
      }
    colnames(tmpdf)[ncol(tmpdf)] <- cov
  }

  out$attributes <- tmpdf

  if(simplify) warning("simplification not included yet")

  class(out) <- "Pairwise"
  return(out)
}




#' Pairwise Comparison
#' @aliases pairwise
#' @description Collates variables into a \code{pairwise} object
#'
#' @usage
#' pairwise(comparisons, weight=NULL,  attributes=NULL)
#'
#' @param comparisons A matrix containing pairwise comparison information
#' @param weight An optional character string giving the column name that contains the weight associated with each observation
#' @param attributes A data frame containing attributes for each observation
#' @param simplify A boolean that currently does nothing but will eventually optionally simplify data down for faster computation
#'
#'
#' @return An object of type \code{pairwise} containing the input variables.
#'
#' @details
#'
#'
#' @examples
#'
#'
#'
#' @export
pairwise <- function(comparisons, weight=NULL, attributes=NULL)
{
  out <- {}
  out$comparisons <- comparisons

  if(is.null(weight))
  {
    out$weight <- rep(1,nrow(out$comparisons))
  }
  else
  {
    out$weight <- weight
  }

  out$attributes <- attributes

  class(out) <- "Pairwise"
  return(out)
}


#' Summarizing Pairwise Comparisons
#'
#' @description Produces a summary of a \code{Pairwise} object.
#'
#' @usage
#' summary(x, numThreads=NULL, includeTies=TRUE, countInconsistency=T, m=3)
#'
#' @param x an object of type \code{pairwise}.
#' @param numThreads The number of threads to use. If \code{numThreads<1}, the method uses that proportion of maximum available threads (rounded up).
#' @param includeTies This really needs to be deleted in the final version
#' @param countInconsistency
#' @param m The length of cycle for measuring inconsistency
#' @details
#'
#' Details here
#'
#' @examples
#'
#' @export
summary.Pairwise <- function(x, numThreads=NULL, includeTies=TRUE, countInconsistency=T, m=3)
{


  proportion_incomparable <- sum(is.na(x$comparisons[upper.tri(x$comparisons)]), na.rm=T)
  proportion_tied <- sum(x$comparisons[upper.tri(x$comparisons)]==0, na.rm=T)

  proportion_incomparable <- proportion_incomparable/length(upper.tri(x$comparisons))
  proportion_tied <- proportion_tied/length(upper.tri(x$comparisons))

  summary_edge <- c(winner=1-proportion_incomparable-proportion_tied,
                    tied=proportion_tied,
                    incomparable=proportion_incomparable)

  summary_weights <- summary(x$weight)

  summary_attributes <- summary(x$attributes)

  if(countInconsistency){

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

    inconsistency <- gpct:::countInconsistency(x$comparisons,x$weight,includeTies=includeTies, m=m)
  }
  else
  {
    inconsistency <- NULL
  }

  out <- list(edge = summary_edge,
              weights = summary_weights,
              attributes = summary_attributes,
              inconsistency = inconsistency
              )

  class(out) <- "summary.Pairwise"

  return(out)
}




