#' Multi-group Win Proportion
#'
#' @description Calculates the win proportion for multiple groups
#'
#' @usage
#' winprop(x, splitTies=FALSE)
#'
#' @param x an object of type \code{pairwise}.
#' @param alpha The acceptable type 1 error used in the test.
#' @param splitTies A boolean specifying how ties should be treated.
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
#' Our paper will go here
#'
#' @export


winProportion <- function(x, splitTies=FALSE, numThreads=NULL)
{

  if(is.numeric(x$group)) x$group <- factor(x$group)
  groupVector <- as.numeric(x$group)
  nGroups <- length(unique(groupVector))

  tallyWeights <- tapply(x$weight,groupVector,sum)
  sumWeights <- sum(tallyWeights)
  prodWeights <- prod(tallyWeights)

  p <- tallyWeights/sumWeights

  if(length(unique(groupVector))<2)
  {
    stop("Need at least two unique groups to continue")
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


  winProp <- genodds:::get_multigroup_winproportion(pairwise = x$comparisons,
                                                    group = groupVector,
                                                    nGroups = nGroups,
                                                    weight = x$weight)

  winProp$tally_dropped <- apply(winProp$tally_dropped,2,sum)
  winProp$tally_split <- apply(winProp$tally_split,2,sum)
  winProp$tally_nWinners <-  apply(winProp$tally_nWinners,2,sum)

  u_dropped <- winProp$tally_dropped/prodWeights
  u_split <- winProp$tally_split/prodWeights


  p_nWin <- winProp$tally_nWinners/sum(winProp$tally_nWinners)



  # This breaks on nGroups>2


  # u_split_null <-  1/nGroups
  # u_dropped_null <- 1/nGroups
  #
  # V_dropped <- sumWeights * (2*nGroups-1) * sum(p*(u_dropped - u_dropped_null)^2 - sum(p*(u_dropped - u_dropped_null))^2)
  # V_split <- sumWeights * (2*nGroups-1) * sum(p*(u_split - u_split_null)^2 - sum(p*(u_split - u_split_null))^2)
  #
  #
  # V_dropped
  # V_split

  # Section 6.2.5 of
  # Lee, A. J. (2019). U-Statistics: Theory and Practice. United States: CRC Press.

  Z_split <-  u_split - (1-p_nWin[1])/nGroups
  Z_dropped <- u_dropped - (p_nWin[2])/nGroups

  Z_split_mean <-  sum(p*Z_split)
  Z_dropped_mean <-  sum(p*Z_dropped)

  V_dropped <- sumWeights * (2*nGroups-1) * sum(p*(Z_dropped - Z_dropped_mean)^2)
  V_split <- sumWeights * (2*nGroups-1) * sum(p*(Z_split - Z_split_mean)^2)

  pval_V_dropped <- pchisq(V_dropped,nGroups-1,lower.tail = FALSE)
  pval_V_split <- pchisq(V_split,nGroups-1,lower.tail = FALSE)

  out <- list(drop=c(prop=u_dropped,Z=Z_dropped,V=V_dropped,p.value=pval_V_dropped),
              split=c(prop=u_split,Z=Z_split,V=V_split,p.value=pval_V_split),
              prop_group = p,
              prop_nWin = p_nWin
  )
  return(out)
}





