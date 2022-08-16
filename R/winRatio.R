#' Win Ratio for Composite Endpoints
#'
#' @description
#' Calculates the (optionally stratified) Win Ratio for Composite Endpoints.
#'
#' @usage
#' winRatio <- function(x, group, stratum=NULL,
#'                     pool_method = c("unweighted","inversevariance","size"),
#'                     alpha=0.05, splitTies=FALSE, numThreads=NULL)
#'
#' @param x an object of type \code{pairwise}.
#' @param group a character string indicating the group variable
#' @param stratum A character vector indicating which covariates should be used for stratification
#' @param pool_method A character vector indicating how each stratum should be weighted
#' @param alpha The acceptable type 1 error used in the test.
#' @param splitTies A boolean specifying how ties should be treated.
#'                  Should be \code{TRUE} for Win Odds,
#'                  or \code{FALSE} for Win Ratio.
#' @param numThreads The number of threads to use. If \code{numThreads<1}, the method uses that proportion of maximum available threads (rounded up).
#' @return A list containing the following:
#' \describe{
#'     \item{pooled_lnodds}{The pooled log(odds).}
#'     \item{pooled_lnconf.int}{(1-\code{alpha})\% Confidence intervals for pooled log(odds).}
#'     \item{pooled_SElnodds}{Standard error of pooled log(odds).}
#'     \item{pooled_SElnnull}{Standard error of pooled log(odds) under the null hypothesis.}
#'     \item{pooled_p}{The p-value of the test of pooled log(odds) = 1.}
#'     \item{pooled_rel_statistic}{Statistic of test that strata odds are equal.}
#'     \item{pooled_rel_p}{p-value for test that strata odds are equal.}
#'     \item{relative_lnodds}{A matrix giving the log of the ratio of odds between strata (generalised relative risk ratio).}
#'     \item{relative_selnodds}{A matrix containing the standard error of the log(relative risk ratio).}
#'     \item{results}{A list containing a summary of each strata measure.}
#'     \item{param.record}{A list containing parameters used in the test.}
#'}
#'
#' @details
#'
#' Need to write this
#'
#' @references
#'
#' Dong, stratified win ratio
#' Dong, generalized analytic etc etc
#'
#' @export
winRatio <- function(x, group, stratum=NULL,
                     pool_method = c("unweighted","inversevariance","size"),
                     alpha=0.05, splitTies=FALSE, numThreads=NULL)
{

  if(is.null(numThreads)){
    numThreads <- Inf
  }
  numThreads <- min(numThreads,RcppParallel::defaultNumThreads())
  if(numThreads < 1){
    numThreads <- ceiling(numThreads * RcppParallel::defaultNumThreads())
  }
  RcppParallel::setThreadOptions(numThreads = numThreads)

  if(prod(x$weight)!=1){
    stop("All weights must be equal to 1.")
  }

  # Excract group and cast to factor
  group <- factor(x$attributes[,group])

  if(length(levels(group))!= 2) stop("group variable must contain 2 unique values")

  if(is.null(stratum)){
    strata <- rep(1,length(group))
  } else {
    if(length(stratum) > 2){
      strata <- apply(x$attributes[,stratum],1,paste,sep="_")
    } else {
      strata <- x$attributes[,stratum]
    }
  }

  # Convert NA values into magic number for C++
  x$comparisons[is.na(x$comparisons)] <- -99999

  lapply(unique(strata),function(i){

    subgroup <- which(strata==i)

    treatLevel <- levels(group)[2]

    winRatioVals <- rep(-1,8)
    gpct:::get_winRatioStats(winRatioVals,
                                pairwise = x$comparisons[subgroup,subgroup],
                                treatList =   which(group[subgroup]==treatLevel)-1,
                                controlList = which(group[subgroup]!=treatLevel)-1,
                                split = splitTies)

    if(min(table(group[subgroup]))<2) winRatioVals[3:8] <- NaN

    names(winRatioVals) <- c("nt","nc","st","sc","cov","stnull","scnull","covnull")

    # Win ratio point estimate
    WR <- winRatioVals["nt"]/winRatioVals["nc"]
    logWR <- log(WR)

    # Null win ratio is 1
    WRnull <- 1
    logWRnull <- 0

    varLogWR <- winRatioVals["st"]/winRatioVals["nt"]^2 +
      winRatioVals["sc"]/winRatioVals["nc"]^2 -
      2*winRatioVals["cov"]/(winRatioVals["nt"]* winRatioVals["nc"])


    n0 <- (winRatioVals["nt"] + winRatioVals["nc"])/2
    varLogWRnull <- (winRatioVals["stnull"] + winRatioVals["scnull"] - 2*winRatioVals["covnull"])/n0^2

    # 2 sided pVal
    pVal <- 2*(1- pnorm(abs(logWR/sqrt(varLogWRnull))))

    # Confidence interval
    confint <- logWR + c(lower=1,upper=-1)*qnorm(alpha/2)*sqrt(varLogWRnull)

    out <- c(
             logWR = unname(logWR),
             seLogWR = unname(sqrt(varLogWR)),
             seLogWRnull = unname(sqrt(varLogWRnull)),
             pVal = unname(pVal),
             confint = confint,
             winRatioVals,
             n=length(subgroup))

    return(out)

  }) -> winRatio_stratified

  if(is.null(stratum)){

    out <- winRatio_stratified

  } else {

    pooledVals <- do.call("rbind",lapply(pool_method,function(i){
      if(i=="inversevariance") {
        weights <- 1/sapply(winRatio_stratified,function(x){x["seLogWR"]^2})
        statistics <- sapply(winRatio_stratified,function(x){x["logWR"]})

        pooled_L <- sum(statistics*weights)/sum(weights)
        pooled_se <- 1/sum(weights)

        pooled_se <- sqrt(pooled_se)

        return(c(pooled_L=pooled_L,pooled_se=pooled_se,p=2*pnorm(abs(pooled_L/pooled_se),lower.tail = FALSE)))

      } else if(i %in% c("size","unweighted")) {

        if(i=="unweighted"){
          weights <- outer(rep(1,5),rep(1,length(winRatio_stratified)))
        } else {
          weights <- outer(rep(1,5), 1/sapply(winRatio_stratified,function(x){unname(x["n"])}))
        }

        # Variances use squared weights
        weights <- weights^outer(c(1,1,2,2,2), rep(1,ncol(weights)))

        statistics <- sapply(winRatio_stratified,function(x){x[c("nt","nc","st","sc","cov")]})

        pooled_statistics <- apply(weights*statistics,1,sum)

        pooled_L <- log(pooled_statistics["nt"]/pooled_statistics["nc"])
        pooled_se <- pooled_statistics["st"]/(pooled_statistics["nt"]^2) +
          pooled_statistics["sc"]/(pooled_statistics["nc"]^2) -
          2* pooled_statistics["cov"]/(pooled_statistics["nt"] * pooled_statistics["nc"])

        pooled_se <- sqrt(pooled_se)

        return(c(pooled_L=pooled_L,pooled_se=pooled_se,p=2*pnorm(abs(pooled_L/pooled_se),lower.tail = FALSE)))
      }
    }))

    rownames(pooledVals) <- pool_method

    names(winRatio_stratified) <- unique(strata)

    out <- list(strata = winRatio_stratified, pooled_statistics = pooledVals)

  }

  return(out)

}
