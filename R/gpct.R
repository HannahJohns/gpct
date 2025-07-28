#' Twice-Generalized Odds Ratios
#'
#' @description Calculates Generalised Pairwise Comparisons for Trend
#'
#' @usage
#' gpct(x, explanatory_var, stratum=NULL,
#'         alpha=0.05, splitTies=FALSE)
#'
#' @param x an object of type \code{pairwise}.
#' @param explanatory_var a character string indicating the explanatory variable
#' @param alpha The acceptable type 1 error used in the test.
#' @param splitTies A boolean specifying how ties should be treated.
#'                  Should be \code{TRUE} for WMW Odds,
#'                  or \code{FALSE} for Agresti's GenOR.
#' @param numThreads The number of threads to use. If \code{numThreads<1}, the method uses that proportion of maximum available threads (rounded up).
#' @param nnt A boolean.
#'            If \code{TRUE}, then print number needed to treat in addition to generalised odds ratios.
#' @param verbose A boolean.
#'                If \code{TRUE}, then print both pooled odds and relative risk ratio matrices
#'                regardless of result of statistical test.
#' @param upper A boolean specifying if the upper triangle
#'              of relative risk ratios should be printed.
#'              If \code{FALSE}, lower triangle is used instead.
#'
#' @return A list with class "\code{gpct}" containing the following:
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
#' For details of Generalized Odds Ratios, see \code{\link{genodds}}.
#'
#' Generalized Generalized Odds Ratios is a generalisation of Agresti's
#' generalized odds ratios (GenOR), capable of handling arbitrary outcomes.
#'
#' This function implements Generalized Generalised Odds Ratios for
#'
#'
#' @references
#'
#'
#' Our paper will go here
#'
#' @export
gpct <- function(x, explanatory_var, stratum=NULL, alpha=0.05, splitTies=FALSE, numThreads=NULL, senull_method="none")
{

  if(is.null(numThreads)){
    numThreads <- Inf
  }
  numThreads <- min(numThreads,RcppParallel::defaultNumThreads())
  if(numThreads < 1){
    numThreads <- ceiling(numThreads * RcppParallel::defaultNumThreads())
  }
  RcppParallel::setThreadOptions(numThreads = numThreads)

  weights <- x$weight
  groupVector <- as.numeric(x$attributes[,explanatory_var])

  tol <- min(diff(sort(unique(groupVector))))/2

  if(length(unique(groupVector))<2){
    stop("Need at least two unique groups to continue")
  }

  # Convert NA values into magic number for C++
  x$comparisons[is.na(x$comparisons)] <- -99999

  if(is.null(stratum)){
    strata <- rep(1,length(groupVector))
  } else {
    if(length(stratum) > 2){
      strata <- apply(x$attributes[,stratum],1,paste,sep="_")
    } else {
      strata <- x$attributes[,stratum]
    }
  }

  id <- 1:nrow(x$comparisons)

  lapply(unique(strata),function(i){

    subgroup <- which(strata==i)

    Rsdt <- gpct:::getRsdt(pairwise = x$comparisons[subgroup,subgroup],
                              group = groupVector[subgroup],
                              weight = weights[subgroup],
                              tol = tol)

    if(splitTies){
      Rsdt[,1] <- Rsdt[,1] + 0.5*Rsdt[,3]
      Rsdt[,2] <- Rsdt[,2] + 0.5*Rsdt[,3]
    }

    # Get Odds/SE/etc

    Pc <- sum(weights[subgroup] * Rsdt[,1])/sum(weights[subgroup])
    Pd <- sum(weights[subgroup] * Rsdt[,2])/sum(weights[subgroup])

    OR <- Pc/Pd
    # SE <- 2/Pd * sqrt(sum( weights * (OR*Rsdt[,2] - Rsdt[,1] )^2 )/sum(weights)^2)
    SE <- 2/(Pd) * sqrt( sum( (weights[subgroup]/sum(weights[subgroup])) * (OR*(Rsdt[,2]) - (Rsdt[,1]) )^2 )/sum(weights[subgroup]))


    probwin_Y <- gpct:::winProb(x$comparisons[subgroup,subgroup],weights[subgroup],population = F)
    probwin_Y <- as.data.frame(probwin_Y)
    colnames(probwin_Y) <- c("prob_win","prob_loss","prob_tie")

    # Log everything

    logOR <- log(OR)
    # logORNull <- log(ORNull)

    logSE <- SE/OR
    # logSENull <- SENull/ORNull
    # logSENull_fast <- SENull_fast/ORNull

    # 2 sided pVal
    pVal <- 2*(1- pnorm(abs(logOR/logSE)))
    # pVal_fast <- 2*(1- pnorm(abs(logOR/logSENull_fast)))

    # Confidence interval
    confint <- logOR + c(lower=1,upper=-1)*qnorm(alpha/2)*logSE
    # confint_fast <- logOR + c(lower=1,upper=-1)*qnorm(alpha/2)*logSENull_fast


    data <- {}
    for(cov in colnames(x$attributes)){
      if(is.null(data)) {
        data <- as.data.frame(x$attributes[subgroup,cov])
      } else {
        data <- cbind(data,as.data.frame(x$attributes[subgroup,cov]))
      }
      colnames(data)[ncol(data)] <- cov
    }

    out <- list(logOR = logOR,
                logSE = logSE,
                pVal = pVal,
                confint = confint,

                data=data,

                probwin_Y=probwin_Y
    )

    return(out)

  }) -> gpct_stratified


  if(is.null(stratum)){

    out <- gpct_stratified[[1]]

  } else {


    weights <- 1/sapply(gpct_stratified,function(x){x$logSE^2})
    statistics <- sapply(gpct_stratified,function(x){x$logOR})

    pooled_L <- sum(statistics*weights)/sum(weights)
    pooled_se <- 1/sum(weights)
    pooled_se <- sqrt(pooled_se)

    pooled_homogeneity_stat <- sum(weights*(statistics - pooled_L)^2)

    pooled_homogeneity_p <- pchisq(pooled_homogeneity_stat, length(weights)-1, lower.tail = FALSE)

    pooledVals <- list(pooled_L=pooled_L,
                       pooled_se=pooled_se,
                       p=2*pnorm(abs(pooled_L/pooled_se),lower.tail = FALSE),
                       pooled_homogeneity_stat = pooled_homogeneity_stat,
                       pooled_homogeneity_p = pooled_homogeneity_p
                       )

    names(gpct_stratified) <- unique(strata)

    out <- list(strata = gpct_stratified, pooled_statistics = pooledVals)

  }

  class(out) <- c("gpct",class(out))

  return(out)
}
