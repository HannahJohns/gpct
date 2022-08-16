#' simulate.survival
#'
#' @description Simulate a time-to-event trial using composite survival data
#'
#' @usage
#' simulate.survival(scales=NULL, rates=NULL, shapes=NULL, censor=NULL, heirarchy=NULL, hr, N, treatRatio = 1,
#'                   arrivalRate = NULL, maxt = NULL, missingAcrualRelax=NULL, minObsTime=0
#'                   )
#'
#' @param scales A vector (optionally named) containing the scale parameter for a series of weibull distributions
#' @param rates An alternative to scale (related by \code{rates=1/scales}).
#' @param shapes A vector containing the shape parameter for a series of weibull distributions
#' @param censor A matrix governing if certain events censor other events. See Details.
#' @param heirarchy An optional named vector giving the ordering of events (first event is most important).
#' @param hr A vector giving the hazard ratio for each event.
#' @param N Total sample size
#' @param treatRatio Ratio of treatment-to-control assignments planned.
#' @param arrivalRate Average number of arrivals per time unit. Assumed Markovian.
#' @param maxt Maximum time to run the trial.
#' @param missingAcrualRelax either "N" or "maxt". Governs behaviour if the trial does not meet recruitment targets.
#' @param minObsTime a numeric vector used if \code{missingAcrualRelax="maxt"}. Gives the minimum amount of time that all participants are followed for.
#' @param arrivalRate The average arrival rate for patients (number per time unit). See details for how to estimate this.
#'
#' @details
#'
#' Shape and Rate are given using the same parameterisation as pweibul, not simsurv.
#'
#' This function is designed to simulate a RCT with survival endpoint data,
#' including realistic censoring of patients who are recruited into the trial
#' at a later point. As such, it is neccesary to direct the method what should happen
#' if the planned recruitment process fails.
#'
#' This is controlled by the \code{missingAcrualRelax} parameter.
#' If \code{missingAcrualRelax="N"}, then the trial stops at \code{maxt} time units
#' regardless of current sample size.
#'
#' If \code{missingAcrualRelax="maxt"}, then the trial continues until the original sample size was met,
#' then continues observation for an extra \code{minObsTime} time units.
#'
#' In practice, this contingency should never be needed. The necessary value of
#' \code{arrivalRate} to accrue all patients on time can be calculated using the following command:
#'
#'  \code{uniroot(f = function(x){maxt-qgamma(tol, shape=N, rate=x)},interval=c(0.1,1000))}
#'
#'  Where \code{tol} is the probability of achieving the accrual target on time, and \code{maxt} and \code{N} are as described above.
#'
#'
#' @export
simulate.survival <- function(scales=NULL, rates=NULL, shapes=NULL, censor=NULL, heirarchy=NULL, hr, N, treatRatio = 1,
                                   arrivalRate = NULL, maxt = NULL, missingAcrualRelax=NULL, minObsTime=0,
                                   numThreads=NULL
                                  )
{


  if(!("simsurv" %in% installed.packages())) stop("simsurv package must be installed to use this function")

  testMethod <- "pwc" # Always run pairwise comparison

  if(!is.null(numThreads)){

    if(length(numThreads)>1)
    {
      if(is.null(names(numThreads))) stop("Must supply names for numThreads if multiple values are supplied")
      unspecified <- numThreads[names(numThreads)==""]
      if(length(unspecified) > 1) stop("Can only supply one unnamed numThreads value.")
    }
    else
    {
      unspecified <- numThreads
    }

    if(length(unspecified) == 1)
    {
      toAdd <- setdiff(testMethod,names(numThreads))
      added_threads <- rep(unspecified,length(toAdd))
      names(added_threads) <- toAdd

    numThreads <- c(numThreads,added_threads)
    numThreads <- numThreads[names(numThreads)!=""]

    }
    else
    {
        if(!(names(numThreads) %in% testMethod &
             testMethod %in% names(numThreads)))
        {
          if("pwc" %in% names(numThreads))
          {
            stop("All testMethods must have number of threads specified")
          }
          else
          {
            stop("All testMethods must have number of threads specified
                 (use \"pwc\" to specify number of threads for calculating the pairwise comparison matrix)")
          }
        }
    }

  }


  if(is.null(scales) + is.null(rates) != 1){
    stop("Must supply exactly one of scales or rates")
  }else if(is.null(scales)){
    scales <- 1/rates
  }

  # Recursive function for correcting simulated data which should be censored.
  # Provided by a DAG, where two colliding censoring events are chosen based on
  # whichever occurs first, and censoring is transitive/propagates through the graph

  if(!is.null(names(scales))){
    varNames <- names(scales)
  }  else {
    varNames <- paste("var",1:length(scales),sep="")
  }


  if(is.null(shapes)){
    shapes <- rep(1,length(varNames))
    names(shapes) <- varNames
  }

  if(is.null(censor)){
    censor <- matrix(0,length(varNames),length(varNames))
    rownames(censor) <- colnames(censor) <- varNames
  }

  if(is.null(heirarchy)){
    heirarchy <- varNames
  }

  if(!is.null(maxt) & !is.null(arrivalRate)){
    if(is.null(missingAcrualRelax)) stop("missingAcrualRelax needed")

    if(!(missingAcrualRelax %in% c("N","maxt")) )
    {
      stop('missingAcrualRelax must be one of "N" or "maxt"')
    }

    if(missingAcrualRelax == "N" & minObsTime > 0) stop("minObsTime is used for maxt relaxation only")

    if(missingAcrualRelax == "maxt" & minObsTime == 0) warning("minObsTime should be greater than zero")

  }

  if(is.null(missingAcrualRelax)){
    # We need a non-null value for this for later checks
    missingAcrualRelax <- "none"
  }


  # Generate random data

  # We will handle censoring later based on input details, assume all events are observed for simsurv

  df <- data.frame(id = 1:N)

  randTable <- c(rep(1,round(N*treatRatio/(treatRatio+1))),rep(0,round(N/(treatRatio+1))))

  # Catch rounding errors
  while(length(randTable)-nrow(df)>=2){
    # Remove a 0 and a 1 to keep even
    randTable <- randTable[2:(length(randTable)-1)]
  }
  randTable <- randTable[order(runif(length(randTable)))]

  if(length(randTable)-nrow(df)==1){
    randTable <- randTable[-1]
  }

  df$trt <- randTable

  for( i in 1:length(scales)){

    eventData <- simsurv::simsurv(dist = "weibull",
                         lambdas = unname(scales[i] ^(-shapes[i])),
                         gammas = unname(shapes[i]),
                         betas = c(trt = unname(log(unlist(hr[i])))),
                         x=df)
    eventData <- eventData[,2:3]

    colnames(eventData) <- paste(varNames[i],colnames(eventData),sep="_")

    eventData

    df <- cbind(df,eventData)
  }

  # Handle censoring driven by competing risks
  correctedTimes <- matrix(0,nrow(df),length(varNames))
  gpct:::correctTimes(correctedTimes,as.matrix(df[paste(varNames,"eventtime",sep="_")]),censor)


  if(correctedTimes[1,1] == -999) error("Cycle detected in censor path!")

  df[,paste(varNames,"status",sep="_")][correctedTimes < as.matrix(df[paste(varNames,"eventtime",sep="_")])] <- 0
  df[,paste(varNames,"eventtime",sep="_")] <- correctedTimes

  if(!is.null(arrivalRate)){

    df <- cbind(df[,1:2],data.frame(arrivalTime = cumsum(rexp(nrow(df),rate=arrivalRate))),df[,-(1:2)])

    flagOverrun <- sum(df$arrivalTime > maxt)>0

    if(missingAcrualRelax == "N" & flagOverrun)
    {
      df <- df[which(df$arrivalTime < maxt),]
    }

    if(missingAcrualRelax == "maxt" & flagOverrun)
    {
      maxt <- max(df$arrivalTime) + minObsTime
    }

    if(!is.null(maxt))
    {
      for( i in 1:length(scales)){

        # Did event occur after maxt?
        endOfTrialCensored <- c(maxt < df$arrivalTime + df[paste(varNames[i],"eventtime",sep="_")])

        df[which(endOfTrialCensored),
           paste(varNames[i],"status",sep="_")] <- 0

        df[which(endOfTrialCensored),
           paste(varNames[i],"eventtime",sep="_")] <- maxt - df[which(endOfTrialCensored),
                                                                "arrivalTime"]
      }
    }
  } else {
    if(!is.null(maxt))
    {
      for( i in 1:length(scales)){

        df[which(maxt < df[,paste(varNames[i],"eventtime",sep="_")]),
           paste(varNames[i],"status",sep="_")] <- 0

        df[which(maxt < df[,paste(varNames[i],"eventtime",sep="_")]),
           paste(varNames[i],"eventtime",sep="_")] <- maxt

      }
    }
  }


  pwc <- gpct:::pairwise_surv(times = paste(varNames[i],"eventtime",sep="_"),
                              statuses = paste(varNames[i],"status",sep="_"),
                              attributes = c("id","trt","arrivalTime"),
                              data = df)

   out <- list(data=df,pwc=pwc)

   return(out)
}




