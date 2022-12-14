% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pairwise_definitions.R
\name{pairwise_surv}
\alias{pairwise_surv}
\title{Pairwise Comparison (composite survival)}
\usage{
pairwise_surv(times, statuses, weight=NULL, attributes=NULL, data, simplify=FALSE, numThreads=1)
}
\arguments{
\item{times}{A character vector giving column names that contain event times.}

\item{statuses}{A character vector giving column names that contain event statuses (1 for event occurred, 0 for censored)}

\item{weight}{An optional character string giving the column name that contains the weight associated with each observation}

\item{attributes}{A character vector giving column names that contain attributes of each patient}

\item{data}{A data frame, containing the specified columns}

\item{numThreads}{The number of threads to use. If \code{numThreads<1}, the method uses that proportion of maximum available threads (rounded up).}

\item{simplify}{A boolean that currently does nothing but will eventually optionally simplify data down for faster computation}
}
\value{
An object of type \code{pairwise} containing the following:
\describe{
    \item{comparisons}{A square matrix containing all-to-all comparisons of outcomes}
    \item{group}{A vector giving the groups associated with each row}
    \item{weight}{The weight of each observation}
    \item{attributes}{A data frame containing attributes for each observation}
 }
}
\description{
Calculates pairwise comparison for ordered survival data
}
\details{
This function appears to perform best when \code{numThreads=1}.
This could be overhead outweighing benefits, or it could be something to do with linking to other functions? idk
}
