% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4.2.5.R
\name{roc.CI.ordinal}
\alias{roc.CI.ordinal}
\title{Function to calculate the full area and partial area and their variances}
\usage{
roc.CI.ordinal(data0, data1, e1, e2)
}
\arguments{
\item{data0}{ordinal test results of nondiseased subjects}

\item{data1}{ordinal test results of diseased subjects}

\item{e1}{lower bound of FPR for partial area}

\item{e2}{upper bound of FPR for partial area}
}
\value{
estimation of full area,partial area and their variances
#example
#data0<-c(38,25,15,19,4)
#data1<-c(1,2,3,14,42)
#roc.CI.ordinal(data0,data1,0,0.2)
}
\description{
Function to calculate the full area and partial area and their variances
}
