% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4.3.2_nonparametric.R
\name{roc.gaussian}
\alias{roc.gaussian}
\title{function to draw smooth ROC curve using gaussian kernel}
\usage{
roc.gaussian(data0, data1)
}
\arguments{
\item{data0}{continuous test results of nondiseased subjects}

\item{data1}{continuous test results of diseased subjects}
}
\value{
a smooth ROC curve
#example
#data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#roc.gaussian(data0,data1)
}
\description{
function to draw smooth ROC curve using gaussian kernel
}
