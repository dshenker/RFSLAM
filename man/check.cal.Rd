% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rf_functions_p1.R
\name{check.cal}
\alias{check.cal}
\title{Check the Calibration of an RFSLAM model}
\usage{
check.cal(predicted.rate, actual.outcomes, rt)
}
\arguments{
\item{predicted.rate}{the event risk predictions}

\item{actual.outcomes}{the true events}

\item{rt}{the risk time values for each patient}
}
\description{
\code{check.cal} checks how well calibrated an RFSLAM model is using the predicted event rates and the true events
}
