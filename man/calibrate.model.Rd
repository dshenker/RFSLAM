% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rf_functions_p1.R
\name{calibrate.model}
\alias{calibrate.model}
\title{Carlibrate RFSLAM model}
\usage{
calibrate.model(p.hat, rf.df.1, target_varname, time_varname)
}
\arguments{
\item{p.hat}{the event risk predictions}

\item{rf.df.1}{the dataframe used for modeling}

\item{target_varname}{name of the variable with the target}

\item{time_varname}{name of the variable with the CPIU count}
}
\description{
\code{calibrate.model} calibrates an RFSLAM model using the true event rates
}
