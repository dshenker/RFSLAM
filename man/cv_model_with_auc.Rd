% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rf_functions_p2.R
\name{cv_model_with_auc}
\alias{cv_model_with_auc}
\title{Perform Cross Validation for RFSLAM model}
\usage{
cv_model_with_auc(
  modeling_df,
  target,
  id_col,
  risk_time_col,
  patient_count_col = "int.n",
  time_col,
  n.folds,
  folds_stratifier,
  drop,
  ntree = 100,
  nodedepth = NULL,
  nsplit = 10
)
}
\arguments{
\item{modeling_df}{dataframe containing data for modeling}

\item{target}{name of the target variable column}

\item{id_col}{name of the column with the patient id's}

\item{risk_time_col}{name of the column with the risk time values}

\item{patient_count_col}{name of the column with the patient counts}

\item{n.folds}{the number of folds for cross validation}

\item{folds_stratifier}{the variable to use for creating the folds}

\item{drop}{the variables that need to be dropped before training the model}

\item{ntree}{number of trees for random forest}

\item{nodedepth}{the node depth for random forest}

\item{nsplit}{the nsplit parameter for random forest}
}
\value{
the average weighted auc value across all of the folds
}
\description{
\code{cv_model_with_auc} performs cross validation and returns the average auc value across the folds for the model
}
