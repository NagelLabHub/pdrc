% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comet_looping.R
\name{extract_model_scores}
\alias{extract_model_scores}
\title{extract_model_scores}
\usage{
extract_model_scores(results_list)
}
\arguments{
\item{results_list}{A named list of results (from loop_comet_data(..., method = "adaptive")).}
}
\value{
A data frame combining model_scores for all samples, with \code{Sample} as the first column.
}
\description{
Extract model_scores tables from all samples in an adaptive fitting result list.
}
