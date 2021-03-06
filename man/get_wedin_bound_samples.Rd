% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rajive_helpfunctions.R
\name{get_wedin_bound_samples}
\alias{get_wedin_bound_samples}
\title{Gets the wedin bounds}
\usage{
get_wedin_bound_samples(X, SVD, signal_rank, num_samples = 1000)
}
\arguments{
\item{X}{Matrix. The data matrix.}

\item{SVD}{List. The SVD decomposition of the matrix. List with entries 'u', 'd', and 'v'from the svd function.}

\item{signal_rank}{Integer.}

\item{num_samples}{Integer. Number of vectors selected for resampling procedure.}
}
\description{
Gets the wedin bounds
}
