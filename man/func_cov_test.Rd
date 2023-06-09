% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_cov_test.R
\name{func_cov_test}
\alias{func_cov_test}
\title{Two sample test for eigendecompositions of functional data}
\usage{
func_cov_test(X, Y, paired = FALSE, n_perm = 1000, pve = 0.99)
}
\arguments{
\item{X}{A matrix of functions}

\item{Y}{A matrix of functions}

\item{paired}{A logical indidating wether you want a paired test}

\item{n_perm}{Number of permutations used to calculate p-value for paired test}

\item{pve}{proportion varaince explained; a value between 0-1}
}
\description{
Function to compare the eigendecomposition of two sets of functional data.
}
