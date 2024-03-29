% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ROKET.R
\name{kernTEST}
\alias{kernTEST}
\title{kernTEST}
\usage{
kernTEST(
  RESI = NULL,
  KK,
  YY = NULL,
  XX = NULL,
  OMNI,
  nPERMS = 1e+05,
  ncores = 1
)
}
\arguments{
\item{RESI}{A numeric vector of null model residuals
\code{names(RESI)} must be set to maintain sample ordering
for survival regression, otherwise set \code{RESI} to \code{NULL}.}

\item{KK}{An array containing double-centered positive semi-definite
kernel matrices. Refer to \code{MiRKAT::D2K()} for transforming 
distance matrices to kernel matrices. The \code{dimnames(KK)[[1]]} and 
\code{dimnames(KK)[[2]]} must match \code{names(RESI)}.
Also set dimnames(KK)[[3]] to keep track of each kernel matrix.}

\item{YY}{A numeric vector of continuous outcomes to be fitted
in a linear model. Defaults to NULL for survival model.}

\item{XX}{A numeric data matrix with first column for intercept,
a column of ones.}

\item{OMNI}{A matrix of zeros and ones. Each column corresponds to a
distance matrix while each row corresponds to an omnibus test. Set
\code{rownames(OMNI)} for labeling outputted p-values and 
\code{colnames(OMNI)} which should match \code{dimnames(KK)[[3]]}.}

\item{nPERMS}{A positive integer to specify the number of
permutation-based p-value calculation}

\item{ncores}{A positive integer for the number of cores/threads
to reduce computational runtime when running for loops}
}
\description{
kernTEST
}
