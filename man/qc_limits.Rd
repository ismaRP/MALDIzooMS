% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{qc_limits}
\alias{qc_limits}
\title{Calculate outlier limits for Atipicality score using \link[robustbase:adjboxStats]{robustbase::adjboxStats}}
\usage{
qc_limits(A)
}
\arguments{
\item{A}{vector of A scores from QC}
}
\value{
Vector with top and bottom limits
}
\description{
Calculate outlier limits for Atipicality score using \link[robustbase:adjboxStats]{robustbase::adjboxStats}
}
