% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{atipicalitySpectra}
\alias{atipicalitySpectra}
\title{Identification of potentially low-quality raw mass spectra}
\usage{
atipicalitySpectra(s, labels = NULL, BPPARAM = bpparam())
}
\arguments{
\item{s}{Spectra object}

\item{labels}{factor to do groupwise calculation of A-score limits for flagging
spectra.}

\item{BPPARAM}{Parallel computing configuration from \link{BiocParallel}. Default is
is \code{\link[BiocParallel:register]{BiocParallel::bpparam()}}}
}
\value{
A list with \code{\link[Spectra:Spectra]{Spectra::Spectra()}} object with \code{QCflag} spectrum
variable added and the atipicality limits for flagging.
}
\description{
Uses \link{atipicality_spectra}, which implements \code{\link[MALDIrppa:screenSpectra]{MALDIrppa::screenSpectra()}} on
peak matrix. It can calculate calculate lower and upper A limits by label.
}
