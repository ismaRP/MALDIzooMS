% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{smooth}
\alias{smooth}
\title{Intensity smoothing spectrapeak matrix}
\usage{
smooth(
  x,
  method = c("MovingAverage", "WeightedMovingAverage", "SavitzkyGolay"),
  hws = 4L,
  k = 3L,
  int_index = 2,
  in_place = FALSE,
  ...
)
}
\arguments{
\item{x}{Peak matrix}

\item{method}{Smoothing method. One of \code{'MovingAverage'}, \code{'WeightedMovingAverage'}
or \code{'SavitzkyGolay'}}

\item{hws}{Half window size}

\item{k}{For Savitzky-Golay filter, this is the polynomial order for the coefficients}

\item{int_index}{Column with the intensity to smooth}

\item{in_place}{Whether smoothing is in place, or a new column with the smoothened
intensity is added}

\item{...}{Parameters passed to other methods. Not in use, required by \link[Spectra:Spectra]{Spectra::addProcessing,Spectra-method}.}
}
\value{
Peak matrix
}
\description{
Intensity smoothing spectrapeak matrix
}
