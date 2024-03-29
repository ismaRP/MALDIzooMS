% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{peaksLocalBG}
\alias{peaksLocalBG}
\title{peaksLocalBG}
\usage{
peaksLocalBG(s, halfWindowSize, mass_range, bg_cutoff, l_cutoff, SNR = 0)
}
\arguments{
\item{s}{\code{\link[MALDIquant]{MassSpectra}} object}

\item{halfWindowSize}{Half-window size parameter for local maximum detection. Passed to \code{\link[MALDIquant]{detectPeaks}}}

\item{mass_range}{Mass window to both sides of a peak to be considered for backgroun modelling}

\item{bg_cutoff}{The peaks within the mass range with intensity below the \code{bg_cutoff} quantile
are considered for background modelling. \code{bg_cutoff=1} keeps all peaks
and \code{bg_cutoff=0.5} would only keep the bottom half.}

\item{l_cutoff}{Likelihood threshold or p-value. Peaks with a probability of being modelled as
background noise higher than this are filtered out.}

\item{SNR}{Signal to noise threshold for peak detection. Passed to \code{\link[MALDIquant]{detectPeaks}}}
}
\description{
peaksLocalBG
}
