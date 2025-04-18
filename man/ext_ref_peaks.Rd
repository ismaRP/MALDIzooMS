% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align_peaks.R
\name{ext_ref_peaks}
\alias{ext_ref_peaks}
\title{Build reference Spectrum from external sequences}
\usage{
ext_ref_peaks(
  sequences,
  mc.cores = 4L,
  gpo_only = F,
  ret.object = c("df", "MassPeaks", "Spectra"),
  non_deam = FALSE
)
}
\arguments{
\item{sequences}{Dataframe with collumns: species, chain, id, sequence}

\item{mc.cores}{Numbre of CPUs for parallel processing}

\item{gpo_only}{Logical. If \code{TRUE}, only the second proline in GPP patterns can be
hydroxilated, becoming GPO. If \code{FALSE}, both can be hyroxylated.}

\item{ret.object}{Type of object to return. It can be:
\itemize{
\item 'df': for a data frame,
\item 'MassPeaks': for a \link[MALDIquant:MassPeaks-class]{MALDIquant::MassPeaks} object
\item 'Spectra': for a \link[Spectra:Spectra]{Spectra::Spectra} object
}}

\item{non_deam}{logical, whether to include deamidated peptides for those
with N or Q aminoacids}
}
\value{
Reference peak list. The type depends on \code{ret.object}:
\link[Spectra:Spectra]{Spectra::Spectra}, \link[MALDIquant:MassPeaks-class]{MALDIquant::MassPeaks} or a data frame
}
\description{
Build reference Spectrum from external sequences
}
