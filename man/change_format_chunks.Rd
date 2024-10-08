% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_spectra.R
\name{change_format_chunks}
\alias{change_format_chunks}
\title{Change MALDI-TOF data format}
\usage{
change_format_chunks(
  spectra_names,
  indir,
  in_ext,
  readf,
  outpath,
  writef,
  sep = "\\t",
  mc.cores = 4,
  nchunks = 80,
  verbose = NULL
)
}
\arguments{
\item{spectra_names}{File names to transform (without extension)}

\item{indir}{Path to data folder}

\item{in_ext}{Extension of input files}

\item{readf}{Input format reading function. One of \code{'fread'}, \code{'table'}
or \code{'mzML'}. \code{'fread'} uses [data.table::fread()],
\code{'table'} uses [utils::read.table()] and \code{mzML} uses [mzR::openMSfile()].}

\item{outpath}{Destination path}

\item{writef}{Output format writing function. one of 'tab' or 'mzML'.
Determines the destination format of the files.}

\item{sep}{Separator for tsv input files. Default is \code{"\t"}}

\item{mc.cores}{Number of cores to use}

\item{nchunks}{Number of chunks to split the files into}

\item{verbose}{Print progress bar}
}
\description{
Read MALDI data in tsv (tab) or mzML format in chunks and export in tsv (tab) or mzML.
It is a wrapper interface around reading and writting functions.
}
