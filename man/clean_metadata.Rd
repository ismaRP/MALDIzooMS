% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/platemap.R
\name{clean_metadata}
\alias{clean_metadata}
\title{Clean spectra metadata}
\usage{
clean_metadata(metadata, folder)
}
\arguments{
\item{metadata}{A data.frame with a replicate spectra per row. It must contain
at least columns sample_name and replicate}

\item{folder}{Folder where spectra data is stored}
}
\value{
data.frame metadata with extra columns \code{file},
\code{spectra_name} and \code{n_replicates}
}
\description{
Remove spectra in metadata that don't exist in data
}
