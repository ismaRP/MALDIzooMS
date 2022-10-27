
#' Title
#'
#' @param x
#' @param n
#'
#' @return n list of vectos
#' @export
#'
#' @examples
chunks = function(x,n){
  split(x, cut(seq_along(x), n, labels = FALSE))
}



#' prepf
#' Pre-load spectra preprocessing function with arguments using a closure
#'
#' @param FUN
#' Function to be loaded with argument in \code{...}
#' @param ...
#' Arguments passed to \code{FUN}
#' @return
#' A closure with the preprocessing function \code{FUN} and preloaded arguments
#' @export
#'
#' @examples
prepFun = function(FUN, ...){
  function(s) {
    FUN(s, ...)
  }
}

#' peptidePseudoClusters
#'
#' Extracts pseudo envelopes from given peptide masses
#'
#' @param x \code{\link[MALDIquant]{MassPeaks}} object
#' @param masses
#' 1D Vector of isotopic masses to search.
#' See examples for how to generate from a vector of single peptide masses
#' @param tol
#' @param n_isopeaks
#' @param min_isopeaks
#'
#' @return data.frame with mass, intensity and s2n values for pseudo-isotopic
#'         clusters
#' @export
#' @importFrom MALDIquant match.closest
#'
#' @examples
#' masses = c(100, 200, 300)
#' n_isopeaks = 5
#' masses = matrix(masses, nrow = n_isopeaks, ncol = length(masses), byrow = TRUE)
#' d = 1.00235
#' masses = masses + (d * 0L:(n_isopeaks - 1L))
#' masses = sort(masses)
#' print(masses)



peptidePseudoClusters = function(x, masses, tol, n_isopeaks, min_isopeaks){
  n_pepts = length(masses)/n_isopeaks

  sel_masses = rep(NA, length(masses))
  sel_int = rep(NA, length(masses))
  sel_snr = rep(NA, length(masses))
  sel_complete = rep(F, length(masses))
  idx = match.closest(x@mass, masses, tolerance = tol*masses)

  sel_complete[idx[!is.na(idx)]] = T

  sel_complete = split(sel_complete,
                       cut(seq_along(sel_complete), n_pepts, labels = FALSE))

  sel_complete = lapply(
    sel_complete,
    function(x){
      s = cumsum(x)
      first_missing = (s[1] == 0) & (s[min_isopeaks + 1] == min_isopeaks)
      last_missing = s[min_isopeaks] == min_isopeaks
      if (!(first_missing | last_missing)){
        return(rep(F, n_isopeaks))
      } else {
        return(x)
      }
    }
  )

  sel_complete = Reduce(c, sel_complete)
  sel_masses[idx[!is.na(idx)]] = x@mass[!is.na(idx)]
  sel_masses[!sel_complete] = NA
  sel_int[idx[!is.na(idx)]] = x@intensity[!is.na(idx)]
  sel_int[!sel_complete] = NA
  sel_snr[idx[!is.na(idx)]] = x@snr[!is.na(idx)]
  sel_snr[!sel_complete] = NA
  return(data.frame(mass=sel_masses, intensity=sel_int, snr=sel_snr))
}



pseudoClusters = function(x, ...){
  prep_args = list(...)
  masses = prep_args$masses
  tol = prep_args$tol
  n_isopeaks = prep_args$n_isopeaks

  # TODO: finish from MALDIquant
}



#' modelLocalBG
#'
#' @param p
#' Matrix, mass-intensity pair
#' @param m
#' Full matrix of mass-intensity pairs
#' @param mass_range
#' Mass window to both sides of \code{p} to be considered for backgroun modelling
#' @param bg_cutoff
#' The peaks within the mass range with intensity below the \code{bg_cutoff} quantile
#' are considered for background modelling. \code{bg_cutoff=1} keeps all peaks
#' and \code{bg_cutoff=0.5} would only keep the bottom half.
#' @param l_cutoff
#' Likelihood threshold or p-value. Peaks with a probability of being modelled as
#' background noise higher than this are filtered out.
#' @return
#' @importFrom MASS fitdistr
#' @importFrom stats dnorm p.adjust pchisq pnorm quantile
#' @importFrom graphics hist
#' @export
#'
#' @examples
modelLocalBG = function(p, full_m, mass_range, bg_cutoff){

  # Get masses around
  bgmask = (full_m[,1] < (p[1] + mass_range)) & (full_m[,1] > (p[1] - mass_range))
  full_m = full_m[bgmask, ]
  # Remove top bg_cutoff
  int_cutoff = quantile(full_m[,2], bg_cutoff)
  bgmask = full_m[,2] < int_cutoff

  if(sum(bgmask) < 15) {
    return(matrix(c(0, NA, NA), nrow=1))
  }

  full_m = full_m[bgmask, ]

  bg = full_m[,2]
  fit = fitdistr(bg, 'normal')

  h = hist(bg, plot=F, breaks = 20)
  binwidth = h$mids[2] - h$mids[1]
  co = h$counts
  ce = dnorm(h$mids, mean= fit$estimate[1], sd = fit$estimate[2]) * length(bg) * binwidth

  chsq = sum((co-ce)^2/ce)
  pval = 1 - pchisq(chsq, length(bg)-1)


  l = pnorm(p[2], mean = fit$estimate[1], sd = fit$estimate[2], lower.tail = F)

  return(matrix(c(l, fit$loglik, pval), nrow=1))

}


#' peaksLocalBG
#'
#' @param s
#' \code{\link[MALDIquant]{MassSpectra}} object
#' @param halfWindowSize
#' Half-window size parameter for local maximum detection. Passed to \code{\link[MALDIquant]{detectPeaks}}
#' @param SNR
#' Signal to noise threshold for peak detection. Passed to \code{\link[MALDIquant]{detectPeaks}}
#' @param mass_range
#' Mass window to both sides of a peak to be considered for backgroun modelling
#' @param bg_cutoff
#' The peaks within the mass range with intensity below the \code{bg_cutoff} quantile
#' are considered for background modelling. \code{bg_cutoff=1} keeps all peaks
#' and \code{bg_cutoff=0.5} would only keep the bottom half.
#' @param l_cutoff
#' Likelihood threshold or p-value. Peaks with a probability of being modelled as
#' background noise higher than this are filtered out.
#' @return
#' @importFrom MALDIquant detectPeaks as.matrix
#' @export
#'
#' @examples
peaksLocalBG = function(s, halfWindowSize, mass_range, bg_cutoff, l_cutoff, SNR=0){

  s = detectPeaks(s, method="SuperSmoother", SNR=SNR,
                  halfWindowSize=halfWindowSize)

  total_peaks = length(s)
  m = as.matrix(s)
  snr_values = s@snr

  l_values = apply(
    X=m, MARGIN=1,
    FUN=modelLocalBG,
    full_m=m, mass_range=mass_range, bg_cutoff=bg_cutoff)
  l_values = t(l_values)

  adj_pval = p.adjust(l_values[,1], 'holm')

  l_values = cbind(l_values, adj_pval)

  colnames(l_values) = c("bg_lik", "fit_loglik", "chsq_pval", "bg_adjlik")

  peaks_mask = l_values[,1] < l_cutoff
  m = m[peaks_mask, ]
  m = cbind(m, snr_values[peaks_mask])
  non_bg_peaks = nrow(m)
  frac_rem = non_bg_peaks/total_peaks

  p = createMassPeaks(
    m[,1], m[,2], m[,3],
    metaData = s@metaData)
  p@metaData$fitting = l_values[peaks_mask, ]
  p@metaData$prepQC$frac_rem = frac_rem
  return(p)
}




#' preprocessData
#'
#' Performs smoothening, baseline removal and peak detection on MALDI samples.
#'
#' @param indir Folder containing spectra.
#' @param readf
#' A string value. choose function to use to read spectra.
#' Currently restricted to one of "fread", "table" or "mzml"
#' @param sep
#' Separator character for input files
#' @param nchunks
#' Number of chunks all the spectra should be divided into for reading and
#' processing. If all spectra is loaded and processed at once, i.e. nchunks=1,
#' it can overload RAM memory. If nchunks>1 data is loaded and processed to the
#' much lighter list of peaks in chunks batches.
#' @param prepf
#' Custom preprocessing function
#' @param spectra
#' Subset of spectra from \code{indir} to analyze. Without extension.
#' If \code{NULL} (default), all spectra is processed.
#' @param ncores
#' Number of cores used for the preprocessing of spectra. The cores will work in
#' parallel with the different spectra within a chunk.
#' @param iocores
#' Number of cores used for I/O operations. For some systems I/O operations involving
#' multiple cores reading or writing at the same time from disk increases time.
#' @param vch
#' Every how many chuncks progress is reported
#' @return A list of objects returned by prepf
#'
#' @importFrom parallel mcmapply detectCores mclapply
#' @export
#' @details
#' When using \code{"fread"} as the reading function, we can specify the separator \code{sep}.
#' \code{"table"} is only advised for cases where the separator is variable and
#' can be one or more spaces or tabs. See \code{\link[utils]{read.table}}.
#' \code{"fread"} uses internally \code{\link[data.table]{fread}}, which is faster.
#' \code{"mzml"} uses internally \code{\link[MALDIquantForeign]{importMzMl}}
#' @examples
preprocessData = function(indir, readf = c("fread", "table", "mzml"), sep=NULL,
                          nchunks = 50, prepf = NULL, spectra = NULL,
                          ncores = NULL, iocores = 1, vch = 5){
  if (is.null(ncores)){
    ncores = detectCores() - 2
  } else {
    ncores = min(detectCores() - 2, ncores)
  }

  cat("Using ", ncores, " cores\n")
  readf = match.arg(readf)
  switch(EXPR=readf,
         "fread" = {
           read_f = function(sep){
             function(f) importTsv(f, sep)
           }
           read_f = read_f(sep)
         },
         "table" = {
           read_f = importTable
         },
         "mzml" = {
           read_f = import_file.MzMl

         }
  )

  spectra_f = list.files(indir)
  # Remove extension
  spectra_f = strsplit(spectra_f, "\\.")
  ext = spectra_f[[1]][2]
  spectra_f = sapply(spectra_f, "[[", 1)
  if (!is.null(spectra)) {
    spectra_f = spectra_f[spectra_f %in% spectra]
  }

  if (nchunks > 1) {
    spectra_chunks = chunks(spectra_f, nchunks)
  } else {
    spectra_chunks = list(spectra_f)
  }
  names(spectra_chunks) = NULL

  peaks = mapply(
    function(x, ch, ext){
      if (ch %% vch == 0){
        cat(sprintf('Chunk %i of %i', ch, nchunks), "\n")
      }
      infiles = file.path(indir, paste0(x, '.', ext))
      l = mclapply(
        infiles,
        read_f,
        mc.cores=iocores,
        mc.silent=T
      )
      names(l) = x
      invisible(mcmapply(
        prepf, l,
        mc.cores=ncores, SIMPLIFY = F
      ))
    },
    spectra_chunks,
    seq_along(spectra_chunks),
    MoreArgs = list(ext=ext),
    SIMPLIFY = F
  )
  # Unlist chunks, so all spectra are in one list of depth=1
  peaks = unlist(peaks, recursive = F, use.names = T)
  return(peaks)
}



