


#' Title
#'
#' @param m
#'
#' @return
#'
#' @examples
as.binary.matrix = function(m){
  stopifnot(is.matrix(m))
  isNA = which(is.na(m))
  m[] = 1L
  m[isNA] = 0L
  mode(m) = "integer"
  m
}


#' Title
#'
#' @param l
#'
#' @return
#'
#' @examples
as.matrix.MassObjectList = function(l){
  mass = unlist(
    lapply(l, function(x)x@mass),
    recursive=FALSE, use.names=FALSE)
  intensity = unlist(lapply(l, function(x)x@intensity),
                     recursive=FALSE, use.names=FALSE)
  uniqueMass = sort.int(unique(mass))
  n = lengths(l)
  r = rep.int(seq_along(l), n)

  i = findInterval(mass, uniqueMass)

  m = matrix(NA_real_, nrow=length(l), ncol=length(uniqueMass),
             dimnames=list(NULL, uniqueMass))
  m[cbind(r, i)] = intensity
  attr(m, "mass") = uniqueMass
  m
}


#' Build internal reference spectrum
#'
#' @param mpl MassPeaksList object
#' @param minFreq
#' @param method
#' @param tolerance
#' @param labels
#'
#' @return
#' @importFrom MALDIquant filterPeaks binPeaks createMassPeaks isMassPeaksList
#' @export
#'
#' @examples
intRefPeaks = function(mpl, method=c('strict', 'relaxed'), minFreq=0.9,
                       tolerance=0.002, labels=NULL){
  method = match.arg(method)

  if (is(mpl, 'Spectra')) {
    # Apply processing queue
    peaks = peaksData(s)
    # Transform into MALDIquant MassPeaksList
    mpl = asMassPeaksList(peaks)
  } else if(!isMassPeaksList(mpl)) {
    stop('mpl must be either a Spectra object or a MassPeaksList')
  }

  if (is.null(labels)){
    referencePeaks = filterPeaks(
      binPeaks(mpl, method=method, tolerance=tolerance),
      minFrequency=minFreq
    )
  } else {
    referencePeaks = filterPeaks(
      binPeaks(mpl, method=method, tolerance=tolerance),
      minFrequency=minFreq, labels=labels
    )
  }


  m = as.binary.matrix(as.matrix.MassObjectList(referencePeaks))

  ## set peak intensity to number of occurrence
  intensity = unname(colMeans(m))

  createMassPeaks(mass=attr(m, "mass"), intensity=intensity)
}



#' Parse sequence
#'
#' @param species
#' @param chain
#' @param id
#' @param sequence
#'
#' @return
#'
#' @importFrom bacollite parse.seq
#' @export
#'
#' @examples
parse_seqs = function(species, chain, id, sequence, gpo_only=F){
  pseq = parse.seq(sequence, max.missed.cleaves = 1, gpo_only=gpo_only)
  pseq$chain = chain
  pseq$species = species
  pseq
}


#' Build reference Spectrum from external sequences
#'
#' @param sequences
#' @param mc.cores
#' @param gpo_only
#' @param ret.object
#'
#' @return
#' @importFrom dplyr bind_rows select
#' @importFrom parallel mcmapply
#' @importFrom tibble as_tibble
#' @export
#'
#' @examples
extRefPeaks = function(sequences, mc.cores=4L, gpo_only=F, ret.object=c('df', 'MassPeaks', 'Spectra')) {

  ret.object = match.arg(ret.object)
  # Digest sequences
  peptides = do.call(
    mcmapply,
    c(list('FUN'=parse_seqs, 'mc.cores'=mc.cores, 'SIMPLIFY'=F,
           'USE.NAMES'=F, 'MoreArgs'=list('gpo_only'=gpo_only)),
      as.list(sequences)))

  peptides = do.call(bind_rows, peptides) %>% as_tibble()

  # peptides = peptides[peptides$nglut==0,]
  peptides = peptides[!duplicated(select(peptides, -species, -seqpos)),] %>%
    select(-species)

  if (ret.object == 'MassPeaks') {
    return(createMassPeaks(mass=sort(peptides$mass1), intensity=rep(1, nrow(peptides))))
  } else if (ret.object == 'Spectra') {
    spd = DataFrame(msLevel=1L, id='Reference')
    spd$intensity = list(rep(1, nrow(peptides)))
    spd$mz = list(sort(peptides$mass1))
    sps = Spectra(spd, backend=MsBackendDataFrame())
    return(sps)
  }

}



#' Align peaks
#'
#' Prior aligning, all processing queue is applied
#' @param s Spectra object. Must be centroided
#' @param minFreq
#' @param tolerance
#' @param labels
#' @param th_peaks
#' @param ...
#'
#' @return A Spectra object with MsBackendDataFrame
#' @importFrom MALDIquant determineWarpingFunctions warpMassPeaks binPeaks
#' @importFrom MALDIquant mass intensity
#' @importFrom Spectra MsBackendDataFrame
#' @export
#'
#' @examples
align_peaks = function(s, tolerance, minFreq=NULL, reference=NULL, labels=NULL, ...){
  args_methods = list(...)

  if (is.null(args_methods$allowNoMatches)) args_methods$allowNoMatches = F
  if (is.null(args_methods$emptyNoMatches)) args_methods$emptyNoMatches = F
  # Apply processing queue
  peaks = peaksData(s)
  # Transform into MALDIquant MassPeaksList
  mpl = asMassPeaksList(peaks)

  if (is.null(reference)){
    if (is.null(minFreq)) minFreq=0.9
    reference = intRefPeaks(
      mpl,
      minFreq = minFreq,
      method = "strict",
      tolerance = tolerance,
      labels = labels
    )
  }

  warpingFunctions = determineWarpingFunctions(
    mpl,
    reference = reference,
    tolerance = tolerance,
    method='lowess',
    allowNoMatches = args_methods$allowNoMatches
  )
  mpl = warpMassPeaks(mpl, warpingFunctions, emptyNoMatches = args_methods$emptyNoMatches)
  # l = binPeaks(l, method = "strict", tolerance = tolerance)
  # l = binPeaks(l, method = "relaxed", tolerance = tolerance)

  # Create new Spectra with MsBackendDataFrame and aligned spectra
  spd = spectraData(sps_mzr)
  spd$mz = lapply(mpl, MALDIquant::mass)
  spd$intensity = lapply(mpl, MALDIquant::intensity)
  peaks_vars = c('mz', 'intensity')
  if ('SNR' %in% colnames(peaks[[1]])) {
    spd$SNR = lapply(peaks, "[", ,3)
    peaks_vars = c(peaks_vars, 'SNR')
  }
  s_aligned = Spectra(
    spd, backend = MsBackendDataFrame(), centroided=TRUE,
    peaksVariables = peaks_vars)
  return(s_aligned)
}



#' Bin peaks across spectra
#'
#' @param s
#' @param ethod
#' @param tolerance
#'
#' @return
#' @export
#' @importFrom MALDIquant mass intensity
#'
#' @examples
bin_peaks = function(s, method=c("strict", "relaxed", "reference"), tolerance=0.002) {
  # Apply processing queue
  peaks = peaksData(s)
  # Transform into MALDIquant MassPeaksList
  mpl = asMassPeaksList(peaks)

  mpl = binPeaks(mpl, method, tolerance)

  spd = spectraData(sps_mzr)
  spd$mz = lapply(mpl, MALDIquant::mass)
  spd$intensity = lapply(mpl, MALDIquant::intensity)
  peaks_vars = c('mz', 'intensity')
  if ('SNR' %in% colnames(peaks[[1]])) {
    spd$SNR = lapply(peaks, "[", ,3)
    peaks_vars = c(peaks_vars, 'SNR')
  }
  s_binned = Spectra(
    spd, backend = MsBackendDataFrame(), centroided=TRUE,
    peaksVariables = peaks_vars)
  return(s_binned)
}



#' Filter peaks
#'
#' @param s
#' @param minFrequency
#' @param minNumber
#' @param labels
#' @param mergeWhitelists
#'
#' @return
#' @export
#'
#' @examples
filter_peaks = function(s, minFrequency, minNumber, labels, mergeWhitelists=FALSE) {
  # Apply processing queue
  peaks = peaksData(s)
  # Transform into MALDIquant MassPeaksList
  mpl = asMassPeaksList(peaks)

  mpl = filterPeaks(mpl,  minFrequency, minNumber, labels, mergeWhitelists=FALSE)

  spd = spectraData(sps_mzr)
  spd$mz = lapply(mpl, MALDIquant::mass)
  spd$intensity = lapply(mpl, MALDIquant::intensity)
  peaks_vars = c('mz', 'intensity')
  if ('SNR' %in% colnames(peaks[[1]])) {
    spd$SNR = lapply(peaks, "[", ,3)
    peaks_vars = c(peaks_vars, 'SNR')
  }
  s_filtered = Spectra(
    spd, backend = MsBackendDataFrame(), centroided=TRUE,
    peaksVariables = peaks_vars)
  return(s_filtered)
}


#' Intensity matrix
#'
#' Turns the list of peaks into an intensity matrix with samples/spectra as rows
#' and mz values as columns.
#'
#' @param s Spectra object
#'
#' @return
#' @export
#' @importFrom MALDIquant intensityMatrix
#'
#' @examples
intensity_matrix = function(s, specvar_to_names) {
  # Apply processing queue
  peaks = peaksData(s)
  names(peaks) = s[[specvar_to_names]]
  # Transform into MALDIquant MassPeaksList
  mpl = asMassPeaksList(peaks)
  m = intensityMatrix(mpl)
  rownames(m) = s[[specvar_to_names]]
  return(m)
}
