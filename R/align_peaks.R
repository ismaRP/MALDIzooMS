


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
int_ref_peaks = function(mpl, method=c('strict', 'relaxed'), minFreq=0.9,
                       tolerance=0.002, labels=NULL, ret.object=c('Spectra','MassPeaks', 'df')){
  method = match.arg(method)
  ret.object = match.arg(ret.object)

  if (is(mpl, 'Spectra')) {
    # Apply processing queue
    peaks = peaksData(mpl)
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

  if (ret.object == 'MassPeaks') {
    return(createMassPeaks(mass=attr(m, "mass"), intensity=intensity))
  } else if (ret.object == 'df') {
    return(data.frame(mass=attr(m, "mass"), intensity=intensity))
  } else if (ret.object == 'Spectra') {
    spd <- DataFrame(
      msLevel = c(1L),
      polarity = c(1L),
      id = c("IntReference"),
      name = c("IntReference"))
    spd$mz = list(attr(m, "mass"))
    spd$intensity = list(intensity)
    return(Spectra(spd))
  } else {
    stop('Return object can only be a MassPeaks or a dataframe (df)')
  }
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
  pseq = parse.seq(sequence, max.missed.cleaves = 0, gpo_only=gpo_only)
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
#' @param non_deam
#'
#' @return
#' @importFrom dplyr bind_rows select
#' @importFrom parallel mcmapply
#' @importFrom tibble as_tibble
#' @export
#'
#' @examples
ext_ref_peaks = function(sequences, mc.cores=4L, gpo_only=F, ret.object=c('df', 'MassPeaks', 'Spectra'),
                         non_deam = FALSE) {

  ret.object = match.arg(ret.object)
  # Digest sequences
  peptides = do.call(
    mcmapply,
    c(list('FUN'=parse_seqs, 'mc.cores'=mc.cores, 'SIMPLIFY'=F,
           'USE.NAMES'=F, 'MoreArgs'=list('gpo_only'=gpo_only)),
      as.list(sequences)))

  peptides = do.call(bind_rows, peptides) %>% as_tibble()
  # peptides = peptides[peptides$nglut==0,]

  peptides = peptides %>% arrange(mass1)
  peptides = peptides[!duplicated(select(peptides, mass1)),] %>%
    select(-species)

  if (non_deam) {
    peptides = peptides %>% filter(nglut == 0)
  }

  if (ret.object == 'MassPeaks') {
    return(createMassPeaks(mass=sort(peptides$mass1), intensity=rep(1, nrow(peptides))))
  } else if (ret.object == 'Spectra') {
    spd = DataFrame(msLevel=1L, id='Reference')
    spd$intensity = list(rep(1, nrow(peptides)))
    spd$mz = list(sort(peptides$mass1))
    sps = Spectra(spd, backend=MsBackendDataFrame())
    return(sps)
  } else if (ret.object == 'df'){
    return(peptides)
  }

}


#' Align Peaks in Spectra
#'
#' This function aligns peaks across spectra using a specified tolerance and optionally a reference peak list.
#'
#' @param s A `Spectra` object containing the spectral data to be aligned.
#' @param tolerance A numeric value specifying the m/z tolerance for peak matching.
#' @param minFreq A numeric value between 0 and 1 representing the minimum frequency of a peak across spectra for it to be considered a reference peak. Defaults to 0.9 if `reference` is `NULL`.
#' @param reference An optional `MassPeaks` object to be used as the reference for alignment. If `NULL`, a reference is created based on `minFreq`.
#' @param labels Optional labels for the reference peaks.
#' @param return_ref Logical. If `TRUE`, the function returns a list with the aligned spectra and the reference peaks. Defaults to `FALSE`.
#' @param return_warp Whether warping functions are returned, added as a Spectra varaible
#' @param ... Additional arguments to control alignment behavior. Supported arguments:
#'   - `allowNoMatches`: Logical, if `TRUE`, allows no matches in warping functions.
#'   - `emptyNoMatches`: Logical, if `TRUE`, inserts empty values for unmatched peaks in the output.
#'
#' @return If `return_ref = FALSE`, returns an aligned `Spectra` object. If `return_ref = TRUE`, returns a list containing the aligned `Spectra` object and the reference peaks.
#'
#' @details The function first processes the peaks in `s`, and if no reference is provided, it generates one based on `minFreq`. The function applies warping functions determined by the `determineWarpingFunctions` function and then aligns the spectra.
#'
#' @examples
#' # Example usage:
#' # align_peaks(spectra, tolerance = 0.1)
#' @importFrom MALDIquant determineWarpingFunctions warpMassPeaks binPeaks
#' @importFrom MALDIquant mass intensity
#' @importFrom Spectra MsBackendDataFrame peaksData
#' @export
#'
align_peaks = function(s, tolerance, minFreq=NULL, reference=NULL, labels=NULL,
                       return_ref=FALSE, return_warp=FALSE,...){
  args_methods = list(...)

  if (is.null(args_methods$allowNoMatches)) args_methods$allowNoMatches = F
  if (is.null(args_methods$emptyNoMatches)) args_methods$emptyNoMatches = F
  # Apply processing queue
  peaks = peaksData(s)
  # Transform into MALDIquant MassPeaksList
  mpl = asMassPeaksList(peaks)

  if (is.null(reference)){
    if (is.null(minFreq)) minFreq=0.9
    reference = int_ref_peaks(
      mpl,
      minFreq = minFreq,
      method = "strict",
      tolerance = tolerance,
      labels = labels,
      ret.object = 'MassPeaks'
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
  spd = spectraData(s)
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
  if (return_warp) s_aligned$warpFunc = warpingFunctions

  if (return_ref)  return(list(s_aligned, reference))
  else return(s_aligned)
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

  spd = spectraData(s)
  spd$mz = lapply(mpl, MALDIquant::mass)
  spd$intensity = lapply(mpl, MALDIquant::intensity)
  spd$peaksCount = sapply(mpl, function(x)length(x@mass))
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


#' Find monoisotopic peaks
#' It is a wrapper on [MALDIquant::monoisotopicPeaks] to work with Spectra objects
#' @param s Spectra object
#' @param ppm tolerance in parts per million
#' @param ...
#' \itemize{
#'    \item minCor
#'    \item distance
#'    \item size
#' }
#' @return
#' @export
#' @importFrom MALDIquant monoisotopicPeaks mass intensity
#' @examples
monoisotopic_peaks = function(s, ppm, ...) {
  # Apply processing queue
  peaks = peaksData(s)
  # Transform into MALDIquant MassPeaksList
  mpl = asMassPeaksList(peaks)
  mpl = monoisotopicPeaks(mpl, tolerance=ppm*1e-6, ...)

  spd = spectraData(s)
  spd$mz = lapply(mpl, MALDIquant::mass)
  spd$intensity = lapply(mpl, MALDIquant::intensity)
  spd$peaksCount = sapply(mpl, function(x)length(x@mass))
  peaks_vars = c('mz', 'intensity')

  if ('SNR' %in% colnames(peaks[[1]])) {
    spd$SNR = lapply(peaks, "[", ,3)
    peaks_vars = c(peaks_vars, 'SNR')
  }

  s_mono = Spectra(
    spd, backend = MsBackendDataFrame(), centroided=TRUE,
    peaksVariables = peaks_vars)
  return(s_mono)
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

  mpl = filterPeaks(mpl,  minFrequency, minNumber, labels, mergeWhitelists=mergeWhitelists)

  spd = spectraData(s)
  spd$mz = lapply(mpl, MALDIquant::mass)
  spd$intensity = lapply(mpl, MALDIquant::intensity)
  spd$peaksCount = sapply(mpl, function(x)length(x@mass))
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


#' Count peaks in each spectra
#'
#' @param s
#'
#' @return A vector of integers with the number of peaks
#' @export
#'
#' @examples
count_peaks = function(s) {
  peaks = peaksData(s)
  return(sapply(peaks, nrow))
}

#' Intensity matrix
#'
#' Turns the list of peaks into an intensity matrix with samples/spectra as rows
#' and mz values as columns.
#'
#' @param s Spectra object
#' @param specvar_to_names
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
  # rownames(m) = s[[specvar_to_names]]
  attr(m, 'mass') = round(attr(m, 'mass'), 2)
  colnames(m) = attr(m, 'mass')
  return(m)
}



#' Get binary matrix
#'
#' @param m intensity matrix
#'
#' @return
#' @export
#'
#' @examples
get_bin_matrix = function(m) {
  stopifnot(is.matrix(m))
  wn = is.na(m)
  m[wn] = 0L
  m[!wn] = 1L
  mode(m) = "integer"
  return(m)
}
