


#' Build internal reference spectrum
#'
#' It comprises a peak binning followed by filtering step, adapted from [MALDIquant::referencePeaks].
#' This function allows to filter peaks groupwise.
#'
#' @param mpl MassPeaksList or Spectra object
#' @param minFreq Minimum frequency of a peak in the list of spectra to be added
#'                to the reference.
#' @param method Method for binning peaks across spectra. 'strict' does not allow
#'               a bin that contains two peaks from the same spectra and would
#'               split it further, while 'relaxed' would allow it.
#' @param tolerance Tolerance deviation value for grouping peaks into the same bin.
#' @param ppm Tolerance deviation in ppm. `tolerance` is ignored if `ppm` is set.
#' @param labels Labels for gorupwise peak filtering.
#' @param ret.object `character`, type of object to return:
#'        'Spectra', 'MassPeaks', or 'df'
#'
#' @return Reference peak list. The type depends on `ret.object`:
#'         [Spectra::Spectra], [MALDIquant::MassPeaks] or a data frame
#' @importFrom MALDIquant filterPeaks binPeaks createMassPeaks isMassPeaksList
#' @export
#'
int_ref_peaks = function(mpl, method=c('strict', 'relaxed'), minFreq=0.9,
                         tolerance=0.002, ppm=NULL, labels=NULL, ret.object=c('Spectra','MassPeaks', 'df')){
  method = match.arg(method)
  ret.object = match.arg(ret.object)

  if (!is.null(ppm)) {
    tolerance = ppm*1e-6
  }

  if (is(mpl, 'Spectra')) {
    # Apply processing queue
    peaks = peaksData(mpl)
    # Transform into MALDIquant MassPeaksList
    mpl = asMassPeaksList(peaks)
  } else if(!isMassPeaksList(mpl)) {
    stop('mpl must be either a Spectra object or a MassPeaksList')
  }

  if (is.null(labels)){
    ref_peaks = filterPeaks(
      binPeaks(mpl, method=method, tolerance=tolerance),
      minFrequency=minFreq
    )
  } else {
    ref_peaks = filterPeaks(
      binPeaks(mpl, method=method, tolerance=tolerance),
      minFrequency=minFreq, labels=labels
    )
  }
  m = get_bin_matrix(intensity_matrix(ref_peaks))
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



#' Parse protein into peptides using
#'
#' @param species Organism name
#' @param chain Protein chain name
#' @param id Protein ID, for example Uniprot ID
#' @param sequence Protein sequence
#'
#' @return Data frame with columns: "seq", "nhyd", "nglut", "mass1", "seqpos",
#'         "missed.cleaves", "chain", "species" and "id"
#'
#' @importFrom bacollite parse.seq
#' @export
#'
parse_seqs = function(species, chain, id, sequence, gpo_only=F){
  pseq = parse.seq(sequence, max.missed.cleaves = 0, gpo_only=gpo_only)
  pseq$chain = chain
  pseq$species = species
  pseq$id = id
  pseq
}


#' Build reference Spectrum from external sequences
#'
#' @param sequences Dataframe with collumns: species, chain, id, sequence
#' @param mc.cores Numbre of CPUs for parallel processing
#' @param gpo_only Logical. If `TRUE`, only the second proline in GPP patterns can be
#'        hydroxilated, becoming GPO. If `FALSE`, both can be hyroxylated.
#' @param ret.object Type of object to return. It can be:
#'  * 'df': for a data frame,
#'  * 'MassPeaks': for a [MALDIquant::MassPeaks] object
#'  * 'Spectra': for a [Spectra::Spectra] object
#' @param non_deam logical, whether to include deamidated peptides for those
#'        with N or Q aminoacids
#'
#' @return Reference peak list. The type depends on `ret.object`:
#'         [Spectra::Spectra], [MALDIquant::MassPeaks] or a data frame
#' @importFrom dplyr bind_rows select
#' @importFrom parallel mcmapply
#' @importFrom tibble as_tibble
#' @export
#'
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
#' This function aligns peaks across spectra using a specified tolerance using a reference peak list.
#' See [MALDIrppa::alignPeaks].
#'
#' @param s A `Spectra` object containing the spectral data to be aligned.
#' @param tolerance A numeric value specifying the m/z tolerance for peak matching.
#'        If `reference` is `NULL`, the same tolerance is used to build an internal reference.
#' @param ppm Tolerance deviation in ppm. `tolerance` is ignored if `ppm` is set.
#' @param minFreq A numeric value between 0 and 1 representing the minimum frequency of a peak across spectra for it to be considered a reference peak. Defaults to 0.9 if `reference` is `NULL`.
#' @param reference An optional [MALDIquant::MassPeaks] object to be used as the reference for alignment. If `NULL`, a reference is created based on `minFreq`.
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
align_peaks = function(s, tolerance=0.002, ppm=NULL, minFreq=NULL, reference=NULL, labels=NULL,
                       return_ref=FALSE, return_warp=FALSE,...){
  args_methods = list(...)

  if (!is.null(ppm)) {
    tolerance = ppm*1e-6
  }
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
#' Wrapper on [MALDIquant::binPeaks] to work on [Spectra::Spectra] objects.
#'
#' @param s Spectra object
#' @param tolerance Tolerance deviation value for grouping peaks into the same bin.
#' @param ppm Tolerance deviation in ppm. `tolerance` is ignored if `ppm` is set.
#' @param ... Paramteres for [MALDIquant::binPeaks], currently:
#'  - `method` Method for binning: \itemize{
#'  \item 'strict': A bin cannot contain peaks from the same spectra. In which case,
#'              it is further split.
#'  \item 'relaxed': This method allows multiple peaks from the same spectea in the same bin.
#'  \item 'reference': The spectra are binned around the peaks of the first spectra in 's'
#'}
#'
#' @return [Spectra::Spectra] object
#' @export
#' @importFrom MALDIquant mass intensity
#'
bin_peaks = function(s, method=c("strict", "relaxed", "reference"),
                     tolerance=0.002, ppm=NULL) {
  if (!is.null(ppm)) {
    tolerance = ppm*1e-6
  }

  # Apply processing queue
  peaks = peaksData(s)
  # Transform into MALDIquant MassPeaksList
  mpl = asMassPeaksList(peaks)

  mpl = binPeaks(mpl, tolerance, method=method)

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
#' @param tolerance tolerance between observed peaks and the expected pseudo cluster
#'        using `distance` between peaks.
#' @param ppm tolerance in parts per million. If use, `tolerance` is ignored.
#' @param ... Arguments passed to [MALDIquant::monoisotopicPeaks]
#'    - 'minCor': `double`, minimum correlation between the experimental pseudo isotopic pattern
#'       and the theoretical one generated by the Poisson model.
#'    - 'distance': mass distance between the peaks in the pattern. Default is 1.00235
#'    - 'size': sizes of isotopic patterns to look for. Using the default,
#'      `3L:10L`, it will look for patterns of size 10 prefferably and decrease it down
#'      to 3.
#'
#' @return [Spectra::Spectra] object with picked monoisotopic peaks
#' @export
#' @importFrom MALDIquant monoisotopicPeaks mass intensity
monoisotopic_peaks = function(s, tolerance, ppm, ...) {
  if (!is.null(ppm)){
    tolerance = ppm * 1e-6
  }
  # Apply processing queue
  peaks = peaksData(s)
  # Transform into MALDIquant MassPeaksList
  mpl = asMassPeaksList(peaks)
  mpl = monoisotopicPeaks(mpl, tolerance=tolerance, ...)

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


#' Filter peaks based on their frequency in the set of spectra
#'
#' It uses [MALDIquant::filterPeaks]
#' @param s [Spectra::Spectra] object
#' @param ... Parameters passed to [MALDIquant::filterPeaks]
#'    - minFrequency Minimum frequency (0 to 1) threshold for a peak to be kept
#'    - minNumber Absolute frquency, similar to `minFrequency`, but with
#'      absolute number of spectra
#'    - labels `factor` for groupwise filtering. `minFrequency` and `minNumber`
#'      are considered within each group.
#'    - mergeWhitelists `logical`, such that when `TRUE` and applying groupwise
#'      filtering with `labels`, peaks that pass the filter in one group, are also kept
#'      or whitelisted in other groups, even if their frequencies are below
#'      `minFrequency` or `minNumber`.
#'
#' @return [Spectra::Spectra] object after peak filtering
#' @export
#'
filter_peaks = function(s, ...) {
  # Apply processing queue
  peaks = peaksData(s)
  # Transform into MALDIquant MassPeaksList
  mpl = asMassPeaksList(peaks)

  mpl = filterPeaks(mpl, ...)

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
#' @param s [Spectra::Spectra] object
#'
#' @return A vector of integers with the number of peaks
#' @export
#'
count_peaks = function(s) {
  peaks = peaksData(s)
  return(sapply(peaks, nrow))
}

#' Intensity matrix
#'
#' Turns the list of peaks into an intensity matrix with samples/spectra as rows
#' and mz values as columns.
#'
#' @param s Spectra or MassPeaks object
#' @param specvar_to_names Spectra variable to use as rownames of the
#'        intensity matrix
#'
#' @return `matrix` with samples as rows and masses as columns, filled with
#'         intensity values
#' @export
#' @importFrom MALDIquant intensityMatrix
#'
intensity_matrix = function(s, specvar_to_names=NULL) {
  # Apply processing queue
  if (class(s) == 'Spectra') {
    peaks = peaksData(s)
    # Transform into MALDIquant MassPeaksList
    mpl = asMassPeaksList(peaks)
    m = intensityMatrix(mpl)
  } else if (isMassPeaksList(s)) {
    m = intensityMatrix(s)
  } else {
    stop('s needs to be a Spectra or a list of MassPeaks objects.')
  }
  # rownames(m) = s[[specvar_to_names]]
  attr(m, 'mass') = round(attr(m, 'mass'), 2)
  colnames(m) = attr(m, 'mass')
  if (!is.null(specvar_to_names)) {
    rownames(m) = s[[specvar_to_names]]
  }
  return(m)
}



#' Get binary matrix
#'
#' @param m intensity matrix
#'
#' @return Binary, 0-1 `matrix`
#' @export
#'
get_bin_matrix = function(m) {
  stopifnot(is.matrix(m))
  wn = is.na(m)
  m[wn] = 0L
  m[!wn] = 1L
  mode(m) = "integer"
  return(m)
}
