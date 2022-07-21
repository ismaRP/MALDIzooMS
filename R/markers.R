
#' Title
#'
#' @param markers
#' @param peaks
#' @param tolerance
#' @param match_tol
#'
#' @return
#' @importFrom MALDIquant createMassPeaks binPeaks monoisotopicPeaks
#' @importFrom MALDIquant intensityMatrix match.closest plot
#' @importFrom dplyr filter
#' @export
#'
#' @examples
pept_fly = function(markers, peaks, tolerance=0.002, match_tol = 0.6){

  ref_peaks = createMassPeaks(
    mass = markers$pymass,
    intensity = rep(1, nrow(markers)))
  # int_ref_peaks = intRefPeaks(peaks, 'strict', 0.8, 0.002)

  # Align peaks
  tolerance = 0.002
  al_peaks = custom_alignPeaks(
    peaks, reference=ref_peaks, tolerance = tolerance)
  # Final bin peaks, relaxed
  bin_peaks = binPeaks(al_peaks, method = 'strict', tolerance = tolerance)
  mip_peaks = monoisotopicPeaks(bin_peaks, minCor=0.9)

  # Get the feature matrix from the peaks using the spectra to interpolate
  intMatrix = intensityMatrix(mip_peaks)
  intMatrix[is.na(intMatrix)] = 0
  masses = round(as.numeric(colnames(intMatrix)), 4)
  colnames(intMatrix) = masses
  binMatrix = intMatrix
  binMatrix[binMatrix > 0] = 1

  pos = match.closest(ref_peaks@mass, masses, tolerance = match_tol)

  frac = colSums(binMatrix[,pos])/nrow(binMatrix)
  frac[is.na(frac)] = 0

  markers = markers %>%
    mutate(
      ms2conf = pept_ms2conf & nhyd_ms2conf & nglut_ms2conf,
      malditof_mass = masses[pos],
      malditof_frac = frac) %>%
    mutate(error_mass = pymass - malditof_mass)


  return(markers)

}



#' Title
#'
#' @param p
#' @param s
#'
#' @return
#' @importFrom bacollite ms_subrange ms_align ms_iso
#' @export
#'
#' @examples
align_pept = function(p, s){
  ts = ms_iso(p$seq, p$nglut, p$nhyd)
  moff = 1.5
  lbl <- min(ts$mass) - moff
  ubl <- max(ts$mass) + moff
  myxlim = c(lbl,ubl)

  subms = ms_subrange(s, lbl, ubl)

  al_cor =  ms_align(ts, subms, myxlim, gauss = 0.2, doplot=F)
  return(al_cor)
}


#' Title
#'
#' @param s
#' @param markers
#'
#' @return
#' @export
#'
#' @examples
align_sample = function(s, markers){
  s = as.matrix(s)
  markers_list = split(markers, seq(nrow(markers)))
  al = lapply(markers_list, align_pept, s)
  return(al)
}



#' Title
#'
#' @param markers
#' @param data_path
#' @param ncores
#' @param metadata
#'
#' @return
#' @export
#'
#' @examples
align_markers = function(markers, data_path, metadata, ncores, iocores){

  prepf = prepFun(align_sample, markers = markers_zooms)

  cor_data = preprocessData(
    indir = data_path, readf = 'mzml',
    nchunks = 2, prepf = prepf, ncores = ncores, iocores = iocores)

  cors = lapply(cor_data, function(x) unlist(lapply(x, function(x) x$cor)))
  cors = do.call(rbind, cors)
  rownames(cors) = sapply(strsplit(rownames(cors), '\\.'), '[[', 1)
  cors = cors[paste0(metadata$sample_name, '_', metadata$replicate), ]


  lags = lapply(cor_data, function(x) unlist(lapply(x, function(x) x$lag)))
  lags = do.call(rbind, lags)
  rownames(lags) = sapply(strsplit(rownames(lags), '\\.'), '[[', 1)
  lags = lags[paste0(metadata$sample_name, '_', metadata$replicate), ]

  return(list(cors, lags))

}



