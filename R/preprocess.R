


#' Quality control spectra from [MALDIrppa::screenSpectra()]
#' Works on peak matrix rather than on [MALDIquant::MassSpectrum]
#' @param x Raw peak matrix with mz and intensities
#'
#' @return numeric, Atipicality score
#' @export
#' @importFrom signal sgolayfilt
#' @importFrom robustbase Qn
#' @importFrom stats median
atipicality_spectra = function(x, ...){
  smax = 100
  nd = 1
  lambda = 0.5
  est = x[,2]/max(x[,2])/smax
  est = sgolayfilt(est, m = nd)
  est = Qn(est)
  med.int = round(median(x[,2]), 4)
  est = (est^lambda)*(1/sqrt((med.int+1)))^(1-lambda)

  return(est)
}


#' Identification of potentially low-quality raw mass spectra
#'
#' Uses [atipicality_spectra], which implements [MALDIrppa::screenSpectra()] on
#' peak matrix
#' @param object [Spectra::Spectra()] object
#'
#' @return A [Spectra::Spectra()] object with \code{QCflag} spectrum variable added
#' @export
#' @importFrom Spectra addProcessing
screen_spectra = function(object){
  qcA = addProcessing(sps_mzr, atipicality_spectra)
  qcA = unlist(peaksData(qcA))

  sps_mzr$Atipicality = qcA
  l = qc_limits(sps_mzr$Atipicality)

  sps_mzr$QCflag = 'QC_pass'
  sps_mzr$QCflag[sps_mzr$Atipicality > l[1] | sps_mzr$Atipicality < l[2]] = 'QC_fail'

  sps_mzr = reset(sps_mzr)

  return(sps_mzr)
}





#' Calculate outlier limits for Atipicality score using [robustbase::adjboxStats]
#'
#' @param A vector of A scores from QC
#'
#' @return Vector with top and bottom limits
#' @export
#' @importFrom robustbase adjboxStats
#'
qc_limits = function(A) {
  threshold = 1.5
  QClimits = adjboxStats(A, coef = threshold)$fence[c(2,1)]
  return(QClimits)
}


#' Place desired intensity column in 2nd position
#'
#' Many methods automatically pick the 2nd column as the intensity to
#' do operations. This function we can swap the intensity columns to have
#' the desired one in the 2nd column
#'
#' @param x Matrix of peaks
#' @param int_index Column index of the intensity
#' @param ... Not in use, required by [Spectra::addProcessing]
#'
#' @return Peak matrix
#' @export
#'
int_col2 = function(x, int_index=4, ...) {
  # Swap intensity
  tmp = x[,2]
  x[, 2] = x[,int_index]
  x[, int_index] = tmp
  # Swap colnames
  tmp_colname = colnames(x)[2]
  colnames(x)[2] = colnames(x)[int_index]
  colnames(x)[int_index] = tmp_colname

  return(x)
}



#' Applies processing queue to Spectra
#'
#' @param s Spectra object
#' @param write_data logical. Whether to write spectra data into a file
#'        This will apply processing, export spectra into mzML file(s) and
#'        reload into a new Spectra object with MsBackendMzR. This is useful
#'        if we want to apply preprocessing and the resulting object is still
#'        too big to keep in memory, either because there are too many spectra
#'        or because we still haven't picked peaks.
#' @param file If `write_data` is `TRUE`, file or files to write spectra.
#'        IF there is one file, all spectra is saved on the same file.
#'        If multiple files, there needs to be one per spectra.
#'
#' @return A Spectra object with either [Spectra::MsBackendMzR] or
#' [Spectra::MsBackendDataFrame] backends
#' @importFrom Spectra export Spectra MsBackendMzR MsBackendDataFrame
#' @export
#'
#' @examples
apply_preprocess = function(s, write_data=FALSE, file=NULL) {

  if (write_data) {
    # Exporting applies all processing steps
    export(object=sps_mzr, backend=MsBackendMzR(), format='mzML',
           file=file)
    # Reload the mzML into a new Spectra object
    s = suppressMessages(
      Spectra(file, source = MsBackendMzR(), centroided = FALSE,
              BPPARAM = MulticoreParam(workers = ncores)))
  } else {
    peaks = peaksData(s)
    spd = spectraData(s)

    spd$mz = lapply(peaks, '[', , 1)
    spd$intensity = lapply(peaks, '[', , 2)
    if ('SNR' %in% colnames(peaks[[1]])) {
      spd$SNR = lapply(peaks, "[", ,3)
    }
    s = Spectra(
      spd, backend = MsBackendDataFrame(), centroided=TRUE,
      peaksVariables = peaks_vars)
  }
  return(s)

}


#' Peak detection
#'
#' It implements [Spectra::pickPeaks()], but allows using intensities in a given
#' column in the peak matrix
#'
#' @param x Peak matrix
#' @param halfWindowSize Half size of the sliding window used to detect local
#'                       maxima.
#' @param method Method for noise estimation. Either \code{"SuperSmoother"} (default)
#'               or \code{"MAD"}
#' @param snr Signal-to-noise ratio threshold above which a peak is picked
#' @param k Parameter for [MsCoreUtils::refineCentroids]. Values to the left and right
#'          of a peak to be considered to refine the peak mz position
#' @param descending Parameter for [MsCoreUtils::refineCentroids]. Only values between
#'                   nearest valleys are used.
#' @param threshold Parameter for [MsCoreUtils::refineCentroids], a proportion.
#'                  Only values above this fraction of the maxima are used.
#' @param int_index Index of the intensity in the peak matrix to be used for peak
#'                  detection.
#' @param add_snr logical, whether to add a column with each peak's SNR
#' @param ... Currently not in use, required by [Spectra::addProcessing].
#'
#' @return Peak matrix, only containing peaks that fullfill all criteria.
#' @importFrom MsCoreUtils localMaxima noise refineCentroids
#' @export
#'
peak_detection = function(x, halfWindowSize = 2L, method = c("MAD", "SuperSmoother"),
                         snr = 0L, k = 0L, descending = FALSE, threshold = 0,
                         int_index=2, add_snr=FALSE, ...){
  n = noise(x[, 1L], x[, int_index], method = method)
  l = localMaxima(x[, int_index], hws = halfWindowSize)
  p = which(l & x[, int_index] > (snr * n))
  n = n[p]
  if (k > 0L) {
    mz = refineCentroids(x = x[, 1L], y = x[, int_index], p = p,
                         k = k, threshold = threshold,
                         descending = descending)
    x = x[p, , drop = FALSE]
    x[,1] = mz
    # cbind(mz = refineCentroids(x = x[, 1L], y = x[, 2L], p = p,
    #                            k = k, threshold = threshold,
    #                            descending = descending),
    #       intensity = x[p, 2L])
  } else {
    x = x[p, , drop = FALSE]
  }
  # Add SNR
  if (add_snr){
    x = cbind(x, x[, int_index]/n)
    colnames(x)[ncol(x)] = 'SNR'
  }
  return(x)
}




#' Baseline correction subtraction
#'
#' @param x Peak matrix
#' @param int_index Index of the intensity to calculate the baseline
#' @param keep_bl Keep basseline in a separate column
#' @param substract_index Index of the intensity column to be substracted.
#' If NULL, there is no substraction
#' @param in_place Replace the intensity in \code{substract_index} with the baseline
#' substracted intensity
#' @param ... Arguments passed to [MsCoreUtils::estimateBaseline()]. Currently \code{'method'}.
#'
#' @return Peak matrix. Depending on the arguments, with the baseline-substracted
#' intensities either in the same or in a new column and with or without the
#' baseline itself.
#' @export
#' @importFrom MsCoreUtils estimateBaseline
baseline_correction = function(x, int_index=2, keep_bl=TRUE,
               substract_index=2, in_place=FALSE, ...) {
  method = list(...)$method
  # iterations = list(...)$iterations
  # decreasing = list(...)$decreasing
  # b = MsCoreUtils::estimateBaseline(
  #   x[,1], x[, int_index], method=method,
  #   iterations=iterations, decreasing=decreasing)
  b = MsCoreUtils::estimateBaseline(
    x[,1], x[, int_index], method=method)
  # Substract baseline
  if (!is.null(substract_index)){
    subs_int = x[, substract_index, drop=FALSE] - b

    if (is.character(substract_index)){
      colnames(subs_int) = paste0(
        substract_index, '_bl_corr_', method)
    } else if (is.integer(substract_index)) {
      colnames(subs_int) = paste0(
        colnames(x)[substract_index], '_bl_corr_', method)
    }

    if (in_place) {
      x[, substract_index] = subs_int
    } else {
      x = cbind(x, subs_int)
    }
  }
  # Keep baseline or not
  if (keep_bl){
    x = cbind(x, b)
    colnames(x)[ncol(x)] = paste0('baseline_', method)
  }
  return(x)
}


#' Intensity smoothing
#'
#' @param x Peak matrix
#' @param method Smoothing method. One of \code{'MovingAverage'}, \code{'WeightedMovingAverage'}
#' or \code{'SavitzkyGolay'}
#' @param hws Half window size
#' @param k For Savitzky-Golay filter, this is the polynomial order for the coefficients
#' @param int_index Column with the intensity to smooth
#' @param in_place Whether smoothing is in place, or a new column with the smoothened
#'                 intensity is added
#' @param ... Parameters passed to other methods. Not in use, required by [Spectra::addProcessing].
#'
#' @return Peak matrix
#' @export
#' @importFrom MsCoreUtils smooth coefMA coefWMA coefSG
smooth = function(x,
                  method=c('MovingAverage', 'WeightedMovingAverage', 'SavitzkyGolay'),
                  hws=4L, k=3L, int_index=2, in_place=FALSE, ...){
  method <- match.arg(method)
  switch(
    method,
    MovingAverage={
     coefs = coefMA(hws)
    },
    WeightedMovingAverage={
     coefs = coefWMA(hws)
    },
    SavitzkyGolay={
     coefs = coefSG(hws, k)
    }
  )
  int_vals = x[, int_index]
  smoothed_int = MsCoreUtils::smooth(int_vals, coefs)
  if (in_place){
    x[, int_index] = smoothed_int
  } else {
    x = cbind(x, smoothed_int)
    colnames(x)[ncol(x)] = paste0('intensity_', method)
  }
  return(x)

}



#' peptidePseudoClusters
#'
#' Extracts pseudo envelopes from given monoisotopic peptide masses
#' It just checks whether the monoisotopic peaks and subsequent peaks
#' at distance of 1.00235 are present. There muse be between \code{min_isopeaks}
#' and \code{n_isopeaks}.
#' It also check there are no gaps. It doesn't check the isotopic envelop shape.
#'
#' @param x Peaks matrix object
#' @param mono_masses Monoisotopic masses
#' @param n_isopeaks Number of isotopic peaks to look for
#' @param min_isopeaks Minimum number of isotopic peaks to consider the envelope
#' @param ... Argumentes passed to [MsCoreUtils::closest]
#'
#' @return Peaks matrix with isotopic masses and NA as place holders for when
#' there is no match
#' @export
#' @importFrom MsCoreUtils closest
#'
#' @examples
#' masses = c(100, 200, 300)
#' n_isopeaks = 5
#' masses = matrix(masses, nrow = n_isopeaks, ncol = length(masses), byrow = TRUE)
#' d = 1.00235
#' masses = masses + (d * 0L:(n_isopeaks - 1L))
#' masses = sort(masses)
#' print(masses)
peptide_pseudo_clusters = function(x, mono_masses, n_isopeaks, min_isopeaks, ...){

  # Create isotopic masses from monoisotopic peptide masses
  npepts = length(mono_masses)
  masses = matrix(mono_masses, nrow = n_isopeaks, ncol = length(mono_masses), byrow = T)
  d = 1.00235
  masses = masses + (d * 0L:(n_isopeaks - 1L))
  masses = sort(masses)
  # Empty matrix with selected masses and intensities
  sel_matrix = matrix(nrow=length(masses), ncol=ncol(x))
  # Run closest function to match x and isotopic masses
  # idx = closest(x[,1], masses, tolerance = x[,1]*tol)
  # idx = closest(x[,1], masses, ppm=250, tolerance=0)
  idx = MsCoreUtils::closest(x[,1], masses, ...)

  # Vector to control for isotopic pattern completness
  sel_complete = rep(F, length(masses))
  sel_complete[idx[!is.na(idx)]] = T
  sel_complete = split(sel_complete,
                       cut(seq_along(sel_complete), npepts, labels = FALSE))
  # If it has a gap or is shorter than min_isopeaks we will just ignore it
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
  # Put selected masses in place in sel_matrix and keep only complete patterns
  sel_matrix[idx[!is.na(idx)], ] = x[!is.na(idx),]
  sel_matrix[!sel_complete, ] = NA
  delta_mass = sel_matrix[,1] - masses
  sel_matrix[,1] = masses
  sel_matrix = cbind(sel_matrix, delta_mass)

  colnames(sel_matrix) = c(colnames(x), 'delta_mass')
  return(sel_matrix)
}




#' Model local backgroung noise around peak
#'
#' @param p
#' Peak as one-row matrix, mass-intensity pair
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
#' @return One-row matrix with p-value, ie., probability of a peak with a higher intensity,
#' normal distribution fitting log-likelihood and p-value.
#' @importFrom MASS fitdistr
#' @importFrom stats dnorm p.adjust pchisq pnorm quantile
#' @importFrom graphics hist
#' @export
#'
model_local_bg = function(p, full_m, mass_range, bg_cutoff){

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


#' Filter a peak list based on likelihood that peaks are above background local
#' noise
#'
#' @param x Peaks matrix
#' @param mass_range
#' Mass window to both sides of a peak to be considered for backgroun modelling
#' @param bg_cutoff
#' The peaks within the mass range with intensity below the \code{bg_cutoff} quantile
#' are considered for background modelling. \code{bg_cutoff=1} keeps all peaks
#' and \code{bg_cutoff=0.5} would only keep the bottom half.
#' @param l_cutoff
#' Likelihood threshold or p-value. Peaks with a probability of being modelled as
#' background noise higher than this are filtered out.
#' @return Peak matrix
#' @export
#'
peaks_local_bg = function(x, mass_range, bg_cutoff, l_cutoff,
                          int_index=2, ...){
  s = cbind(x[,1, drop=F], x[,int_index, drop=F])
  total_peaks = nrow(s)
  # m = as.matrix(s)

  l_values = apply(
    X=s, MARGIN=1,
    FUN=model_local_bg,
    full_m=s, mass_range=mass_range, bg_cutoff=bg_cutoff)
  l_values = t(l_values)

  adj_pval = p.adjust(l_values[,1], 'holm')

  l_values = cbind(l_values, adj_pval)

  colnames(l_values) = c("bg_lik", "fit_loglik", "chsq_pval", "bg_adjlik")

  peaks_mask = l_values[,1] < l_cutoff
  x = x[peaks_mask, ]

  non_bg_peaks = nrow(x)
  frac_rem = non_bg_peaks/total_peaks

  return(x)
}


