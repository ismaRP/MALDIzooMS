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
model_local_bg = function(p, bg){
  fit = fitdistr(bg, 'normal')

  n = length(bg)
  sd0 = sqrt((n-1)/n)*sd(bg)
  mx = mean(bg)
  msd = sd0/sqrt(n)
  # loglik = sum(dnorm(bg, mx, sd0, log=TRUE))
  # Chi-square test
  # h = hist(bg, plot=F, breaks = 20)
  # binwidth = h$mids[2] - h$mids[1]
  # co = h$counts
  # ce = dnorm(h$mids, mean= fit$estimate[1], sd = fit$estimate[2]) * length(bg) * binwidth
  # chsq = sum((co-ce)^2/ce)
  # pval = 1 - pchisq(chsq, length(bg)-1)

  l = pnorm(p, mean = mx, sd = sd0, lower.tail = F)

  return(c(l, NA, msd))

}



# peaks_local_bg = function(x, mass_range, bg_cutoff, l_cutoff
#                           int_index=2L, ...){
#   # x = as.matrix(x)
#   min_bg_peaks = 15
#   l_values = matrix(NA, nrow=nrow(x), ncol=4)
#
#   for (i in 1:nrow(x)) {
#     p = x[i,]
#     bg_mask = abs(p[[1]] - x[,1]) <= mass_range
#     # bg_mask = (x[,1] < (p[1] + mass_range)) & (x[,1] > (p[1] - mass_range))
#     bg = x[bg_mask, int_index]
#     int_cutoff = quantile(bg, bg_cutoff)
#     bg = bg[bg < int_cutoff]
#     n = length(bg)
#     if (n < min_bg_peaks) {
#       ldata = c(0, NA, NA)
#     } else {
#       print(p)
#       print(bg)
#       print(int_index)
#       ldata = model_local_bg(p[[int_index]], bg)
#     }
#     l_values[i, ] = c(ldata, n)
#   }
#
#   # adj_pval = p.adjust(l_values[,1], 'holm')
#   # l_values = cbind(l_values, adj_pval)
#
#   colnames(l_values) = c("bg_lik", "fit_loglik", "fit_error", 'n_bg')
#
#   peaks_mask = l_values[,1] < l_cutoff
#   x = cbind(x[peaks_mask, ], l_values[peaks_mask,])
#
#   return(x)
# }



#' Filter a peak list based on likelihood that peaks are above background local
#' noise
#'
#' @param x Peaks matrix
#' @param mass_range
#' Mass window to both sides of a peak to be considered for backgroun modelling
#' @param bg_cutoff
#' Numeric between 0 and 1.
#' The peaks within the mass range with intensity below the \code{bg_cutoff}*100 quantile
#' are considered for background modelling. \code{bg_cutoff=1} keeps all peaks
#' and \code{bg_cutoff=0.5} would only keep the bottom half.
#' @param l_cutoff
#' Likelihood threshold or p-value. Peaks with a probability of being modelled as
#' background noise higher than this are filtered out.
#' @param int_index Column where to get intensities from matrix `x`
#' @return Peak matrix
#' @export
#'
peaks_local_bg = function(x, mass_range, bg_cutoff, l_cutoff,
                          int_index=2L, ...){
  # x = as.matrix(x)
  min_bg_peaks = 15
  masses = x[,1]
  intensities = x[, int_index]

  l_values = mapply(
    function(m, i) {
      neighbor_peaks = intensities[which(abs(m-masses)<mass_range)]
      n = length(neighbor_peaks)
      if (n < min_bg_peaks) {
        return(0)
      }
      int_cutoff = quantile(neighbor_peaks, bg_cutoff)
      bg = neighbor_peaks[neighbor_peaks <= int_cutoff]
      # Checks for special cases
      std_dev = sd(bg)
      if (is.na(std_dev)) {
        return(0)
      } else if (std_dev == 0) {
        if (i < int_cutoff) {
          # The current peak is equal to a constant background
          return(1)
        } else {
          # The current peak is above this constant noise level
          return(0)
        }
      }
      sd0 = sqrt((n-1)/n)*std_dev
      mx = mean(bg)
      l = pnorm(i, mean = mx, sd = sd0, lower.tail = F)
      return(l)
    },
    masses, intensities)

  # adj_pval = p.adjust(l_values[,1], 'holm')
  # l_values = cbind(l_values, adj_pval)
  peaks_mask = l_values < l_cutoff
  cn = c(colnames(x), 'bg_likelihood')
  x = cbind(x[peaks_mask, ], l_values[peaks_mask])
  colnames(x) = cn

  return(x)
}


#' Filter a peak list based on likelihood that peaks are above background local
#' noise
#'
#' @param s Spectra object
#' @param calc_frac_kept logical. Whether to calculate fraction of peaks kept
#' after peak local noise model filtering
#' @param ... arguments to [peaks_local_bg]
#'
#' @return
#' @export
#' @importFrom Spectra peaksData
#'
#' @examples
peakLocalNoiseModel = function(s, apply_queue=FALSE, calc_frac=FALSE,...) {

  if (!'centroided' %in% spectraVariables(s) | any(s$centroided == FALSE)) {
    stop('Spectra are in profile mode. peakLocalNoise only works on centroided spectra.')
  }

  if (calc_frac) {
    bef_lgs = lengths(s)
  }

  dots_args = list(...)

  pl_args = dots_args[c('mass_range', 'bg_cutoff', 'l_cutoff', 'int_index')]
  pl_args = pl_args[!sapply(pl_args, is.null)]
  pl_args[['object']] = s
  pl_args[['FUN']] = MALDIzooMS::peaks_local_bg
  s = do.call(addProcessing, pl_args)
  s$centroided = TRUE

  if (apply_queue) {
    ap_args = dots_args[c('in_memory', 'write_data', 'file', 'BPPARAM')]
    ap_args = ap_args[!sapply(ap_args, is.null)]
    ap_args[['s']] = s
    s = do.call(apply_preprocess, ap_args)
    if (calc_frac) aft_lgs = lengths(s)
    s$fracPeaksKept = aft_lgs/bef_lgs
  } else if (!apply_queue & calc_frac) {
    warning('To calculate the fraction of peaks left after apply peak local
             noise model, the processing queue will be applied')
    aft_lgs = lengths(s)
    s$fracPeaksKept = aft_lgs/bef_lgs
  }
  return(s)
}


#' Filter by SNR threshold
#'
#' @param x peaks matrix
#' @param snr signal-to-noise threshold
#' @param ...
#'
#' @return peaks matrix
#' @export
#'
#' @examples
filter_snr = function(x, snr, ...) {
  filt = x[['SNR']] > snr
  x = x[filt,]
  return(x)
}

#' Filter peaks in Spectra object by SNR
#'
#' @param s Spectra object
#' @param snr signal-to-noise threshold
#'
#' @return Spectra object
#' @export
#' @importFrom Spectra addProcessing
#'
#' @examples
filterSNR = function(s, snr) {
  s = addProcessing(s, filter_snr, snr=snr)
  s = apply_preprocess(
    s, in_memory=TRUE, write_data=FALSE)
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


#' Peak detection on Spectra object
#'
#' @param s Spectra object
#' @param apply_queue logical, whether to apply the processing queue
#' @param ... Arguments passed to [peak_detection] and [apply_preprocess]
#'
#' @return A Spectra object
#' @importFrom Spectra addProcessing
#' @export
#'
#' @examples
peakDetection = function(s, apply_queue=FALSE, ...) {
  dots_args = list(...)

  pd_args = dots_args[c('halfWindowSize', 'method',
                        'snr', 'k', 'descending', 'threshold',
                        'int_index', 'add_snr')]
  pd_args = pd_args[!sapply(pd_args, is.null)]
  pd_args[['object']] = s
  pd_args[['FUN']] = MALDIzooMS::peak_detection
  s = do.call(addProcessing, pd_args)
  # s = addProcessing(
  #   s, MALDIzooMS::peak_detection, ...)
  s$centroided = TRUE

  if (apply_queue) {
    ap_args = dots_args[c('in_memory', 'write_data', 'file', 'BPPARAM')]
    ap_args = ap_args[!sapply(ap_args, is.null)]
    ap_args[['s']] = s
    s = do.call(apply_preprocess, ap_args)
    # s = apply_preprocess(s, in_memory = TRUE, write_data=FALSE)
  }
  return(s)
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
  idx = MsCoreUtils::closest(x[,1], masses, duplicates='closest', ...)

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




