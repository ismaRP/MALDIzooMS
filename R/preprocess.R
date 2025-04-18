

#' Print message only if verbose
#'
#' @param msg Message to print
#' @param verbose logical, whether to print the message or not
#'
#' @export
#'
#' @examples
#' verbose = TRUE
#' print_progress('', verbose)
print_progress = function(msg, verbose) {
  if (verbose) message(msg)
}

#' Plot Spectra in different preprocessing stages
#'
#' @param l List of Spectra objects. Each Spectra is at a different state in the
#' preprocessing.
#' @param p Spectra object with peaks
#' @param scale_factor Vector of scale factors used to normalize peak with respect
#' to the rest of the spectra.
#' @param mzrange mz range to limit plot
#' @param prep_labels (optional) character vector with labels for each preprocessing
#'                    step to be plotted
#' @param n_col Number of columns in the plot
#'
#' @export
#' @importFrom Spectra concatenateSpectra
#'
plot_preprocess = function(l, p=NULL, scale_factor=1, mzrange=c(800, 4000),
                           prep_labels=c('Raw', 'Smoothed', 'Baseline corrected', 'Peaks'),
                           n_col=2) {
  # l = unlist(unname(list(unname(l))))
  l = sapply(l, apply_preprocess, write_data = FALSE, USE.NAMES = FALSE)

  spectra_lgs = lengths(l)
  spec_ids = lapply(l, function(x) x$spectrumId)
  if (!is.null(p)) {
    p = apply_preprocess(p, write_data = FALSE)
    spectra_lgs = c(spectra_lgs, length(p))
    spec_ids[[length(spec_ids)+1]] = p$spectrumId
    if (any(p$normalized == TRUE) | any(peaksData(p[1])[[1]][,2] < 0)) {
      # De-normalize to put in the same scale
      p$invScaleFactor = 1/scale_factor
      p = applyProcessing(normalizeSpectra(p, scale_f='invScaleFactor'))
    }
  }
  # Check lengths are the same
  if (!all(spectra_lgs == spectra_lgs[1])) {
    stop('All Spectra objects must be the same length')
  }
  # Check they have the same spectrumIds
  if (!all(sapply(spec_ids[-1], function(x) all(x==spec_ids[[1]])))) {
    stop("All Spectra objects must contain the same spectrumId's")
  }

  colors_prep = palette.colors(n=length(l)+1, palette = 'R4')
  colors_prep = colors_prep[2:(length(l)+1)]
  n_row = ceiling(spectra_lgs[1]/n_col)
  par(mfrow=c(n_row, n_col), mar=c(4, 4, 2, 1))
  for (i in 1:length(l[[1]])) {
    prep_i = Spectra::concatenateSpectra(sapply(l, '[', i))
    prep_i = filterMzRange(prep_i, mz=mzrange)

    plotSpectraOverlay(
      prep_i, col=alpha(colors_prep, 0.7),
      main=prep_i$spectrumId[1], type='l')
    if (!is.null(p)) plotSpectra(p[i], add=T, col='black')
  }
  legend('topright', legend=prep_labels, col=c(colors_prep, 'black'),
         lty=1, lwd=2)
}

#' Quality control spectra from [MALDIrppa::screenSpectra()]
#' Works on peak matrix rather than on [MALDIquant::MassSpectrum]
#' @param x Raw peak matrix with mz and intensities
#' @param ... Currently not in use, required by [Spectra::addProcessing,Spectra-method].
#'
#' @return `numeric`, Atipicality score
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
#' peak matrix. It can calculate calculate lower and upper A limits by label.
#'
#' @param s Spectra object
#' @param BPPARAM Parallel computing configuration from [BiocParallel]. Default is
#'        is [BiocParallel::bpparam()]
#' @param labels factor to do groupwise calculation of A-score limits for flagging
#'        spectra.
#'
#' @return A list with [Spectra::Spectra()] object with \code{QCflag} spectrum
#'         variable added and the atipicality limits for flagging.
#' @export
#' @importFrom Spectra addProcessing
#' @importFrom BiocParallel register bpparam
atipicalitySpectra = function(s, labels=NULL, BPPARAM = bpparam()){
  qcA = addProcessing(s, atipicality_spectra)

  register(BPPARAM)

  qcA = unlist(peaksData(qcA, BPPARAM=BPPARAM))

  s$Atipicality = qcA

  if (!is.null(labels)) {
    QCmin = vector('numeric', length(s))
    QCmax = vector('numeric', length(s))
    if (!is.factor(labels)) labels = factor(labels)
    for (l in levels(labels)) {
      idx = which(labels == l)
      qcl = qc_limits(s$Atipicality[idx])
      QCmin[idx] = qcl[1]
      QCmax[idx] = qcl[2]
    }
  } else {
    QClimits = qc_limits(s$Atipicality)
    QCmin = QClimits[1]
    QCmax = QClimits[2]
  }

  s$QCflag = 'QCpass'
  s$QCmin = QCmin
  s$QCmax = QCmax
  qcpass = sapply(
    seq_along(s$Atipicality),
    function(i) {
      (s$QCmin[i] < s$Atipicality[i] & s$Atipicality[i] < s$QCmax[i])
    })
  s$QCflag[!qcpass] = 'QCfail'

  s = reset(s)

  return(s)
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
  QClimits = adjboxStats(A, coef = threshold)$fence
  return(QClimits)
}


#' Applies processing queue to Spectra
#'
#' @param s Spectra object
#' @param in_memory logical. `in_memory` is `TRUE` keeps the peaks data in memory using
#'        a [Spectra::MsBackendDataFrame]. If `FALSE`, peaks data are kept on disk
#'        un a MzML file using the [Spectra::MsBackendMzR]. If `FALSE`, `write_data` is
#'        overridden and the spectra is written to file.
#' @param write_data logical. Whether to write spectra data into a file
#'        This will apply processing, export spectra into mzML file(s) and
#'        reload into a new Spectra object with MsBackendMzR. This is useful
#'        if we want to apply preprocessing and the resulting object is still
#'        too big to keep in memory, either because there are too many spectra
#'        or because we still haven't picked peaks.
#' @param file If `write_data` is `TRUE`, or `in_memory` is `FALSE`,
#'        required file or files to write spectra.
#'        IF there is one file, all spectra is saved on the same file.
#'        If multiple files, there needs to be one per spectra.
#' @param ncores Number of cores for parallel processing of the newly created
#'        Spectra object. It can also be specified using ncores in the metadata
#'        slot of `s`. If both are specified, the value of the function argument
#'        takes precedence.
#'        If both are `NULL`, it is set to the maximum available cores - 2
#'        with `parallel::detectCores() - 2`. See [parallel::detectCores()]
#' @param verbose Whether to print progress and progress bars.
#' @param BPPARAM bpparam for the newly created Spectra object. Default is
#'        [BiocParallel::MulticoreParam] with `ncores`.
#'
#' @return A Spectra object with either [Spectra::MsBackendMzR] (`in_memory=FALSE`) or
#'         [Spectra::MsBackendMemory] (`ìn_memory=TRUE`) backends
#' @importFrom Spectra export Spectra MsBackendMzR MsBackendDataFrame
#' @importFrom BiocParallel MulticoreParam multicoreWorkers SerialParam SnowParam
#' @importFrom BiocParallel register
#' @export
#'
apply_preprocess = function(
    s, in_memory=TRUE,
    write_data=FALSE,
    file=NULL,
    ncores=NULL,
    verbose=FALSE,
    BPPARAM=NULL) {

  if (is.null(BPPARAM)) {
    if (is.null(ncores)) {
      ncores = multicoreWorkers()
    }
    if (ncores == 1) {
      BPPARAM = SerialParam(progressbar = verbose)
    } else if (.Platform$OS.type == "windows") {
      BPPARAM = SnowParam(workers=ncores, progressbar = verbose)
    } else {
      BPPARAM = MulticoreParam(workers=ncores, progressbar = verbose)
    }
  }
  register(BPPARAM)
  pcs = processingChunkSize(s)
  if (in_memory) { # Apply processing and keep in memory
    print_progress('Running processing queue...', verbose)
    peaks = peaksData(s, columns=peaksVariables(s), BPPARAM=BPPARAM)
    # peaks = peaksData(s)
    spd = spectraData(s)
    peaks_vars = colnames(peaks[[1]])
    for (i in 1:ncol(peaks[[1]])) {
      spd[[peaks_vars[i]]] = lapply(peaks, '[', , i)
    }
    s = Spectra(
      spd, backend = MsBackendMemory(),
      peaksVariables = peaks_vars, metadata=s@metadata,
      BPPARAM = BPPARAM)
    s$peaksCount = count_peaks(s)
    processingChunkSize(s) = pcs

    if (write_data) {
      print_progress('Exporting mzml files...', verbose=verbose)
      export(object=s, backend=MsBackendMzR(), format='mzML',
             file=file, BPPARAM=BPPARAM)
    }
  } else { # Keep data on disk and write to file
    if (write_data == FALSE) {
      warning('For in_memory=FALSE data needs to be written to file')
    }
    if (is.null(file)) {
      stop(
        paste0('For in_memory=FALSE a file to keep the spectra on disk needs',
        'to be provided'))
    }
    # Exporting runs the processing queue
    print_progress('Running processing queue and exportin mzml files...',
                   verbose=verbose)
    export(object=s, backend=MsBackendMzR(), format='mzML',
           file=file, BPPARAM=BPPARAM)
    # Reload the file into a MsBackendMzR
    s = suppressMessages(
      Spectra(file, source = MsBackendMzR(), centroided = FALSE,
              BPPARAM = MulticoreParam(workers = ncores)))
    processingChunkSize(s) = pcs
  }
  return(s)

}


#' Concatenation of Spectra
#'
#' Prior concatenation, it applies the processing queue and changes to a
#' MsBackendDataFrame
#'
#' @param x List of [Spectra::Spectra] object
#' @param ... Parameters for [apply_preprocess]
#'
#' @return A single [Spectra::Spectra] object
#' @export
#'
concatenate_spectra = function(x, ...) {
  x = sapply(x, apply_preprocess, ...)
  x = concatenateSpectra(x)
  return(x)
}





#' Baseline correction subtraction
#'
#' @param x Peaks matrix
#' @param int_index Index of the intensity to calculate the baseline
#' @param keep_bl Keep basseline in a separate column
#' @param substract_index Index of the intensity column to be substracted.
#'        If NULL, there is no substraction
#' @param in_place Replace the intensity in \code{substract_index} with the baseline
#'        substracted intensity
#' @param ... Arguments passed to [MsCoreUtils::estimateBaseline()]. Currently `method`.
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
    x[,1], x[, int_index], ...)
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


#' Apply baseline correction to Spectra object
#'
#' @param s Spectra object
#' @param ... Arguments passed to [baseline_correction]
#'
#' @return Spectra object
#' @export
#'
baselineCorrection = function(s, ...) {
  s = addProcessing(
    s, baseline_correction, ...)
  return(s)
}


#' Intensity smoothing spectrapeak matrix
#'
#' @param x Peak matrix
#' @param method Smoothing method. One of \code{'MovingAverage'}, \code{'WeightedMovingAverage'}
#' or \code{'SavitzkyGolay'}
#' @param hws Half window size
#' @param k For Savitzky-Golay filter, this is the polynomial order for the coefficients
#' @param int_index Column with the intensity to smooth
#' @param in_place Whether smoothing is in place, or a new column with the smoothened
#'                 intensity is added
#' @param ... Parameters passed to other methods. Not in use, required by [Spectra::addProcessing,Spectra-method].
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


#' Smooth spectra
#'
#' @param s [Spectra::Spectra] object
#' @param ... Parameters passed to [smooth]
#'
#' @return A [Spectra::Spectra] object
#' @export
#'
smoothSpectra = function(s, ...) {
  s = addProcessing(
    s, MALDIzooMS::smooth, ...)
  s$smoothed = TRUE
  return(s)
}




#' Normalize peaks matrix
#'
#' @param x Peaks matrix
#' @param func Function to be applied to the intensities
#' @param scaleFactor Scaling factor
#' @param ... Currently not in use, required by [Spectra::addProcessing,Spectra-method].
#'
#' @return Peaks matrix
#' @export
#'
normalize_spectra = function(x, scaleFactor=1, func=NULL, ...) {
  if (is.null(func)) {
    x[,2] = x[,2] / scaleFactor
  } else {
    x[,2] = x[,2] / func(x[,2])
  }
  return(x)
}



#' Normalize spectra object
#'
#' @param s [Spectra::Spectra] object
#' @param scale_f spectra variable to use as normalization factor or a function
#' to be applied to spectra intensities
#'
#' @return [Spectra::Spectra] object
#' @export
#'
normalizeSpectra = function(s, scale_f) {
  if (is.function(scale_f)) {
    s = addProcessing(s, normalize_spectra, func=scale_f)
  } else if (is.character(scale_f) & scale_f %in% colnames(s@backend@spectraData)) {
    s$scaleFactor = s[[scale_f]]
    s = addProcessing(s, normalize_spectra, spectraVariables='scaleFactor')
  } else {
    stop('"scale_f" must be either a spectra variable or a function')
  }
  s$normalized = TRUE
  return(s)
}

