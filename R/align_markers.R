
#' Title
#'
#' @param markers
#' @param peaks
#' @param tolerance
#' @param match_tol
#' @param match_tol_unit
#' @param peaksby
#' @param aug_deam
#'
#' @return
#' @importFrom MALDIquant createMassPeaks binPeaks monoisotopicPeaks
#' @importFrom MALDIquant intensityMatrix match.closest plot
#' @importFrom dplyr filter
#' @importFrom tidyr %>%
#' @export
#'
#' @examples
pept_fly = function(markers, peaks, peaksby=NULL, tolerance=0.002,
                    match_tol = 0.6, match_tol_unit=c('Da', 'ppm'), aug_deam = F){

  match_tol_unit = match.arg(match_tol_unit)

  if (match_tol_unit == 'ppm'){
    match_tol = match_tol * markers$mass1/1e6
  }

  ref_peaks = createMassPeaks(
    mass = markers$mass1,
    intensity = rep(1, nrow(markers)))
  # int_ref_peaks = intRefPeaks(peaks, 'strict', 0.8, 0.002)

  if (is.null(peaksby)) {
    peaksby = factor(rep('fake', length(peaks)))
    markers[['fake']] = TRUE
  } else{
    if (length(peaksby) != length(peaks)) {
      stop("peaksby length and peaks must have the same length")
    }
    peaksby = factor(peaksby)
  }

  llp = sort(levels(peaksby))
  names(llp) = llp
  idxp <- lapply(llp, function(x)which(peaksby == x))

  counts_mat = matrix(0, ncol = length(llp), nrow = nrow(markers))
  colnames(counts_mat) = llp

  totals_mat = matrix(0, ncol = length(llp), nrow = nrow(markers))
  colnames(totals_mat) = llp

  # Align peaks
  al_peaks = custom_alignPeaks(
    peaks, reference=ref_peaks, tolerance = tolerance,
    allowNoMatches=T, emptyNoMatches=T)
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

  for (i in seq_along(llp)){
    j = idxp[[i]] # indices of peaks to select
    l = llp[i] # taxid that we are doing
    tpos = markers[[l]]
    counts = colSums(binMatrix[j, pos])
    counts[is.na(counts)] = 0
    counts_mat[tpos, l] = counts[tpos]
    totals_mat[tpos, l] = nrow(binMatrix[j,])
  }

  frac = rowSums(counts_mat) / rowSums(totals_mat)
  frac[is.na(frac)] = 0

  markers = markers %>%
    mutate(
      # ms2conf = pept_ms2conf & nhyd_ms2conf & nglut_ms2conf,
      malditof_mass = masses[pos],
      malditof_frac = frac) %>%
    mutate(error_mass = mass1 - malditof_mass) %>%
    mutate(error_ppm = error_mass/mass1 * 1e6)

  if (aug_deam) {
    msg_aug = paste0(
      'Setting fraction of different deamidated ',
      'versions of peptides to the maximum found'
    )
    message(msg_aug)
    markers = augment_deam(markers)
  }
  if (is.null(peaksby)) {
    markers = markers %>% select(-fake)
  }
  return(markers)

}



#' Title
#'
#' @param frac
#' @param ndeam
#'
#' @return
#'
#' @examples
.check_deam = function(frac, ndeam) {
  if (max(ndeam) == min(ndeam)) {
    return(frac)
  } else {
    return(max(frac))
  }
}

#' Title
#'
#' @param markers
#'
#' @return
#' @importFrom dplyr group_by mutate
#' @importFrom tidyr %>%
#' @export
#'
#' @examples
augment_deam = function(markers){
  m_by_seq = markers %>% group_by(seq, nhyd) %>%
    mutate(malditof_maxfrac = .check_deam(malditof_frac, ndeam)) %>%
    ungroup()
  return(m_by_seq)
}



#' Title
#'
#' @param ts
#' @param data
#' @param txlim
#' @param myby
#' @param gauss
#'
#' @return
#' @importFrom stats approx ksmooth
#' @export
#'
#' @examples
ccf_data = function(ts, data, txlim, myby, gauss) {

  #create an interpolation (isodists is accurate to 2 decimal places)

  #Generate the x values
  xout = seq(from = txlim[1], to = txlim[2], by = myby)
  #Resample against xout.
  yii = approx(x=data[,1], y=data[,2], xout=xout, method="linear", rule = 2)

  nmax = max(yii$y)

  yii$y = (yii$y-min(yii$y))/(nmax-min(yii$y))
  #yii$y = yii$y/max(yii$y)
  #Now resample the theoretical data:

  yri = data.frame(
    x = xout,
    y = rep(0, length(xout)),
    p = F
  )
  #go through each peak
  for(i in 1:length(ts$prob)){
    idx = which.min(abs(yri$x - ts$mass[i]))
    yri$y[idx] = ts$prob[i]
    yri$p[idx] = T
  }

  #Apply gaussian smoothing if set
  if(!is.na(gauss)){
    yrii = ksmooth(yri$x,yri$y,"normal",bandwidth = gauss)
    yrii$y = yrii$y / max(yrii$y)
    yri$y = yrii$y
  }

  return(data.frame(x=xout, yri=yri$y, yii=yii$y, p=yri$p))
}


#' Title
#'
#' @param txlim
#' @param laglim
#' @param rs_data
#' @param myby
#'
#' @return
#' @importFrom stats ccf
#' @export
#'
#' @examples
max_cc = function(rs_data, laglim, myby, mass1){

  mylagmax = laglim/myby

  ccd = ccf(rs_data$yri, rs_data$yii, plot=F, lag.max = mylagmax)

  cor = ccd$acf[,,1]
  lag = ccd$lag[,,1]

  max_idx = which.max(cor)
  lag = round(lag[max_idx] * myby, 3)
  err_ppm = lag/mass1 * 1e6

  out = data.frame(
    cor = round(cor[max_idx], 3),
    lag = lag,
    err_ppm = err_ppm
  )

  return(out)

}



#' Title
#'
#' @param s MassSpectrum object
#' @param markers markers data frame
#' @param myby
#' @param gauss
#' @param laglim
#' @param ...
#'
#' @return
#' @importFrom bacollite ms_subrange ms_iso
#' @export
#'
#' @examples
align_sample = function(s, markers, myby, gauss, laglim, ...){
  s = smoothIntensity(s, ...)
  s = as.matrix(s)
  markers = markers[, c('seq', 'ndeam', 'nhyd', 'pept_id')]
  markers_list = split(markers, seq(nrow(markers)))
  al = lapply(
    markers_list,
    function(p, s, myby, gauss, laglim){

      ts = ms_iso(p$seq, p$ndeam, p$nhyd)
      moff = 1.5
      lbl = min(ts$mass) - moff
      ubl = max(ts$mass) + moff
      myxlim = c(lbl, ubl)
      subms = ms_subrange(s, lbl, ubl)

      rs_data = ccf_data(ts, subms, myxlim, myby=myby, gauss=gauss)

      al_cor =  max_cc(rs_data, myby=myby, laglim=laglim, mass1=ts$mass[1])
      al_cor$pept_id = p$pept_id
      return(al_cor)
    },
    s, myby, gauss, laglim)
  return(al)
}



#' Title
#'
#' @param markers
#' @param data_path
#' @param ncores
#' @param metadata
#' @param nchunks
#' @param iocores
#' @param vch
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
get_max_cor = function(markers, data_path, metadata,
                       nchunks, ncores, iocores, vch = 5, ...){

  prepf = prepFun(align_sample, markers = markers_zooms, ...)
  spectra = paste0(metadata$sample_name, "_", metadata$replicate)
  cor_data = preprocessData(
    indir = data_path, readf = 'mzml', nchunks = nchunks, prepf = prepf,
    spectra = spectra, ncores = ncores, iocores = iocores, vch=vch)

  cor_data = mapply(
    function(x, n) {
      d = mapply(
        function(x, n){
          # x$pept = n
          x
        }, x, names(x), SIMPLIFY = FALSE, USE.NAMES = FALSE)
      d = do.call(bind_rows, d)
      # s = strsplit(strsplit(n, "\\.")[[1]][1], "_")[[1]]
      s = strsplit(n, "_")[[1]]
      d$sample_name = s[1]
      d$replicate = s[2]
      d
    }, cor_data, names(cor_data), SIMPLIFY = FALSE)
  cor_data = do.call(bind_rows, cor_data)
  return(cor_data)

}



#' Title
#'
#' @param s
#' @param markers
#' @param myby
#' @param gauss
#' @param laglim
#' @param ...
#'
#' @return
#' @importFrom bacollite ms_subrange ms_iso
#' @importFrom MALDIquant smoothIntensity
#' @export
#'
#' @examples
align_sample_full = function(s, markers, myby, gauss, laglim, ...) {
  # Smoothing
  s = smoothIntensity(s, ...)
  s = as.matrix(s)
  markers = markers[, c('seq', 'ndeam', 'nhyd', 'pept_id')]
  markers_list = split(markers, seq(nrow(markers)))
  al = lapply(
    markers_list,
    function(p, s, myby, gauss, laglim){

      ts = ms_iso(p$seq, p$ndeam, p$nhyd)
      moff = 1.5
      lbl = min(ts$mass) - moff
      ubl = max(ts$mass) + moff
      myxlim = c(lbl,ubl)
      subms = ms_subrange(s, lbl, ubl)

      rs_data = ccf_data(ts, subms, myxlim, myby=myby, gauss=gauss)
      rs_data$pept_id= p$pept_id

      al_cor =  max_cc(rs_data, myby=myby, laglim=laglim, mass1=ts$mass[1])
      al_cor$pept_id = p$pept_id

      rs_data$x_lag = rs_data$x + al_cor$lag

      return(list(cor = al_cor, al = rs_data))
    },
    s, myby, gauss, laglim)
  return(al)
}



#' Title
#'
#' @param markers
#' @param data_path
#' @param metadata
#' @param nchunks
#' @param ncores
#' @param iocores
#' @param vch
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
align_markers = function(markers, data_path, metadata,
                         nchunks, ncores, iocores, vch = 5, ...){

  prepf = prepFun(align_sample_full, markers = markers_zooms, ...)

  spectra = paste0(metadata$sample_name, "_", metadata$replicate)
  alignment_data = preprocessData(
    indir = data_path, readf = 'mzml', nchunks = nchunks, prepf = prepf,
    spectra = spectra, ncores = ncores, iocores = iocores, vch=vch)

  cor_data = mapply(
    function(x, n) {
      d = mapply(
        function(x, n){
          # x$cor$pept = n
          x$cor
        }, x, names(x), SIMPLIFY = FALSE, USE.NAMES = FALSE)
      d = do.call(bind_rows, d)
      # s = strsplit(strsplit(n, "\\.")[[1]][1], "_")[[1]]
      s = strsplit(n, "_")[[1]]
      d$sample_name = s[1]
      d$replicate = s[2]
      d
    }, alignment_data, names(alignment_data), SIMPLIFY = FALSE)
  cor_data = do.call(bind_rows, cor_data)

  xy_data = mapply(
    function(x, n) {
      d = mapply(
        function(x, n){
          # x$al$pept = n
          x$al
        }, x, names(x), SIMPLIFY = FALSE, USE.NAMES = FALSE)
      d = do.call(bind_rows, d)
      s = strsplit(strsplit(n, "\\.")[[1]][1], "_")[[1]]
      d$sample_name = s[1]
      d$replicate = s[2]
      d
    }, alignment_data, names(alignment_data), SIMPLIFY = FALSE)

  xy_data = do.call(bind_rows, xy_data)

  return(list(cor_data = cor_data, xy_data = xy_data))
}


