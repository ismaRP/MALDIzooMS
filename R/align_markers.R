
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
#' @importFrom tidyr %>%
#' @export
#'
#' @examples
pept_fly = function(markers, peaks, tolerance=0.002, match_tol = 0.6, aug_deam = F){

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

  if (aug_deam) {
    print('Augmenting fraction of deamidation versions of peptides')
    markers = augment_deam(markers)
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
.check_deam = function(frac, nglut) {
  if (max(nglut) == min(nglut)) {
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
    mutate(malditof_maxfrac = .check_deam(malditof_frac, nglut)) %>%
    ungroup()
  return(m_by_seq)
}


#' Title
#'
#' @param ts
#' @param data
#' @param txlim
#' @param laglim
#' @param myby
#' @param gauss
#' @param normlim
#' @param ccylim
#'
#' @return
#' @export
#'
#' @examples
ts_align = function(ts, data, txlim, laglim, myby, gauss=NA, normlim=NA, ccylim=c(-0.1,0.5)){

  #STEPSIZE & MAX LAG
  if(is.na(myby)) myby = 0.005 #125
  #mylagmax gives the 'reverse scaling' of the stepsize - useful when comparing etc.
  # mylagmax = 1/myby
  mylagmax = laglim/myby

  #create an interpolation (isodists is accurate to 2 decimal places)

  #Generate the x values
  xout = seq(from = txlim[1], to = txlim[2], by = myby)
  #Resample against xout.
  yii = approx(x=data[,1], y=data[,2], xout=xout, method="linear", rule = 2)
  #renormalise this segment
  #yii$y = yii$y/max(yii$y)

  nmax = max(yii$y)

  #If all the values are the same then 1: normalisation will fail and 2: no correlation - so return zeros:
  ry = range(yii$y)
  if(ry[1] == ry[2]){
    out = data.frame(
      cor = 0,
      lag = 0
    )
    return(out)
  }

  yii$y = (yii$y-min(yii$y))/(nmax-min(yii$y))

  #Now resample the theoretical data:
  # yri = approx(x=ts$mass,y=ts$prob,xout=xout, method="linear", rule = 2)
  # #set yvals to zero
  # yri$y[] =0
  yri = list(
    x = xout,
    y = rep(0, length(xout))
  )
  #go through each peak
  for(i in 1:length(ts$prob)){
    idx = which.min(abs(yri$x-ts$mass[i]))
    yri$y[idx] = ts$prob[i]
  }

  #Apply gaussian smoothing if set
  if(!is.na(gauss)){
    yrii = ksmooth(yri$x,yri$y,"normal",bandwidth = gauss)
    yrii$y = yrii$y / max(yrii$y)
    yri$y = yrii$y
  }

  ######################################################################################
  #Now we can do the cross-correlation:                                 4*mylagmax
  ccd = ccf(yri$y, yii$y, plot=F, lag.max = mylagmax)
  #message(sprintf("Length yr = %d, len yi = %d",length(yr),length(yi)))
  #plot(yr,yi,type="l",ylim=c(0,max(1000,max(yi)))  )


  cor = ccd$acf[,,1]
  lag = ccd$lag[,,1]

  max_idx = which.max(cor)

  out = data.frame(
    cor = res_max$cor,
    lag = res_max$lag * myby
  )

  return(out)

}


#' Title
#'
#' @param p
#' @param s
#'
#' @return
#' @importFrom bacollite ms_subrange ms_iso
#' @export
#'
#' @examples
align_pept = function(p, s){
  ts = ms_iso(p$seq, p$nglut, p$nhyd)
  moff = 1.5
  lbl = min(ts$mass) - moff
  ubl = max(ts$mass) + moff
  myxlim = c(lbl,ubl)

  subms = ms_subrange(s, lbl, ubl)

  al_cor =  ts_align(ts, subms, myxlim, gauss = 0.2, doplot=F)
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



