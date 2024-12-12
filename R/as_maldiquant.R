
#' Transform list of peak matrices into a list of [MALDIquant::MassPeaks]
#'
#' @param pl List of peak matrices
#'
#' @return List of [MALDIquant::MassPeaks] objects
#' @export
#' @importFrom MALDIquant createMassPeaks
#'
asMassPeaksList = function(pl) {
  maldiquant_pl = lapply(
    pl, function(x) createMassPeaks(x[,1], x[,2])
  )
  names(maldiquant_pl) = names(pl)
  return(maldiquant_pl)
}
