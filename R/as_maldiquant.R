
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
    pl, function(x){
      mz = x[,1]
      int = x[,2]
      createMassPeaks(mz, int)
    })
  names(maldiquant_pl) = names(pl)
  return(maldiquant_pl)
}
