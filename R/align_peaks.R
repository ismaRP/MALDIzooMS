


#' Title
#'
#' @param m
#'
#' @return
#'
#' @examples
as.binary.matrix = function(m){
  stopifnot(is.matrix(m))
  isNA = which(is.na(m))
  m[] = 1L
  m[isNA] = 0L
  mode(m) = "integer"
  m
}


#' Title
#'
#' @param l
#'
#' @return
#'
#' @examples
as.matrix.MassObjectList = function(l){
  mass = unlist(lapply(l, function(x)x@mass),
                recursive=FALSE, use.names=FALSE)
  intensity = unlist(lapply(l, function(x)x@intensity),
                     recursive=FALSE, use.names=FALSE)
  uniqueMass = sort.int(unique(mass))
  n = lengths(l)
  r = rep.int(seq_along(l), n)

  i = findInterval(mass, uniqueMass)

  m = matrix(NA_real_, nrow=length(l), ncol=length(uniqueMass),
             dimnames=list(NULL, uniqueMass))
  m[cbind(r, i)] = intensity
  attr(m, "mass") = uniqueMass
  m
}


#' Title
#'
#' @param l
#' @param minFreq
#' @param method
#' @param tolerance
#' @param labels
#'
#' @return
#' @importFrom MALDIquant filterPeaks binPeaks createMassPeaks
#' @export
#'
#' @examples
intRefPeaks = function(l, method=c('strict', 'relaxed'), minFreq=0.9,
                       tolerance=0.002, labels=NULL){
  method = match.arg(method)
  if (is.null(labels)){
    referencePeaks = filterPeaks(
      binPeaks(l, method=method, tolerance=tolerance),
      minFrequency=minFreq
    )
  } else {
    referencePeaks = filterPeaks(
      binPeaks(l, method=method, tolerance=tolerance),
      minFrequency=minFreq, labels=labels
    )
  }


  m = as.binary.matrix(as.matrix.MassObjectList(referencePeaks))

  ## set peak intensity to number of occurrence
  intensity = unname(colMeans(m))

  createMassPeaks(mass=attr(m, "mass"), intensity=intensity)
}



#' Title
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
  pseq = parse.seq(sequence, max.missed.cleaves = 1, gpo_only=gpo_only)
  pseq$chain = chain
  pseq$species = species
  pseq
}


#' Title
#'
#' @param sequences
#' @param mc.cores
#' @param gpo_only
#'
#' @return
#' @importFrom dplyr bind_rows select
#' @importFrom parallel mcmapply
#' @importFrom tibble as_tibble
#' @export
#'
#' @examples
extRefPeaks = function(sequences, mc.cores=4L, gpo_only=F, masspeaks=T) {
  # Digest sequences
  peptides = do.call(
    mcmapply,
    c(list('FUN'=parse_seqs, 'mc.cores'=mc.cores, 'SIMPLIFY'=F,
           'USE.NAMES'=F, 'MoreArgs'=list('gpo_only'=gpo_only)),
      as.list(sequences)))

  peptides = do.call(bind_rows, peptides) %>% as_tibble()

  # peptides = peptides[peptides$nglut==0,]
  peptides = peptides[!duplicated(select(peptides, -species, -seqpos)),] %>%
    select(-species)

  if (masspeaks){
    return(createMassPeaks(mass=sort(peptides$mass1), intensity=rep(1, nrow(peptides))))
  } else {
    return(peptides)
  }

}



#' Title
#'
#' @param l
#' @param minFreq
#' @param tolerance
#' @param labels
#' @param th_peaks
#' @param ...
#'
#' @return
#' @importFrom MALDIquant determineWarpingFunctions warpMassPeaks binPeaks
#' @export
#'
#' @examples
custom_alignPeaks = function(l, tolerance, minFreq=NULL, reference=NULL, labels=NULL, ...){
  args_methods = list(...)

  if (is.null(args_methods$allowNoMatches)) args_methods$allowNoMatches = F
  if (is.null(args_methods$emptyNoMatches)) args_methods$emptyNoMatches = F

  if (is.null(reference)){
    if (is.null(minFreq)) minFreq=0.9
    reference = intRefPeaks(
      l,
      minFreq = minFreq,
      method = "strict",
      tolerance = tolerance,
      labels = labels
    )
  }


  warpingFunctions = determineWarpingFunctions(
    l,
    reference = reference,
    tolerance = tolerance,
    method='lowess',
    allowNoMatches = args_methods$allowNoMatches
  )
  l = warpMassPeaks(l, warpingFunctions, emptyNoMatches = args_methods$emptyNoMatches)
  # l = binPeaks(l, method = "strict", tolerance = tolerance)
  # l = binPeaks(l, method = "relaxed", tolerance = tolerance)

  return(l)
}


