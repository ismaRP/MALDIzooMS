

#' Title
#'
#' @param l
#' @param labels
#' @param method
#' @param ignore.na
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
mergeReplicates <- function(l, labels, names, method=c("mean", "median", "sum"),
                           ignore.na=TRUE, ...) {

  ## test arguments
  MALDIquant:::.stopIfNotIsMassPeaksList(l)
  method <- match.arg(method)

  fun <- switch(method,
                "mean" = {
                  colMeans
                },
                "median" = {
                  .colMedians
                },
                "sum" = {
                  colSums
                }
  )

  .doByLabels(l=l, labels=labels, names=names, FUN=.mergeMassPeaks, fun=fun,
              ignore.na=ignore.na, ...)
}




#' Title
#'
#' @param l
#' @param labels
#' @param names
#' @param FUN
#' @param ...
#' @param mc.cores
#'
#' @return
#' @importFrom parallel mcmapply
#'
#' @examples
.doByLabels <- function(l, labels, names, FUN, ..., mc.cores=1L) {

  ## test parameters
  MALDIquant:::.stopIfNotIsMassObjectList(l)

  FUN <- match.fun(FUN)

  if (!missing(labels)) {
    ## drop unused levels and turn argument into factor
    if (is.factor(labels)) {
      labels <- droplevels(labels)
    } else {
      ## preserve order in labels
      labels <- factor(labels, levels=unique(labels))
    }

    if (missing(names)){
      names = lapply(split(as.vector(labels), labels),
        FUN=function(ll){
          nl = length(ll)
          return(paste(ll, 1:nl, sep = '_'))
        }
      )
    } else {
      if (!is.list(names)){
        names = split(names, labels)
      }
    }

    if (length(labels) != length(l)) {
      stop("For each item in ", sQuote("l"), " there must be a label in ",
           sQuote("labels"), "!")
    }

    ## replace tapply by split to preserve order
    tmp <- mcmapply(
      FUN=FUN,
      split(unlist(l), labels),
      split(as.vector(labels), labels),
      names,
      MoreArgs=list(...), mc.cores=mc.cores)

    k <- unlist(tmp)

    if (length(k) != length(tmp)) {
      k <- unsplit(tmp, labels)
    }
  } else {
    k <- FUN(l, ...)
  }

  k
}



#' Title
#'
#' @param l
#' @param ll
#' @param ln
#' @param fun
#' @param ignore.na
#'
#' @return
#' @importFrom MALDIquant createMassPeaks
#'
#' @examples
.mergeMassPeaks <- function(l, ll, ln,  fun=colMeans, ignore.na=TRUE) {

  fun <- match.fun(fun)

  ## create a matrix which could merged
  m <- MALDIquant:::.as.matrix.MassObjectList(l)

  mass <- attr(m, "mass")

  ## avoid named intensity/snr slot
  colnames(m) <- NULL

  isNA <- is.na(m)
  if (!ignore.na) {
    m[isNA] <- 0L
  }

  ## merge intensities
  intensity <- fun(m, na.rm=TRUE)

  ## merge snr
  for (i in seq_along(l)) {
    m[i, !isNA[i, ]] <- l[[i]]@snr
  }
  snr <- fun(m, na.rm=TRUE)

  ## merge metaData
  metaData <- .mergeMetaData(lapply(l, function(x)x@metaData), ll, ln)

  createMassPeaks(mass=mass, intensity=intensity, snr=snr, metaData=metaData)
}



#' Title
#'
#' @param m
#' @param ll
#' @param ln
#'
#' @return
#' @importFrom dplyr bind_rows
#'
#' @examples
.mergeMetaData <- function(m, ll, ln) {

  .flat <- function(x)unname(unlist(x))
  nm <- names(m[[1L]])
  names(nm) <- nm
  lapply(nm,
    function(n) {
      cur <- m[[1L]][[n]]
      allm <- lapply(m, "[[", n)
      len <- lengths(allm)
      if (all(sapply(allm, is.data.frame))){
        allm = mapply(
          function(x, lab, nam) {
            x$sample = lab
            x$replicate = nam
            x
          },
          allm, ll, ln, SIMPLIFY = F
        )
        return(do.call(bind_rows, allm))
      } else {
        if (!all(length(cur) == len) || !all(.flat(cur) == .flat(allm))) {
          if (!is.list(cur)) {
            all <- unlist(allm)
          }
          return(unname(allm))
        } else {
          return(cur)
        }
      }
    }
  )
}


