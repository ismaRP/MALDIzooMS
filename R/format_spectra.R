#' Read MALDI data in a given format in chunks and export in a different one
#'
#' @param indir
#' @param readf
#' @param nchunks
#' @param writef
#' @param spectra_names
#' @param mc.cores
#' @param outpath
#' @param in_fmt
#'
#' @return
#' @export
#'
#' @examples
change_format_chunks = function(spectra_names, indir, in_fmt, readf, outpath, writef, mc.cores=4, nchunks = 80){
  switch(EXPR=readf,
         "fread" = {
           read_f = importTsv
         },
         "table" = {
           read_f = importTable
         },
         "mzML" = {
           read_f = import_file.MzMl
         }
  )
  switch(EXPR=writef,
         "tab" = {
           fmt = "tab"
           rw_f = rw_chunk_tsv
         },
         "mzML" = {
           fmt = "mzML"
           rw_f = rw_chunk_mzml
         }
  )
  isFile = !isTRUE(file.info(outpath)$isdir)
  if (isFile & (nchunks > 1 || writef!='mzML')){
    stop("1 file write mode is only supported for mzML format and 1 chunk")
  }

  # spectra_f = list.files(indir)
  # Filter out empty files
  filter_empty = lapply(
    paste0(spectra_names, '.', in_fmt),
    check_empty,
    indir
  )
  filter_empty = unlist(filter_empty)
  spectra_names = spectra_names[filter_empty]

  if (nchunks > 1) spectra_chunks = chunks(spectra_names, nchunks)
  else spectra_chunks = list(spectra_names)

  invisible(mcmapply(
    rw_f,
    spectra_chunks,
    seq_along(spectra_chunks),
    MoreArgs=list(indir, read_f, outpath, in_fmt, fmt, nchunks), mc.cores=mc.cores
  ))

}


#' Title
#'
#' @param x sample names
#' @param indir
#' @param read_f
#' @param fmt
#' @param outpath
#' @param in_fmt
#'
#' @return
#' @importFrom MALDIquantForeign exportMzMl
#' @export
#'
#' @examples
rw_chunk_mzml = function(x, ch, indir, read_f, outpath, in_fmt, fmt, nchunks) {
  if (ch %% 5 == 0 | ch == nchunks){
    cat(sprintf('Chunk %i of %i', ch, nchunks), "\n")
  }
  # Create infiles
  infiles = paste0(x, '.', in_fmt)
  infiles = file.path(indir, infiles)
  l = lapply(infiles, read_f)
  l = mapply(
    function(s, n){
      s@metaData$id = n
      s
    }, l, x)

  # Write to outpath
  # outfiles = sub(pattern="\\.[[:alnum:]]+?$|(/|\\\\)+[^.\\\\/]+$",
  #                replacement="", x=x)
  isFile = !isTRUE(file.info(outpath)$isdir)
  if (!isFile) {
    outfiles = paste0(x, '.', fmt)
    outfiles = file.path(outpath, outfiles)
    invisible(mapply(exportMzMl, l, outfiles, MoreArgs = list(force=T)))
  } else {
    exportMzMl(l, path=outpath, force=T)
  }

}

#' Title
#'
#' @param x
#' @param indir
#' @param read_f
#' @param outdir
#' @param fmt
#'
#' @return
#' @export
#'
#' @examples
rw_chunk_tsv = function(x, ch, indir, read_f, outpath, in_fmt, fmt, nchunks) {
  if (ch %% 5 == 0 | ch == nchunks){
    cat(sprintf('Chunk %i of %i', ch, nchunks), "\n")
  }
  infiles = paste0(x, '.', in_fmt)
  infiles = file.path(infiles, in_fmt)
  l = lapply(infiles, read_f)
  l = mapply(
    function(s, n){
      s@metaData$id = n
      s
    }, l, x)

  # outfiles = sub(pattern="\\.[[:alnum:]]+?$|(/|\\\\)+[^.\\\\/]+$",
  #                replacement="", x=x)
  outfiles = paste0(x, '.', fmt)
  outfiles = file.path(outpath, outfiles)
  exportTsv(l, path=outfiles, force=T)
}

##### IMPORT FUNCTIONS READ FILE BY FILE

#' Title
#'
#' @param f
#'
#' @return
#' @importFrom MALDIquant createMassSpectrum
#' @importFrom data.table fread
#' @importFrom tibble tibble
#' @export
#'
importTsv = function(f) {
  s = tibble(fread(f, colClasses=c("numeric", "numeric"), sep="\t"))
  s = createMassSpectrum(
    mass=s[[1]],
    intensity=s[[2]],
    metaData=list(file=f)
  )
  return(s)
}


#' Title
#'
#' @param f
#'
#' @return
#' @importFrom MALDIquant createMassSpectrum
#' @importFrom utils read.table
#' @export
#'
importTable = function(f){
  s = read.table(f)
  s = createMassSpectrum(
    mass=s[[1]],
    intensity=s[[2]],
    metaData=list(file=f)
  )
  return(s)
}

#' Title
#'
#' @param f
#'
#' @return
#' @importFrom MALDIquantForeign importMzMl
#' @export
#'
#' @examples
import_file.MzMl = function(f) {
  s = importMzMl(f, verbose=F)
  return(s[[1]])
}


##### EXPORT FUNCTIONS WRITE LISTS OF FILES
#' Title
#'
#' @param l
#' @param path
#'
#' @return
#' @importFrom data.table fwrite
#' @export
#'
#' @examples
exportTsv = function(l, path) {
  mapply(
    function(x, f) fwrite(list(x@mass, x@instensity), f),
    l, path
  )
}



chunks = function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

check_empty = function(x, indir){
  f = file.path(indir, x)
  t = readLines(f,1)
  if (length(t)>0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}





