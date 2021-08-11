#' Read MALDI data in a given format in chunks and export in a different one
#'
#' @param indir
#' @param readf
#' @param outdir
#' @param writef
#' @param chunks
#'
#' @return
#' @export
#'
#' @examples
change_format_chunks = function(indir, readf, outdir, writef, nchunks = 80){
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
           fmt = ".tab"
           rw_f = rw_chunk_tsv
         },
         "mzML" = {
           fmt = ".mzML"
           rw_f = rw_chunk_mzml
         }
  )
  spectra_f = list.files(indir)

  # Filter out empty files
  filter_empty = lapply(
    spectra_f,
    check_empty
  )
  filter_empty = unlist(filter_empty)
  spectra_f = spectra_f[filter_empty]

  if (nchunks > 1) spectra_chunks = chunks(spectra_f, chunks)
  else spectra_chunks = list(spectra_f)


  invissible(lapply(
    spectra_chunks,
    rw_f,
    indir, read_f, outdir, fmt
  ))

}


#' Title
#'
#' @param x
#' @param indir
#' @param read_f
#' @param outdir
#' @param write_f
#' @param fmt
#'
#' @return
#' @importFrom MALDIquantForeign importMzMl
#' @export
#'
#' @examples
rw_chunk_mzml = function(x, indir, read_f, outdir, fmt) {
  infiles = file.path(indir, x)
  l = lapply(infiles, read_f)
  exportMzMl(l, path=outdir)
}

#' Title
#'
#' @param x
#' @param indir
#' @param read_f
#' @param outdir
#' @param write_f
#' @param fmt
#'
#' @return
#' @export
#'
#' @examples
rw_chunk_tsv = function(x, indir, read_f, outdir, fmt) {
  infiles = file.path(indir, x)
  l = lapply(infiles, read_f)

  outfiles = sub(pattern="\\.[[:alnum:]]+?$|(/|\\\\)+[^.\\\\/]+$",
                 replacement="", x=x)
  outfiles = paste0(outfiles, '.', fmt)
  outfiles = file.path(outdir, outfiles)

  exportTsv(l, path=outfiles)
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
  s = importMzMl(f)
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





