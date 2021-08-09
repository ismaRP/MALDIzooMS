

importTsv = function(f) {
  s = read_tsv(f, col_names=c('mass', 'intensity'), col_types='dd')
  s = createMassSpectrum(
    mass=s[[1]],
    intensity=s[[2]],
    metaData=list(file=f)
  )
  return(s)
}


importTable = function(f){
  s = read.table(f)
  s = createMassSpectrum(
    mass=s[[1]],
    intensity=s[[2]],
    metaData=list(file=f)
  )
  return(s)
}

import_file.MzMl = function(f) {
  s = importMzMl(f)
  return(s[[1]])
}


chunks = function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

check_empty = function(x, indir){
  f = file.path(inndir, x)
  t = readLines(f,1)
  if (length(t)>0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


rw_chunk = function(x, indir, read_f, outdir, write_f, fmt) {
  infiles = file.path(indir, x)
  l = lapply(infiles, read_f)
  # Change extension
  # x = sub(pattern="\\.[[:alnum:]]+?$|(/|\\\\)+[^.\\\\/]+$",
  #         replacement="", x=x)
  # x = paste0(x, '.', fmt)
  # outfiles = file.path(opt$outdir, x)
  write_f(l, path=outdir)
}


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
change_format_chunks = function(indir, readf, outdir, writef, chunks = 80){
  switch(EXPR=readf,
         "readr" = {
           read_f = importTsv
         },
         "table" = {
           read_f = importTable
         },
         "mzml" = {
           read_f = import_file.MzMl
         }
  )
  switch(EXPR=writef,
         "tab" = {
           fmt = ".tab"
           write_f = exportTab
         },
         "txt" = {
           fmt = ".tab"
           write_f = exportTab
         },
         "mzml" = {
           fmt = ".mzML"
           write_f = exportImzMl
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

  spectra_chunks = chunks(spectra_f, chunks)


  invissible(lapply(
    spectra_chunks,
    rw_chunk,
    indir, readf, outdir, writef, fmt
  ))
}


#' Read MALDI data in a given format at once and export in a different one
#'
#' @param indir
#' @param readf
#' @param outdir
#' @param writef
#'
#' @return
#' @export
#'
#' @examples
change_format = function(indir, readf, outdir, writef){

  switch(EXPR=readf,
         "readr" = {
           fmt = ".tab"
           read_f = importTab
         },
         "table" = {
           fmt = ".tab"
           read_f = importTab
         },
         "mzml" = {
           fmt = ".mzML"
           read_f = importMzMl
         }
  )
  switch(EXPR=writef,
         "tab" = {
           fmt = ".tab"
           write_f = exportTab
         },
         "txt" = {
           fmt = ".tab"
           write_f = exportTab
         },
         "mzml" = {
           fmt = ".mzML"
           write_f = exportImzMl
         }
  )
  l = read_f(path=indir)
  write_f(l, path=outdir)

}
