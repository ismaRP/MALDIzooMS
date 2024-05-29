

initialize_mzml_header = function() {
  header_df = data.frame(
    seqNum=1,
    acquisitionNum=0,
    msLevel=1,
    polarity=-1,
    peaksCount=0,
    totIonCurrent=0,
    retentionTime=NaN,
    basePeakMZ=0,
    basePeakIntensity=0,
    collisionEnergy=NaN,
    ionisationEnergy=NaN,
    lowMZ=0,
    highMZ=0,
    precursorScanNum=NaN,
    precursorMZ=NaN,
    precursorCharge=NaN,
    precursorIntensity=NaN,
    mergedScan=NaN,
    mergedResultScanNum=NaN,
    mergedResultStartScanNum=NaN,
    mergedResultEndScanNum=NaN,
    injectionTime=NaN,
    filterString='',
    spectrumId='',
    centroided=FALSE,
    ionMobilityDriftTime=NaN,
    isolationWindowTargetMZ=NaN,
    isolationWindowLowerOffset=NaN,
    isolationWindowUpperOffset=NaN,
    scanWindowLowerLimit=NaN,
    scanWindowUpperLimit=NaN
  )
  return(header_df)
}


#' Read MALDI data in a given format in chunks and export in a different one
#'
#' @param indir
#' @param readf One of 'fread', 'table' or 'mzML'
#' @param nchunks
#' @param writef One of 'tab' or 'mzML'
#' @param spectra_names
#' @param mc.cores
#' @param outpath
#' @param in_fmt
#'
#' @return
#' @export
#'
#' @examples
change_format_chunks = function(spectra_names, indir, in_fmt, readf, outpath,
                                writef, sep='\t', mc.cores=4, nchunks = 80,
                                verbose=NULL){
  switch(EXPR=readf,
         "fread" = {
           read_f = function(sep){
             function(x){
               importTsv(f=x, sep=sep)
             }
           }
           read_f = read_f(sep=sep)
         },
         "table" = {
           read_f = function(sep){
             function(x){
               importTable(f=x, sep=sep)
             }
           }
           read_f = read_f(sep=sep)
         },
         "mzML" = {
           read_f = function(f) peaks(openMSfile(f))
         }
  )
  switch(EXPR=writef,
         "tab" = {
           rw_f = rw_chunk_tsv
         },
         "mzML" = {
           header_df = initialize_mzml_header()
           rw_f = rw_chunk_mzml
         }
  )

  isFile = !isTRUE(file.info(outpath)$isdir)
  if (isFile & (nchunks > 1 || writef!='mzML')){
    stop("1 file write mode is only supported for mzML format and 1 chunk")
  }

  # spectra_f = list.files(indir)
  # Filter out empty files
  # filter_empty = lapply(
  #   paste0(spectra_names, '.', in_fmt),
  #   check_empty,
  #   indir
  # )
  # filter_empty = unlist(filter_empty)
  # spectra_names = spectra_names[filter_empty]

  if (nchunks > 1) spectra_chunks = chunks(spectra_names, nchunks)
  else spectra_chunks = list(spectra_names)

  if (verbose) {
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = nchunks, # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar
  } else {
    pb = NULL
  }

  invisible(mcmapply(
    rw_f,
    spectra_chunks,
    seq_along(spectra_chunks),
    MoreArgs=list(indir=indir, read_f=read_f, outpath=outpath, in_fmt=in_fmt, pb=pb),
    mc.cores=mc.cores
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
#' @importFrom mzR writeMSData
#' @importFrom dplyr bind_rows
#' @importFrom fs is_file
#' @export
#'
#' @examples
rw_chunk_mzml = function(x, ch, indir, read_f, outpath, in_fmt, pb) {
  # Create infiles
  infiles = paste0(x, '.', in_fmt)
  infiles = file.path(indir, infiles)
  l = lapply(infiles, read_f)

  header_template = initialize_mzml_header()
  # Write to outpath
  # outfiles = sub(pattern="\\.[[:alnum:]]+?$|(/|\\\\)+[^.\\\\/]+$",
  #                replacement="", x=x)
  isFile = fs::is_file(outpath)
  if (!isFile) { # save multiple  files in outpath
    invisible(mapply(
      export_mzR, l, x,
      MoreArgs = list(outpath=outpath, header_template=header_template)))
  } else {
    headers_df = bind_rows(mapply(
      generate_header, l, x, MoreArgs=list(header_template = header_template)
    ))
    writeMSData(l, file=outpath, header=headers_df,
                backend='pwiz', ouformat='mzml')
  }

  if (!is.null(pb)){
    setTxtProgressBar(pb, ch)
  }
}


generate_header = function(x, id, header_template){
  header_template$peaksCount = nrow(x)
  header_template$lowMZ = min(x[,1])
  header_template$highMZ = max(x[,1])
  header_template$totIonCurrent = sum(x[,2])
  header_template$spectrumId = id
  return(header_template)
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
rw_chunk_tsv = function(x, ch, indir, read_f, outpath, in_fmt, pb) {
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
  outfiles = paste0(x, '.tab')
  outfiles = file.path(outpath, outfiles)
  exportTsv(l, path=outfiles)
  if (!is.null(pb)){
    setTxtProgressBar(pb, ch)
  }
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
importTsv = function(f, sep="\t") {
  s = as.matrix(
    fread(f, colClasses=c("numeric", "numeric"), sep=sep,
          col.names=c('mz', 'intensity')))
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
importTable = function(f, sep=''){
  s = read.table(f, col.names = c('mz', 'intensity'))
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
  invisible(mapply(
    function(x, f) fwrite(x, f, sep="\t"),
    l, path
  ))
}

#' Export a matrix into mzML using mzR
#'
#' @param x Matrix with mz and intensities
#' @param id Spectra ID
#' @param template_header
#'
#' @return
#' @export
#' @importFrom mzR writeMSData
#'
#' @examples
export_mzR = function(x, id, outpath, header_template) {
  header = generate_header(x, id, header_template)
  outfile = paste0(id, '.mzML')
  outfile = file.path(outpath, outfile)
  writeMSData(list(x), header=header, file=outfile,
              backend='pwiz', outformat='mzml')
}


chunks = function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))


#' Title
#'
#' @param x
#' @param indir
#'
#' @return
#' @export
#'
#' @examples
check_empty = function(x, indir){
  f = file.path(indir, x)
  t = readLines(f,1)
  if (length(t)>0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}





