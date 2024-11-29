

#' Returns mzML header data for [mzR]
#'
#' @return A one row data.frame
#'
initialize_mzml_header = function() {
  header_df = data.frame(
    seqNum = 1,
    acquisitionNum = 0,
    msLevel = 1,
    polarity = -1,
    peaksCount = 0,
    totIonCurrent = 0,
    retentionTime = NaN,
    basePeakMZ = 0,
    basePeakIntensity = 0,
    collisionEnergy = NaN,
    ionisationEnergy = NaN,
    lowMZ = 0,
    highMZ = 0,
    precursorScanNum = NaN,
    precursorMZ = NaN,
    precursorCharge = NaN,
    precursorIntensity = NaN,
    mergedScan = NaN,
    mergedResultScanNum = NaN,
    mergedResultStartScanNum = NaN,
    mergedResultEndScanNum = NaN,
    injectionTime = NaN,
    filterString = '',
    spectrumId = '',
    centroided = FALSE,
    ionMobilityDriftTime = NaN,
    isolationWindowTargetMZ = NaN,
    isolationWindowLowerOffset = NaN,
    isolationWindowUpperOffset = NaN,
    scanWindowLowerLimit = NaN,
    scanWindowUpperLimit = NaN
  )
  return(header_df)
}


#' Change MALDI-TOF data format
#'
#' Read MALDI data in tsv (tab) or mzML format in chunks and export in tsv (tab) or mzML.
#' It is a wrapper interface around reading and writting functions.
#'
#' @param indir Path to data folder
#' @param readf Input format reading function. One of \code{'fread'}, \code{'table'}
#'              or \code{'mzML'}. \code{'fread'} uses [data.table::fread()],
#'              \code{'table'} uses [utils::read.table()] and \code{mzML} uses [mzR::openMSfile()].
#' @param nchunks Number of chunks to split the files into
#' @param writef Output format writing function. one of 'tab' or 'mzML'.
#'               Determines the destination format of the files.
#' @param spectra_names File names to transform (without extension)
#' @param mc.cores Number of cores to use
#' @param outpath Destination path
#' @param in_ext Extension of input files
#' @param sep Separator for tsv input files. Default is \code{"\t"}
#' @param verbose Print progress bar
#'
#' @return NULL
#' @importFrom mzR openMSfile
#' @export
#'
change_format_chunks = function(spectra_names, indir, in_ext, readf, outpath,
                                writef, sep='\t', mc.cores=4, nchunks = 80,
                                verbose=FALSE){
  switch(EXPR = readf,
         "fread" = {
           read_f = function(sep){
             function(x){
               import_tsv(f = x, sep = sep)
             }
           }
           read_f = read_f(sep = sep)
         },
         "table" = {
           read_f = function(sep){
             function(x){
               import_table(f = x, sep = sep)
             }
           }
           read_f = read_f(sep = sep)
         },
         "mzML" = {
           read_f = function(f) peaks(openMSfile(f))
         }
  )
  switch(EXPR = writef,
         "tab" = {
           rw_f = rw_chunk_tsv
         },
         "mzML" = {
           header_df = initialize_mzml_header()
           rw_f = rw_chunk_mzml
         }
  )

  isFile = !isTRUE(file.info(outpath)$isdir)
  if (isFile & (nchunks > 1 || writef != 'mzML')) {
    stop("1 file write mode is only supported for mzML format and 1 chunk")
  }

  # spectra_f = list.files(indir)
  # Filter out empty files
  # filter_empty = lapply(
  #   paste0(spectra_names, '.', in_ext),
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
    MoreArgs = list(indir = indir, read_f = read_f, outpath = outpath,
                    in_ext = in_ext, pb = pb),
    mc.cores = mc.cores
  ))

}


#'
#' Reads data in writes in mzML format
#'
#' Internal use only. Use [change_format_chunks] instead
#' @importFrom mzR writeMSData
#' @importFrom dplyr bind_rows
#' @importFrom fs is_file
#'
rw_chunk_mzml = function(x, ch, indir, read_f, outpath, in_ext, pb) {
  # Create infiles
  infiles = paste0(x, '.', in_ext)
  infiles = file.path(indir, infiles)
  l = lapply(infiles, read_f)

  header_template = initialize_mzml_header()
  # Write to outpath
  # outfiles = sub(pattern="\\.[[:alnum:]]+?$|(/|\\\\)+[^.\\\\/]+$",
  #                replacement="", x=x)
  isFile = !isTRUE(file.info(outpath)$isdir)
  if (!isFile) { # save multiple  files in outpath
    cat(sprintf('Saving spectra in multiple mzML files in %s folder', outpath))
    invisible(mapply(
      export_mzml, l, x, 1,
      MoreArgs = list(outpath = outpath, header_template = header_template)))
  } else {
    cat(sprintf('Saving all spectra into %s', outpath))
    idx = seq(1:length(l))
    headers_df = bind_rows(mapply(
      generate_header, l, x, idx, MoreArgs = list(header_template = header_template),
      SIMPLIFY = FALSE
    ))
    writeMSData(l, file = outpath, header = headers_df,
                backend = 'pwiz', outformat = 'mzml')
  }

  if (!is.null(pb)) {
    setTxtProgressBar(pb, ch)
  }
}


generate_header = function(x, id, idx, header_template){
  header_template$seqNum = idx
  header_template$peaksCount = nrow(x)
  header_template$lowMZ = min(x[,1])
  header_template$highMZ = max(x[,1])
  header_template$totIonCurrent = sum(x[,2])
  header_template$spectrumId = id
  max_idx = which.max(x[,2])
  header_template$basePeakMZ = x[max_idx,1]
  header_template$basePeakIntensity = x[max_idx,2]
  return(header_template)
}




#' Reads data in writes in mzML format
#'
#' Internal use only. Use [change_format_chunks] instead
rw_chunk_tsv = function(x, ch, indir, read_f, outpath, in_ext, pb) {
  infiles = paste0(x, '.', in_ext)
  infiles = file.path(infiles, in_ext)
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
  export_tsv(l, path = outfiles)
  if (!is.null(pb)) {
    setTxtProgressBar(pb, ch)
  }
}

##### IMPORT FUNCTIONS READ FILE BY FILE

#' Import data in tsv format
#'
#' Uses [data.table::fread]
#'
#' @param f Path to file
#' @param sep Separator. Default is \code{"\t"}
#'
#' @return Matrix with mz and intensities
#' @importFrom MALDIquant createMassSpectrum
#' @importFrom data.table fread
#' @importFrom tibble tibble
#' @export
#'
import_tsv = function(f, sep="\t") {
  s = as.matrix(
    fread(f, colClasses = c("numeric", "numeric"), sep = sep,
          col.names = c('mz', 'intensity')))
  return(s)
}


#' Import data in tsv format
#'
#' Uses [utils::read.table]
#'
#' @param f Path to file
#' @param sep Separator. Default is \code{''}, which for [utils::read.table()] can
#'            be one or more white spaces or tabs. If you're certain the separator
#'            is exactly \code{"\t"}, [import_tsv()] is a faster option.
#'
#' @return A data.frame with mz and intensity values
#' @importFrom MALDIquant createMassSpectrum
#' @importFrom utils read.table
#' @export
#'
import_table = function(f, sep = ''){
  s = read.table(f, col.names = c('mz', 'intensity'))
  return(s)
}

##### EXPORT FUNCTIONS WRITE LISTS OF FILES
#' Export list of spectra in multiple \code{"tsv"} files
#'
#' @param l List of spectra
#' @param path List of paths to write each spectra in \code{l}.
#'
#' @importFrom data.table fwrite
#' @export
#'
export_tsv = function(l, path) {
  invisible(mapply(
    function(x, f) fwrite(x, f, sep = "\t"),
    l, path
  ))
}

#' Export a matrix into mzML using [mzR::writeMSdata]
#'
#' @param x Matrix with mz and intensities
#' @param id Spectra ID
#' @param outpath Path to output folder
#' @param header_template Header template to complete with spectra data and id
#'
#' @export
#' @importFrom mzR writeMSData
#'
export_mzml = function(x, id, idx, outpath, header_template) {
  header = generate_header(x, id, idx, header_template)
  outfile = paste0(id, '.mzML')
  outfile = file.path(outpath, outfile)
  writeMSData(list(x), header = header, file = outfile,
              backend = 'pwiz', outformat = 'mzml')
}


chunks = function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))


#' Check is spectra file is empty
#'
#' @param x File name
#' @param indir Directory to data
#'
#' @return \code{TRUE} if file has any data, \code{FALSE} if it's empty
#' @export
#'
check_empty = function(x, indir){
  f = file.path(indir, x)
  return(file.size(f) == 0L)
}





