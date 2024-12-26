

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
#' @param spectra_in File names or paths to transform (with or without extension)
#' @param indir Path to data folder. If `NULL`, it assumes `spectra_in` contains whole
#'              paths to files.
#' @param readf Input format reading function. One of \code{'fread'}, \code{'table'}
#'              or \code{'mzML'}. \code{'fread'} uses [data.table::fread()],
#'              \code{'table'} uses [utils::read.table()] and \code{mzML} uses [mzR::openMSfile()].
#' @param spectra_out Names of the spectra to go in the output. These can either constitute
#'                    filenames and/or spectrumId in a mzML file.
#' @param outpath Destination path
#' @param nchunks Number of chunks to split the files into
#' @param writef Output format writing function. one of 'tab' or 'mzML'.
#'               Determines the destination format of the files.
#'
#' @param ncores Number of cores to use
#'
#' @param in_ext Extension of input files
#' @param sep Separator for tsv input files. Default is \code{"\t"}
#' @param verbose Print progress bar
#'
#' @return NULL
#' @importFrom mzR openMSfile
#' @importFrom BiocParallel MulticoreParam SnowParam SerialParam bpmapply bplapply bpoptions
#' @importFrom fs is_file is_dir path_file path_dir file_exists path
#' @export
#'
change_format_chunks = function(spectra_in, spectra_out=NULL, outpath=NULL, readf=NULL,
                                indir=NULL, in_ext=NULL,
                                writef=NULL, sep='\t', ncores=4, nchunks = 80,
                                verbose=FALSE){
  if (is.null(readf)) stop('readf must be defined')
  switch(EXPR = readf,
         "fread" = {
           read_f = function(sep){
             function(x){
               import_tsv(f = x, sep = sep)
             }
           }
           read_f = read_f(sep = sep)
           auto_in_ext = '.csv'
         },
         "table" = {
           read_f = function(sep){
             function(x){
               import_table(f = x, sep = sep)
             }
           }
           read_f = read_f(sep = sep)
           auto_in_ext = '.txt'
         },
         "mzML" = {
           read_f = function(f) peaks(openMSfile(f))
           auto_in_ext = '.mzML'
         },
         {stop("The read function readf must be 'fread', 'table' or 'mzML'")}
  )
  if (is.null(writef)) stop('writef must be defined')
  switch(EXPR = writef,
         "tab" = {
           rw_f = rw_chunk_tsv
         },
         "mzML" = {
           header_df = initialize_mzml_header()
           rw_f = rw_chunk_mzml
         },
         {stop("The write function writef must be either 'tab' or 'mzML'")}
  )

  # From the given input arguments, we will create:
  #   - in_names: Spectra names without file extension
  #   - in_files : paths to each spectra file
  if (is.null(indir) & (all(file_exists(spectra_in)) & length(spectra_in) > 1)) {
    # We assume spectra_in contains full paths
    if (any(is_dir(spectra_in))) {
      stop(
        paste0(
          'Either provide full paths in spectra_in, or provide indir with a path\n',
          'and indir with an path to the file spectra')
      )
    }
    in_files = spectra_in
    in_names = path_file(spectra_in)
    in_names = sapply(strsplit(in_names, '\\.'), '[', 1, USE.NAMES = F)
  } else if (is.null(indir) & any(!file_exists(spectra_in))) {
    stop("Files in spectra_in don't exists, did you mean to provide indir path?")
  } else if (length(indir) == 1 & is_dir(indir)) {
    # We assume spectra_in contains spectra names
    if (any(grepl('\\..+?$', spectra_in))) {
      # If spectra_in have extensions, remove them, extract name and plug back
      ext_files = sapply(
        strsplit(spectra_in, '\\.'),
        function(x) {
          if (length(x)>1) x[length(x)] else ''
        })
      ext_files = paste0('.', ext_files)
      in_names = path_file(spectra_in)
      path_spec = path_dir(spectra_in)
      in_names = sapply(strsplit(in_names, '\\.'), '[', 1, USE.NAMES = F)
      in_files = path(indir, path_spec, paste0(in_names, ext_files))
    } else {
      if (is.null(in_ext)) {
        cat(sprintf('Using %s as input files extension', auto_in_ext))
        in_ext = auto_in_ext
      }
      in_names = path_file(spectra_in)
      in_files = path(indir, paste0(spectra_in, '.', in_ext))
    }
  } else if (length(spectra_in) == 1 & file_exists(spectra_in)) {
    # We assume a single spectra file is given
  } else {
    stop('Please provide indir path and or spectra_in')
  }


  if (all(!file_exists(in_files))) {
    print(in_files)
    stop('None of the files exists, are in_ext and indir right?')
  } else if (any(!file_exists(in_files))) {
    stop("Some of the files don't exist, are in_ext and indir right?")
  }

  if (length(outpath) > 1) {
    stop(
      'Only use outpath for a general directory.\n',
      'If you need to specify different directories per file, use spectra_out\n')
  }

  # From the given input arguments, we will create:
  #   - out_names: Spectra names without file extension
  #   - out_files: path/s to each spectra
  if (is.null(outpath)) {
    message('Multiple files, 1 per spectra.\nUsing spectra_out paths.\n')
    # Multiple output files, 1 per spectra
    many_to_one=FALSE
    if (!is.null(spectra_out)) {
      # We assume spectra_out contains full paths
      out_files = spectra_out
      out_names = path_file(spectra_out)
      if (any(grepl('\\..+?$', out_names))) {
        out_names = sapply(strsplit(out_names, '\\.'), '[[', 1, USE.NAMES = F)
        add_ext = FALSE
      } else {
        add_ext = TRUE
      }
    } else {
      stop('Provide output files either in spectra_out or outpath')
    }
  } else if (is_dir(outpath)) {
    message('Writing one spectra per file.\n')
    # Multiple output files, 1 per spectra
    many_to_one=FALSE
    ext_files = ''
    add_ext = TRUE # Will ask writing function to add its own extension to files
    if (is.null(spectra_out)) {
      # Get spectra names from spectra_in
      message('Using spectra_in names.\n')
      out_names = in_names
    } else if (any(grepl('/', spectra_out)) | any(grepl('\\..+?$', spectra_out))) {
      message('Joining outpath with spectra_out paths.\n')
      out_dir = path_dir(spectra_out)
      outpath = path(outpath, out_dir)
      out_names = path_file(spectra_out)
      out_names = sapply(strsplit(out_names, '\\.'), '[[', 1, USE.NAMES = F)
      ext_files = sapply(
        strsplit(spectra_out, '\\.'),
        function(x) if (length(x)>1) x[length(x)] else ''
      )
      ext_files = paste0('.', ext_files)
      add_ext = FALSE # outdir extensions override writing function extensions
    } else {
      out_names = spectra_out
    }
    out_files = path(outpath, paste0(out_names, ext_files))
  } else { # outpath is assumed a file. Check if folder exists?
    # 1 file with multiple spectra
    many_to_one = TRUE
    parent_f = path_dir(outpath)
    if (!dir_exists(parent_f)){
      stop(sptrinf('%s does not exist', parent_f))
    }
    message(sprintf('Saving all spectra in a single file %s\n', outpath))
    add_ext = FALSE
    if (!grepl('\\..+?$', outpath)) {
      warning('Output file has no extension. It will be automatically added based',
              ' on writing format function.')
      add_ext = TRUE
    }

    if (is.null(spectra_out)) {
      out_names = in_names
    } else {
      out_names = spectra_out
    }
    out_files = outpath
  }

  if (many_to_one & (nchunks > 1 | writef != 'mzML')) {
    stop('Aggregating all spectra into one file is only available for ',
         'nchunks = 1 and writef = "mzML"')
  }

  if (is.null(ncores)) {
    if (many_to_one) {
      ncores=1
    } else {
      ncores = detectCores() - 2
    }
  }

  if (ncores == 1) {
    BPPARAM = SerialParam(progressbar=FALSE)
  } else if (.Platform$OS.type == "windows") {
    BPPARAM = SnowParam(workers=ncores, progressbar=FALSE)
  } else {
    BPPARAM = MulticoreParam(workers=ncores, progressbar=FALSE)
  }

  # Get processing chunks
  if (nchunks > 1) {
    spectra_chunks = chunks(seq_along(in_names), nchunks)
  }
  else {
    spectra_chunks = list(seq_along(in_names))
  }

  invisible(bpmapply(
    rw_f,
    spectra_chunks,
    seq_along(spectra_chunks),
    MoreArgs = list(in_files = in_files, out_files = out_files, out_names = out_names,
                    read_f = read_f, many_to_one=many_to_one,
                    add_ext = add_ext, BPPARAM=BPPARAM),
    BPPARAM = SerialParam(progressbar = verbose)
  ))
}


#' Reads data in writes in mzML format
#'
#' Internal use only. Use [change_format_chunks] instead
#' @importFrom mzR writeMSData
#' @importFrom dplyr bind_rows
#' @importFrom fs is_file
#'
rw_chunk_mzml = function(x, ch, in_files, out_files, out_names,
                         read_f, many_to_one, add_ext, BPPARAM) {
  # Create infiles
  if (!many_to_one) out_files = out_files[x]
  in_files = in_files[x]
  out_names = out_names[x]
  if (add_ext) out_files = paste0(out_files, '.mzML')

  l = bplapply(in_files, read_f, BPPARAM=BPPARAM)

  header_template = initialize_mzml_header()

  header_df = bpmapply(
    generate_header,
    l, out_names, rep(1, length(l)),
    MoreArgs = list(header_template = header_template),
    SIMPLIFY = FALSE, BPPARAM=BPPARAM
  )

  if (!many_to_one) { # save multiple  files in outpath
    invisible(bpmapply(
      function(x, outpath, header) {
        writeMSData(list(x), header = header, file = outpath,
                    backend = 'pwiz', outformat = 'mzml')
      },
      l, out_files, header_df,
      BPPARAM = BPPARAM
    ))
  } else {
    header_df = do.call(bind_rows, header_df)
    header_df$seqNum = 1:length(l)
    writeMSData(l, file = out_files, header = header_df,
                backend = 'pwiz', outformat = 'mzml')
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




#' Reads data in writes in tsv format
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
#' @param idx Scan number of index
#' @param outpath Path to output folder
#' @param header_template Header template to complete with spectra data and id
#'
#' @export
#' @importFrom mzR writeMSData
#'
export_mzml = function(x, outpath, header) {
  writeMSData(list(x), header = header, file = outpath,
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





