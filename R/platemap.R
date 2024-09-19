#' Arrange samples from different platemaps into a common folder
#'
#' @description
#' ZooMS MALDI spectra from the same plate are usually stored in a separate
#' folder with the plate coordinate as file name.
#' Using a table that maps plate ID and coordinates to sample name and replicate,
#' this function copies all the spectra into the same folder with appropriate
#' sample name and replicate number.
#'
#' @param platemap
#' Data frame mapping plate and coordinate to sample name and replicate.
#' It can contain the following columns:
#' \itemize{
#'   \item folder: main folder of group of plates.
#'   \item subfolder: subfolder from a specific plate where spetra is stored.
#'   \item basename: prefix of samples in the folder that is preceeded by coordinates.
#'   \item sample_name: (preferably) unique sample name.
#'   \item spot1, spot2, spot3: (for wide format) coordinates of each of the replicates.
#'         There can be more than 3 replicates, but the column name always need
#'         to have this format: "spot#"
#'   \item spot: (for long format), coordinates of the replicate.
#'   \item date (optionl): the date each plate was spotted or analyzed.
#' }
#'
#' @param format
#' Whether the platemap is in "wide" or "long" format.
#' In the wide format, each row is one sample and the coordinates of the
#' replicates are in different columns as spot1, spot2, ...
#' In the long format, each row is an individual replicate. In this case, the
#' platemap must contain a column called "replicate" with the replicate number
#' and the column with the coordinates should be called "spot".
#' Long format assumes no duplicates.
#' @param basepath
#' Path to plates
#' @param ext
#' File extension. E.g. mzML, csv, txt, tab
#' @param keep_dupl
#' Logical. Should duplicate samples/replicates be kept? If so, a suffix is added.
#' The suffix has this form: "__#".
#' @param keep_incomplete_tripl
#' Logical. When samples come in triplicates, whether incomplete triplicates
#' should be skipped.
#' @param outfolder
#' Destination path for renamed spectra.
#' @param dry_run
#' Don't copy spectra, just return an updated platemap, with possible duplicates
#' included, in long format and with input and output filenames
#'
#' @return
#' A data.frame with the modified platemap in long format, so each replicate is
#' now in 1 row. Duplicated samples, if kept, have the "__#" suffix.
#' Duplicated samples are flagged with \code{dupl}. All duplicates but the last
#' one, by \code{date}, in the platemap are flagged.
#' Column \code{dupl_all} will flag all duplicated samples.
#' Path for each sample to the original plate location is included in a column.
#' @importFrom dplyr arrange filter mutate group_by
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' platemap = read_csv('path/to/platemap.csv')
#' platemap_1 = collect_triplicates(
#'   platemap,
#'   basepath = 'path/to/platefolder', ext = "txt",
#'   outfolder = 'path/to/samples',
#'   keep_dupl = T, keep_incomplete_tripl = T, dry_run=F)
#' }
collect_triplicates = function(platemap, format=c('wide', 'long'), basepath="", ext=NULL, outfolder=NULL,
                               keep_dupl=T, keep_incomplete_tripl=T, dry_run=F) {

  format = match.arg(format)

  if (any(is.null(basepath), is.null(ext), is.null(outfolder)) &
      !dry_run) {
    stop("Specifiy path for platemaps, extension and output folder")
  }

  if (format == 'wide'){
    # Gather platemap
    platemap = platemap %>%
      gather("replicate", "spot", c(spot1, spot2, spot3)) %>%
      arrange(folder, subfolder, sample_name) %>%
      mutate(replicate = as.numeric(substr(replicate, 5, 5))) %>%
      mutate(inpath = file.path(basepath, folder,
                                paste0(subfolder, '.', ext),
                                paste0(basename, '_', spot, '.', ext))) %>%
      mutate(exists = file.exists(inpath))
    if (any(duplicated(paste0(platemap$sample_name, '_', platemap$replicate)), na.rm = T)) {
      if (keep_dupl == F) {
        warning("Duplicated samples, only the last one will be considered")
      } else {
        warning("Duplicated samples, adding suffix to duplicates to make unique sample names")
        duplidx = with(platemap,
                       ave(seq_along(sample_name), paste0(sample_name, '_', replicate), FUN=seq_along))
        dupl_all = with(platemap,
                        duplicated(paste0(sample_name, '_', replicate)) | duplicated(paste0(sample_name, '_', replicate), fromLast = T))
        dupl = with(platemap,
                    duplicated(paste0(sample_name, '_', replicate), fromLast = T))
        platemap = platemap %>%
          mutate(dupl_all = dupl_all, dupl = dupl) %>%
          mutate(
            sample_name = ifelse(dupl_all,
                                 paste0(sample_name, "__", duplidx),
                                 sample_name)
          )
      }
    }
  } else { # format == 'long'
    platemap = platemap %>%
      mutate(inpath = file.path(basepath, folder,
                                paste0(subfolder, '.', ext),
                                paste0(basename, '_', spot, '.', ext))) %>%
      mutate(exists = file.exists(inpath))
  }

  if (!dry_run){
    if (!keep_incomplete_tripl) {
      platemap_files = platemap %>% group_by(sample_name) %>%
        mutate(complete = all(exists)) %>% ungroup() %>%
        filter(complete)
    } else {
      platemap_files = platemap %>% filter(exists)
    }

    platemap_files = platemap_files %>%
      mutate(outpath = file.path(outfolder,
                                 paste0(sample_name, "_", replicate, ".", ext)))

    invisible(file.copy(platemap_files$inpath,
                        platemap_files$outpath))
  }

  return(platemap)

}


#' Extract spectra name from a full path with extensions
#'
#' Spectra name has the format "samplename_replicate"
#'
#' @param files Paths to the individual spectra files including the file name and extension
#'
#' @return A vector with sample names
#' @export
#' @importFrom fs path_file
#'
get_spectra_name = function(files){
  files = strsplit(fs::path_file(files), '\\.')
  spectra_name = sapply(files, '[[', 1)

  return(spectra_name)
}



#' Separate spectra names into sample name and replicate.
#'
#' Spectra names have the format sample_replicate.
#'
#' @param spectra_names Spectra names
#' @param sep Separator of sample name and replicate. Default is '_'
#'
#' @return \code{separate_sample_replicate} returns a new data.frame with columns
#' sample and replicate
#' @export
#'
separate_sample_replicate = function(spectra_names, sep='_'){
  parts = strsplit(spectra_names, '_')
  samples = sapply(parts, function(x) paste0(x[-length(x)], collapse='_'))
  replicates = sapply(parts, function(x) x[[length(x)]])

  return(data.frame(sample=samples, replicate=replicates))
  # ((?:[a-zA-Z1-9]+_*[a-zA-Z1-9]*)*?)(_)(\d)$

}


#' @rdname separate_sample_replicate
#' @param sp data.frame with data, in which the spectra name is in 'spectra_name'
#'           column
#' @param sample_sep Separator found within a sample name
#' @param repl_seq Separator between sample name and replicate number
#'
#' @return \code{separate_sample_replicate_df} returns the data.frame in \code{sp}
#' with sample and replicate columns
#' @export
#'
separate_sample_replicate_df = function(sp, sample_sep='_', repl_sep='_'){
  sample_sep = paste0(sample_sep, collapse='')
  separated = separate_wider_regex(
    sp, spectra_name,
    patterns=c('sample'=sprintf('(?:[a-zA-Z1-9]+(?:[%s][a-zA-Z1-9])*)*?', sample_sep),
               repl_sep,
               'replicate'='\\d$'),
    cols_remove=FALSE)
  return(separated)
}


#' Clean spectra metadata
#'
#' Remove spectra in metadata that don't exist in data
#'
#' @param metadata A data.frame with a replicate spectra per row. It must contain
#' at least columns sample_name and replicate
#' @param folder Folder where spectra data is stored
#'
#' @return data.frame metadata with extra columns \code{file},
#' \code{spectra_name} and \code{n_replicates}
#' @export
#' @importFrom dplyr group_by mutate ungroup
#'
clean_metadata = function(metadata, folder) {

  spectra_files = dir(folder)
  spectra_names = sapply(strsplit(spectra_files, '\\.'), '[[', 1)

  data = data.frame(
    spectra_name = spectra_names,
    file = spectra_files
  )

  metadata$spectra_name = paste0(
    metadata$sample_name, '_', metadata$replicate)

  idx = match(metadata$spectra_name, data$spectra_name)

  metadata$file = data[idx, 'file']
  metadata = metadata[!is.na(metadata$file),]

  metadata = metadata %>% group_by(sample_name) %>%
    mutate(n_replicates = n()) %>% ungroup()

  return(metadata)

}




