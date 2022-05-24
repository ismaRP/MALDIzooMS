#' Arrange samples from different platemaps into a common folder
#'
#' @description
#' ZooMS MALDI spectra from the same plate are usually stored in a separate
#' folder with the plate coordinate as file name.
#' Using a table that maps plateID and coordinate to sample name and replicate,
#' this function copies all the spectra into the same folder with appropriate
#' sample name and replicate number.
#'
#' @param platemap
#' Data frame mapping plate and coordinate to sample name and replicate.
#' It must contain the following columns:
#' \itemize{
#'   \item folder: main folder of group of plates.
#'   \item subfolder: subfolder from a specific plate where spetra is stored.
#'   \item basename: prefix of samples in the folder that is preceeded by coordinates.
#'   \item sample_name: (preferably) unique sample name.
#'   \item spot1, spot2, spot3: coordinates of each of the replicates.
#'   \item date: the date each plate was spotted or analyzed.
#'   There can be more than 3 replicates, but it always need to have this name: "spot#"
#' }
#' @param basepath
#' Path to plates
#' @param ext
#' File extension. E.g. mzML, csv, txt, tab
#' @param keep_dupl
#' Logical. Should duplicate samples be kept? If so, an suffix is added.
#' The suffix has this form: "__#".
#' @param keep_incomplete_tripl
#' Logical. When samples come in triplicates, whether incomplete triplicates
#' should be removed
#' @param outfolder
#' Destination path for renamed spectra
#' @param dry_run
#' Don't copy spectra, just return an updated platemap, with possible duplicates
#' included, gathered and with input and output filenames
#'
#' @return
#' A data frame with the modified platemap. Each replicate is now 1 row. Duplicated
#' samples have the "__#" suffix. Duplicated samples are flagged: \code{dupl} will flag
#' all duplicates but the last one analyzed by \code{date} in the platemap. \code{dupl_all}
#' will flag all duplicated samples.
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


collect_triplicates = function(platemap, basepath="", ext=NULL, outfolder=NULL,
                               keep_dupl=T, keep_incomplete_tripl=T, dry_run=F) {

  if (any(is.null(basepath), is.null(ext), is.null(outfolder)) &
      !dry_run) {
    stop("Specifiy path for platemaps, extension and output folder")
  }

  platemap = platemap %>%
    arrange(date, sample_name)

  if (any(duplicated(platemap$sample_name), na.rm = T)) {
    if (keep_dupl == F) {
      warning("Duplicated samples, only the last one will be considered")
    } else {
      warning("Duplicated samples, adding suffix to duplicates to make unique sample names")
      duplidx = with(platemap,
        ave(seq_along(sample_name), sample_name, FUN=seq_along))
      dupl_all = with(platemap,
        duplicated(sample_name) | duplicated(sample_name, fromLast = T))
      dupl = with(platemap,
        duplicated(sample_name, fromLast = T))
      platemap = platemap %>%
        mutate(dupl_all = dupl_all, dupl = dupl) %>%
        mutate(
          sample_name = ifelse(dupl_all,
                               paste0(sample_name, "__", duplidx),
                               sample_name)
      )
    }
  }

  # Gather platemap
  platemap = platemap %>%
    gather("spot", "coordinates", c(spot1, spot2, spot3)) %>%
    arrange(folder, subfolder, sample_name) %>%
    mutate(replicate = as.numeric(substr(spot, 5, 5))) %>%
    mutate(inpath = file.path(basepath, folder,
                              paste0(subfolder, '.', ext),
                              paste0(basename, '_', coordinates, '.', ext))) %>%
    mutate(exists = file.exists(inpath))

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

