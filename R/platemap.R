#' Title
#'
#' @param platemap
#' @param basepath
#' @param ext
#' @param keep_dupl
#' @param keep_incomplete_tripl
#' @param outfolder
#' @param dry_run
#' Don't copy spectra, just return an updated platemap, with possible duplicates
#' included, gathered and with input and output filenames
#'
#'
#' @return
#' @importFrom dplyr arrange filter mutate group_by
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
collect_triplicates = function(platemap, basepath, ext, outfolder,
                               keep_dupl, keep_incomplete_tripl, dry_run=F) {

  platemap = platemap %>%
    arrange(basename, sample_name)

  if (any(duplicated(platemap$sample_name), na.rm = T)) {
    if (keep_dupl == F) {
      warning("Duplicated samples, only the last one will be considered")
    } else {
      warning("Duplicated samples, adding suffix to duplicates to make unique sample names")
      duplidx = with(platemap,
        ave(seq_along(sample_name), sample_name, FUN=seq_along))
      dupl_all = with(platemap,
        duplicated(sample_name) | duplicated(sample_name, fromLast = T))

      platemap = platemap %>%
        mutate(duplicated = dupl_all) %>%
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
    mutate(replicate = substr(spot, 5, 5)) %>%
    mutate(inpath = file.path(basepath, folder,
                              paste0(subfolder, '.', ext),
                              paste0(basename, '_', coordinates, '.', ext))) %>%
    mutate(exists = file.exists(inpath))

  if (!keep_incomplete_tripl) {
    platemap = platemap %>% group_by(sample_name) %>%
      mutate(complete = all(exists)) %>% ungroup() %>%
      filter(complete)
  } else {
    platemap = platemap %>% filter(exists)
  }

  platemap = platemap %>%
    mutate(outpath = file.path(outfolder,
                               paste0(sample_name, "_", replicate, ".", ext)))

  if (!dry_run){
    invisible(file.copy(platemap$inpath, platemap$outpath))
  }

  return(platemap)

}

