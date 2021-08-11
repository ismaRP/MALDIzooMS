#' Title
#'
#' @param platemap
#' @param basepath
#'
#' @return
#' @importFrom dplyr arrange filter mutate group_by
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
collect_triplicates = function(file, basepath, ext, readf, outpath, writef,
                               keep_dupl, keep_incomplete_tripl) {

  platemap = platemap %>%
    arrange(basename, sample_name)

  if (any(duplicated(platemap$sample_name), na.rm = T)) {
    if (keep_dupl == F) {
      warning("Duplicated samples, only the last one will be considered")
    } else {
      duplidx = with(platemap,
        ave(seq_along(sample_name), sample_name, FUN=seq_along))
      dupl_all = with(platemap,
        duplicated(sample_name) | duplicated(sample_name, fromLast = T))
      platemap %>% mutate(
        sample_name = ifelse(dupl_all, paste0(sample_name, "__", duplidx),
                             sample_name)
      )
    }
  }


  # Gather platemap
  platemap = platemap %>%
    gather("spot", "coordinates", c(spot1, spot2, spot3)) %>%
    arrange(folder, subfolder, quire) %>%
    mutate(replicate = substr(spot, 5, 5)) %>%
    mutate(inpath = file.path(path, folder, paste0(subfolder, '.', ext),
                              paste0(basename, '_', spot, '.', ext))) %>%
    mutate(exists = file.exists(inpath))


  if (!keep_incomplete_tripl) {
    platemap = platemap %>% group_by(sample_mame) %>%
      mutate(complete = all(exists)) %>%
      filter(complete)
  } else {
    platemap = platemap %>% filter(exists)
  }
  platemap = platemap %>%
    mutate(outpath = file.path(outpath, paste0(sample_name, "_", replicate)))

  file.copy(platemap$inpath, platemap$outpath)

}

