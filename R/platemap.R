#' Title
#'
#' @param platemap
#' @param basepath
#'
#' @return
#' @importFrom readr read_csv
#' @export
#'
#' @examples
collect_triplicates = function(file, basepath, readf, outpath, writef,
                               keep_dupl) {

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
           write_f = exportMzMl
         }
  )

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

  platemap = filter_nonexsist(platemap)

  a = do.call(
    mapply,
    c(list('FUN'=read_from_plate, MoreArgs=list(path=basepath, outpath=outpath),
           'SIMPLIFY'=F, 'USE.NAMES'=F), as.list(platemap))
  )

}


move_samples = function() {

}

filter_nonexsist = function(platemap, keep_incomplete_tripl){
  platemap = platemap %>%
    mutate(
      path1 = file.path(path, folder, paste0(subfolder, '.txt'),
                        paste0(basename, '_', spot1, '.txt')),
      path2 = file.path(path, folder, paste0(subfolder, '.txt'),
                        paste0(basename, '_', spot2, '.txt')),
      path3 = file.path(path, folder, paste0(subfolder, '.txt'),
                        paste0(basename, '_', spot3, '.txt'))
    ) %>%
    mutate(
      exists1 = file.exists(path1),
      exists2 = file.exists(path2),
      exists3 = file.exists(path3),
    )

  if (keep_incomplete_tripl == T) {
    platemap = platemap %>% filter( exists1 | exists2 | exists3 )
  } else {
    platemap = platemap %>% filter( exists1, exists2, exists3 )
  }

  return(platemap)
}
