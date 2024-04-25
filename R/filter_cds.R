#' @title Filter a CDS
#' @description This function provides a pipe-friendly method to filter cds objects.
#' @param cds The CDS to filter.
#' @param cells Optional:  a tibble of cell metadata for the cells you wish to keep.  Use bb_cellmeta(). Default: 'all'
#' @param genes Optional:  a tibble of gene metadata for the genes you wish to keep.  Use bb_rowmeta().  Default: 'all'
#' @return A filtered CDS
#' @import tidyverse
#' @importFrom cli cli_abort
#' @rdname filter_cds
#' @export
filter_cds <- function(cds, cells = "all", genes = "all") {
  if (all(cells == "all")) {
    cells <- bb_cellmeta(cds) %>%
      pull(cell_id)
  } else {
    cells <- cells %>%
      pull(cell_id)
    if (length(cells) == 0) cli::cli_abort("Your filtering removed all cells.  Aborting.")
  }
  if (all(genes == "all")) {
    genes <- bb_rowmeta(cds) %>%
      pull(feature_id)
  }  else {
    genes <- genes %>%
      pull(feature_id)
    if (length(genes) == 0) cli::cli_abort("Your filtering removed all genes.  Aborting.")
  }
  result <- cds[genes, cells]
  return(result)
}
