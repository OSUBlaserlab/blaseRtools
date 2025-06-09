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

  if (is.character(cells) && length(cells) == 1 && cells == "all") {
    cells <- bb_cellmeta(cds) %>%
      pull(cell_id)
  } else if (is.data.frame(cells) && nrow(cells) > 0) {
    cells <- cells %>%
      pull(cell_id)
  } else if (is.data.frame(cells) && nrow(cells) == 0) {
      cli::cli_abort("Your filtering removed all cells.  Aborting.")

  } else {
    cli::cli_abort("You need to supply a dataframe/tibble for the cells argument.")
  }

  if (is.character(genes) && length(genes) == 1 && genes == "all") {
    genes <- bb_rowmeta(cds) %>%
      pull(feature_id)
  } else if (is.data.frame(genes) && nrow(genes) > 0) {
    genes <- genes %>%
      pull(feature_id)
  } else if (is.data.frame(genes) && nrow(genes) == 0) {
    cli::cli_abort("Your filtering removed all genes.  Aborting.")

  } else {
    cli::cli_abort("You need to supply a dataframe/tibble for the genes argument.")
  }

  # result <- list(genes, cells)
  result <- cds[genes, cells]
  return(result)
}
