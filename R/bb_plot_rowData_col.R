#' A helper function to generate a data frame in the proper form for aggregate expression plotting with bb_gene_umap.
#'
#' Use as argument for the "gene_or_genes" parameter for bb_gene_umap.
#'
#' @param cds CDS from which to extract the gene metadata.  Should be the same cds as the enclosing function.
#' @param rowData_col Gene metadata column to aggregate by.
#' @param filter_in Subset of values to focus on.  Each will become a facet in the final plot. Default is to keep everything except NA values.
#' @param filter_out Option to filter out any unwanted values. Default is to not filter out anything.
#' @return A data frame in the format needed to pass into bb_gene_umap.
#' @export
#' @import tidyverse monocle3
#' @examples
bb_plot_rowData_col <- function(cds,
                                rowData_col,
                                filter_in = NULL,
                                filter_out = NULL) {
  rowData(cds) %>%
    as_tibble() %>%
    filter(str_detect(!!sym(rowData_col),ifelse(is.null(filter_in),".*",filter_in))) %>%
    filter(str_detect(!!sym(rowData_col),ifelse(is.null(filter_out),"things you don't want",filter_out), negate = TRUE)) %>%
    select(id, rowData_col)
  
  
}
