#' Align a CDS object according to a metadata variable 
#' 
#' @param cds A cell data set object to align
#' @param align_by A metadata column to align by  
#' @param n_cores Number of cores to reduce dimensions by.  
#' @return A modified cell data set object with aligned dimensions and new metadata columns holding prealignment umap coordinates.
#' @export
#' @import tidyverse monocle3
bb_align <- function(cds, align_by, n_cores = 8) {
  cds_aligned <- align_cds(cds = cds, alignment_group = align_by)
  prealignment_dims <- plot_cells(cds)[["data"]] %>% as_tibble() %>% select(data_dim_1,data_dim_2)
  colData(cds_aligned)$prealignment_dim1 <- prealignment_dims$data_dim_1
  colData(cds_aligned)$prealignment_dim2 <- prealignment_dims$data_dim_2
  cds_aligned <- reduce_dimension(cds_aligned, cores = n_cores)
  return(cds_aligned)
}
