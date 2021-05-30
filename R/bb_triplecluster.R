#' A function to generate clusters from scRNA-seq data
#'
#' Based on Monocle3's Partitions, Leiden, and Louvain clustering methods.  Implemented mostly with default values.
#'
#' @param cds_string A cell data set object
#' @param n_top_markers Number of top markers to identify per cell group
#' @param outfile Name of a csv file to hold the top marker results.
#' @param n_cores Number of processor cores to use
#' @return A modified cell data set object with leiden, louvain, and partition metadata columns.  Outputs csv files with top markers into specified directory, or current working directory if NULL.
#' @export
#' @import tidyverse monocle3
#' @examples
bb_triplecluster <- function(cds,
                             n_top_markers = 50,
                             outfile,
                             n_cores = 8) {
  
  
  # cluster the cells----------------------------------------------------------
  cds <-
    cluster_cells(cds, verbose = TRUE, cluster_method = "leiden")
  cds_louvain <-
    cluster_cells(cds, verbose = TRUE, cluster_method = "louvain")
  
  # make explicit columns----------------------------------------------
  colData(cds)$leiden <- clusters(cds)
  colData(cds)$partition <- partitions(cds)
  colData(cds)$louvain <- clusters(cds_louvain)
  
  #calculate top markers------------------------------------------------------------
  leiden_top_markers <-
    top_markers(
      cds,
      group_cells_by = "leiden",
      reference_cells = 1000,
      cores = n_cores,
      genes_to_test_per_group = n_top_markers
    ) %>%
    mutate(cluster_method = "leiden") %>%
    relocate(cluster_method, .before = cell_group) %>%
    mutate(cell_group = paste0("leiden ", cell_group))
  
  partition_top_markers <-
    top_markers(
      cds,
      group_cells_by = "partition",
      reference_cells = 1000,
      cores = n_cores,
      genes_to_test_per_group = n_top_markers
    ) %>%
    mutate(cluster_method = "partition") %>%
    relocate(cluster_method, .before = cell_group) %>%
    mutate(cell_group = paste0("partition ", cell_group))
  
  louvain_top_markers <-
    top_markers(
      cds,
      group_cells_by = "louvain",
      reference_cells = 1000,
      cores = n_cores,
      genes_to_test_per_group = n_top_markers
    ) %>%
    mutate(cluster_method = "louvain") %>%
    relocate(cluster_method, .before = cell_group) %>%
    mutate(cell_group = paste0("louvain ", cell_group))
  
  # output and return the results-------------------------------------------------
  all_top_markers <- bind_rows(list(
    leiden_top_markers,
    partition_top_markers,
    louvain_top_markers
  )) %>%
    mutate(cluster_method = factor(cluster_method, levels = c("partition", "leiden", "louvain"))) %>%
    arrange(cluster_method) %>%
    write_csv(outfile)
  return(cds)
  
  
}
