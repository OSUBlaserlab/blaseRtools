#' @title A function to generate clusters from scRNA-seq data
#' @description Based on Monocle3's Partitions, Leiden, and Louvain clustering methods.  Implemented mostly with default values.  Seurat objects will be converted to cell_data_set objects for the clustering.  The function produces a list of top markers for each cluster type and returns these assignments to the original object as new cell metadata columnts.
#' @param obj A Seurat or cell_data_set object
#' @param n_top_markers Number of top markers to identify per cell group, Default: 50
#' @param outfile Name of a csv file to hold the top marker results.  If null, will place "top_markers.csv" in the working directory, Default: NULL
#' @param n_cores Number of processor cores to use, Default: 8
#' @param cds Provided for backwards compatibility for existing code.  If a value is supplied it will be transferred to obj and a warning message will be emitted, Default: NULL
#' @return A modified Seurat or cell_data_set object
#' @rdname bb_triplecluster
#' @export
#' @importFrom monocle3 cluster_cells clusters partitions top_markers
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr mutate relocate bind_rows arrange select
#' @importFrom readr write_csv
bb_triplecluster <- function(obj,
                             n_top_markers = 50,
                             outfile = NULL,
                             n_cores = 8,
                             cds = NULL) {
  cds_warn(cds)
  obj_stop(obj)

  if ("Seurat" %in% class(obj)) {
    cds <- as.cell_data_set(obj)
  } else {
    cds <- obj
  }

  # cluster the cells----------------------------------------------------------
  cds <-
    monocle3::cluster_cells(cds, verbose = TRUE, cluster_method = "leiden")
  cds_louvain <-
    monocle3::cluster_cells(cds, verbose = TRUE, cluster_method = "louvain")

  # make explicit columns----------------------------------------------
  SummarizedExperiment::colData(cds)$leiden <-
    monocle3::clusters(cds)
  SummarizedExperiment::colData(cds)$partition <-
    monocle3::partitions(cds)
  SummarizedExperiment::colData(cds)$louvain <-
    monocle3::clusters(cds_louvain)

  #calculate top markers------------------------------------------------------------
  leiden_top_markers <-
    monocle3::top_markers(
      cds,
      group_cells_by = "leiden",
      reference_cells = 1000,
      cores = n_cores,
      genes_to_test_per_group = n_top_markers
    ) %>%
    dplyr::mutate(cluster_method = "leiden") %>%
    dplyr::relocate(cluster_method, .before = cell_group) %>%
    dplyr::mutate(cell_group = paste0("leiden ", cell_group))

  partition_top_markers <-
    monocle3::top_markers(
      cds,
      group_cells_by = "partition",
      reference_cells = 1000,
      cores = n_cores,
      genes_to_test_per_group = n_top_markers
    ) %>%
    dplyr::mutate(cluster_method = "partition") %>%
    dplyr::relocate(cluster_method, .before = cell_group) %>%
    dplyr::mutate(cell_group = paste0("partition ", cell_group))

  louvain_top_markers <-
    monocle3::top_markers(
      cds,
      group_cells_by = "louvain",
      reference_cells = 1000,
      cores = n_cores,
      genes_to_test_per_group = n_top_markers
    ) %>%
    dplyr::mutate(cluster_method = "louvain") %>%
    dplyr::relocate(cluster_method, .before = cell_group) %>%
    dplyr::mutate(cell_group = paste0("louvain ", cell_group))

  # output and return the results-------------------------------------------------
  all_top_markers <- dplyr::bind_rows(list(
    leiden_top_markers,
    partition_top_markers,
    louvain_top_markers
  )) %>%
    dplyr::mutate(cluster_method = factor(cluster_method, levels = c("partition", "leiden", "louvain"))) %>%
    dplyr::arrange(cluster_method)
  if (!is.null(outfile)) {
    readr::write_csv(x = all_top_markers, file = outfile)
  } else {
    readr::write_csv(x = all_top_markers, file = "top_markers.csv")
  }

  obj <-
    bb_tbl_to_coldata(obj,
                      min_tbl = bb_cellmeta(cds) |> dplyr::select(cell_id, louvain, leiden, partition))

  return(obj)


}
