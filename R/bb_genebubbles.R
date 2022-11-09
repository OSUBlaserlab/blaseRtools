#' @title Create a Gene Bubble/Dot Plot
#' @description This is a very data-dense plot and is the recommended way for showing expression of single markers/genes by cell group.  By default, this function will return an unfaceted ggplot with cell groups on the X axis and genes on the Y axis with dot size representing proportion of cells in the cell group expressing a gene and color scale representing per-cell expression.
#'
#' But it also may be of interest to add aesthetic variables such as facets or additional color scales.  There are two ways this function will facilitate that.  First, you can supply a vector of cell groups to the cell_grouping argument and a the cells will be grouped by the composite value of these factors.  Usually if you are doing this, you also will want to have access to the components of this composite variable to facet by.  So you can supply "data" to the return_value argument to get a tibble.  From there you can modify as necessary and generate a ggplot assigning aesthetics and scales as desired and using geom_point.
#' @param obj A Seurat or cell_data_set object.
#' @param genes Gene or genes to plot.
#' @param cell_grouping Cell metadata column to group cells by.  Supply more than one in a vector to generate a composite variable.
#' @param experiment_type Experiment data to plot.  Usually will be either "Gene Expression" or "Antibody Capture", Default: 'Gene Expression'
#' @param scale_expr Whether to scale expression by gene, Default: TRUE
#' @param gene_ordering By default, genes will be ordered by a clustering algorithm.  Supply "as_supplied" to plot the genes in the order supplied to the "genes" argument , Default: c("bicluster", "as_supplied")
#' @param group_ordering By default, cell groups will be ordered by a clustering algorithm.  Supply "as_supplied" to plot the cell groups in the order supplied to "cell_grouping", Default: c("bicluster", "as_supplied")
#' @param return_value Whether to return a plot or data in tibble form, Default: c("plot", "data")
#' @return A ggplot or a tibble
#' @rdname bb_genebubbles
#' @export
#' @importFrom tidyr unite pivot_longer
#' @importFrom tidyselect matches
#' @importFrom dplyr select left_join mutate pull
#' @importFrom tibble as_tibble tibble
#' @importFrom cli cli_alert_warning
#' @importFrom ggplot2 ggplot geom_point scale_size_area scale_color_viridis_c labs
bb_genebubbles <- function(obj,
                           genes,
                           cell_grouping,
                           experiment_type = "Gene Expression",
                           scale_expr = TRUE,
                           gene_ordering = c("bicluster", "as_supplied"),
                           group_ordering = c("bicluster", "as_supplied"),
                           return_value = c("plot", "data")) {
  blaseRtools:::obj_stop(obj)
  gene_ordering <- match.arg(gene_ordering)
  group_ordering <- match.arg(group_ordering)
  return_value <- match.arg(return_value)
  e_t <- experiment_type

  if (length(dim(genes)) > 1)
    stop("This visualization cannot be used for aggregated genes.")
  ids <- blaseRtools:::get_gene_ids(obj, genes, et = e_t)

  # if this cell grouping is not present in the obj, then create it
  cg <- cell_grouping
  min_tbl <- bb_cellmeta(obj) |>
    tidyr::unite(col = "cell_grouping",
                 tidyselect::matches(paste("^", cell_grouping, "$", sep = "")),
                 remove = FALSE) |>
    dplyr::select(cell_id, cell_grouping)
  obj <- bb_tbl_to_coldata(obj, min_tbl)

  mat <- bb_aggregate(
    obj,
    cell_group_df = bb_cellmeta(obj) |>
      dplyr::select(cell_id, cell_grouping),
    experiment_type = e_t
  )
  mat <- mat[rownames(mat) %in% ids,]
  if (class(mat) == "numeric") {
    mat <- as.matrix(mat)
    colnames(mat) <- ids
  } else {
    mat <- t(as.matrix(mat))
  }
  if (scale_expr)
    mat <- scale(mat)
  bin_mat <- bb_aggregate(
    obj,
    norm_method = "binary",
    cell_group_df = bb_cellmeta(obj) |>
      dplyr::select(cell_id, cell_grouping),
    experiment_type = e_t
  )
  bin_mat <- bin_mat[rownames(bin_mat) %in% ids,]
  if (class(bin_mat) == "numeric") {
    bin_mat <- as.matrix(bin_mat)
    colnames(bin_mat) <- ids
  } else {
    bin_mat <- t(as.matrix(bin_mat))
  }
  # put in dataframe
  dat <- tibble::as_tibble(mat, rownames = "group") |>
    tidyr::pivot_longer(-group, names_to = "feature_id", values_to = "expression")
  bin_dat <- tibble::as_tibble(bin_mat, rownames = "group") |>
    tidyr::pivot_longer(-group, names_to = "feature_id", values_to = "proportion")
  plot_data <-
    dplyr::left_join(dat, bin_dat, by = c("group", "feature_id"))

  # add on the gene short names
  rm <- bb_rowmeta(obj, experiment_type = e_t)
  plot_data <-
    dplyr::left_join(plot_data,
                     rm,
                     by = "feature_id")

  # reorder the genes and groups
  if (gene_ordering == "as_supplied") {
    plot_data <-
      dplyr::mutate(plot_data,
                    gene_short_name = factor(gene_short_name, levels = genes))
  } else if (length(genes) == 1) {
    cli::cli_alert_warning("Groups are not clustered when a single gene is supplied.")
  } else {
    gene_order <-
      tibble::tibble(feature_id = levels(bicluster_bubbles(mat)$genes)) |>
      dplyr::left_join(bb_rowmeta(obj, experiment_type = e_t), by = "feature_id") |>
      dplyr::pull(gene_short_name) |>
      unique()
    plot_data <-
      dplyr::mutate(plot_data,
                    gene_short_name = factor(gene_short_name, levels = gene_order))
  }

  if (group_ordering == "bicluster" && length(genes) > 1) {
    plot_data <-
      dplyr::mutate(plot_data, group = factor(group, levels = levels(bicluster_bubbles(mat)$groups)))
  }
  # now join the cell metadata back on in case you want it later
  # return(list(plot_data, bb_cellmeta(obj)))
  group_meta <- bb_cellmeta(obj) |>
    select(matches(paste0("^", c(cg, "cell_grouping"), "$"))) |>
    distinct()
  plot_data <- left_join(plot_data, group_meta, by = c("group" = "cell_grouping"))

  plot <-
    ggplot2::ggplot(
      plot_data,
      mapping = aes(
        x = group,
        y = gene_short_name,
        color = expression,
        size = proportion
      )
    ) +
    ggplot2::geom_point(pch = 16) +
    ggplot2::scale_size_area() +
    ggplot2::scale_color_viridis_c(option = ifelse(experiment_type == "Antibody Capture", "A", "D")) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      color = ifelse(experiment_type == "Antibody Capture", "Binding", "Expression"),
      size = "Proportion"
    )

  if (return_value == "plot") {
    return(plot)
  } else {
    return(plot_data)
  }



}



#' @importFrom stats as.dist cor
#' @importFrom pheatmap pheatmap
bicluster_bubbles <- function(x) {
  row_dist <- stats::as.dist((1 - stats::cor(t(x))) / 2)
  row_dist[is.na(row_dist)] <- 1
  col_dist <- stats::as.dist((1 - stats::cor(x)) / 2)
  col_dist[is.na(col_dist)] <- 1
  ph <-
    pheatmap::pheatmap(
      x,
      useRaster = T,
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      show_rownames = F,
      show_colnames = F,
      clustering_distance_cols = col_dist,
      clustering_distance_rows = row_dist,
      clustering_method = "ward.D2",
      silent = TRUE,
      filename = NA
    )
  genes <-
    factor(colnames(x), levels = colnames(x)[ph$tree_col$order])
  groups <-
    factor(row.names(x), levels = row.names(x)[ph$tree_row$order])
  return(list(groups = groups, genes = genes))
}
