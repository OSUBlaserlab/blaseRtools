#' @title Make a Heatmap of Aggregated Gene Expressionf from a CDS
#' @description Supply a cds subset and aggregation values and get a heatmap.  Plot the return value using cowplot::plot_grid().  This function wraps several complicated functions from ComplexHeatmap and tries to apply default values for the most common use case:  plotting top markers from a cds with cells grouped arbitrarily (usually by cluster of some type).  This function provides the option to aggregate by genes as well in order to plot gene modules.  For cell and gene aggregation, provide a name (in the form of a string) of the corresponding metadata column and the aggregation will be performed.  There are many aesthetic parameters for this function and even more available for the internal ComplexHeatmap::Heatmap() function.  The best way to adjust parameters not provided here is to scroll the popup box that appears when you type ComplexHeatmap::Heatmap() in RStudio and pass them in via the ellipsis.  The default is to put genes in columns and cells in rows.  You can flip this behavior by setting flip_axis = TRUE.  This will require adjustment of some aesthetic parameters.  Complex annotations (beyond labeling cherry-picked genes) are not currently supported.
#' @param cds_subset The subset of cells and genes you want to plot as a heatmap.  Best approach is to pipe the cds through filter_cds() and into this function.
#' @param cellmeta_col The name of a cell metadata column to aggregate cells by; one of cellmeta_col and rowmeta_col must not be NULL, Default: NULL
#' @param rowmeta_col The name of a row metadata column to aggregate cells by; one of cellmeta_col and rowmeta_col must not be NULL, Default: NULL
#' @param heatmap_highlights A vector of gene names to highlight using anno_mark(), Default: NULL
#' @param three_colors A vector of colors for the main color scale, Default: c("blue4", "ivory", "red3")
#' @param flip_axis Logical; whether to plot genes as rows (TRUE) or columns (FALSE), Default: FALSE
#' @param name Name of the main color scale, Default: NULL
#' @param heatmap_legend_param Graphical parameters for the main heatmap legend, Default: list(title_gp = gpar(fontface = "plain", fontsize = 9), grid_width = unit(0.14,
#'    "in"), labels_gp = gpar(fontsize = 8))
#' @param row_dend_width Row dendrogram width, Default: unit(5, "mm")
#' @param column_dend_height Column dendrogram height, Default: unit(5, "mm")
#' @param column_dend_side Side on which to plot the column dendrogram, Default: 'bottom'
#' @param show_row_names Logical; whether or not to show rownames, Default: T
#' @param row_names_gp Graphical parameters for the row names, Default: gpar(fontsize = 9)
#' @param show_column_names Logical; whether or not to show column names, Default: F
#' @param row_dend_gp Graphical parameters for the row dendrogram, Default: gpar(lwd = 0.5)
#' @param column_dend_gp Graphical parameters for teh column dendrogram, Default: gpar(lwd = 0.5)
#' @param row_title Row title text, Default: NULL
#' @param column_title Column title text, Default: NULL
#' @param padding Padding between gene names on the heatmap highlights, Default: 1.5
#' @param labels_rot Rotation of the heatmap highlight labels, Default: 45
#' @param ... Optional arguments to pass to ComplexHeatmap::Heatmap()
#' @return A complex heatmap in the form of a gtree.
#' @rdname bb_cds_heatmap
#' @import circlize ComplexHeatmap tidyverse monocle3
#' @export
bb_cds_heatmap <- function(cds_subset,
                           cellmeta_col = NULL,
                           rowmeta_col = NULL,
                           heatmap_highlights = NULL,
                           three_colors = c("blue4", "ivory", "red3"),
                           flip_axis = FALSE,
                           name = NULL,
                           heatmap_legend_param = list(
                             title_gp = gpar(fontface = "plain", fontsize = 9),
                             grid_width = unit(0.14, "in"),
                             labels_gp = gpar(fontsize = 8)
                           ),
                           row_dend_width = unit(5, "mm"),
                           column_dend_height = unit(5, "mm"),
                           column_dend_side = "bottom",
                           show_row_names = T,
                           row_names_gp = gpar(fontsize = 9),
                           show_column_names = F,
                           row_dend_gp = gpar(lwd = 0.5),
                           column_dend_gp = gpar(lwd = 0.5),
                           row_title = NULL,
                           column_title = NULL,
                           padding = 1.5,
                           labels_rot = 45,
                           ...) {
  # generate the cell and gene groupings from cds metadata
  if (!is.null(cellmeta_col)) {
    cell_groupings <- bb_cellmeta(cds_subset) |>
      select(cell_id, matches(cellmeta_col))
  } else {
    cell_groupings <- NULL
  }

  if (!is.null(rowmeta_col)) {
    gene_groupings <- bb_rowmeta(cds_subset) |>
      select(feature_id, matches(rowmeta_col))
  } else {
    gene_groupings <- NULL
  }

  # aggregate gene expression

  agg_mat <- aggregate_gene_expression(
    cds = cds_subset,
    cell_group_df = cell_groupings,
    gene_group_df = gene_groupings,
    norm_method = "log",
    pseudocount = 1
  )

  # convert from sparse to regular matrix
  agg_mat <- as.matrix(agg_mat)

  # fix the rownames

  if (is.null(rowmeta_col)) {
    rownames(agg_mat) <-
      left_join(tibble(feature_id = rownames(agg_mat)),
                bb_rowmeta(cds_subset)) |>
      pull(gene_short_name)
  }

  # transpose and then put all of the genes (columns) on the same scale
  agg_mat <- scale(t(agg_mat))

  # make the heatmap color scale
  col_fun <-
    colorRamp2(breaks = c(min(agg_mat),
                          0,
                          max(agg_mat)),
               colors = three_colors)

  # make the annotation object
  if (!is.null(heatmap_highlights)) {
    heatmap_anno_df <-
      map(
        .x = heatmap_highlights,
        .f = function(x) {
          index <- which(colnames(agg_mat) == x)
          return(index)
        }
      ) %>% set_names(heatmap_highlights) %>%
      bind_cols() %>%
      pivot_longer(everything()) %>%
      as.data.frame()

    heatmap_gene_anno <- HeatmapAnnotation(
      foo = anno_mark(
        at = heatmap_anno_df$value,
        labels = heatmap_anno_df$name,
        labels_gp = gpar(fontsize = 8),
        padding = padding,
        labels_rot = labels_rot
      ),
      which = ifelse(flip_axis, "row", "column")
    )

  }
  # make the heatmap finally

  Heatmap_normal <-  function(anno, ...) {
    Heatmap(
      matrix = agg_mat,
      col = col_fun,
      name = name,
      heatmap_legend_param = heatmap_legend_param,
      row_dend_width = row_dend_width,
      column_dend_height = column_dend_height,
      column_dend_side = column_dend_side,
      show_row_names = show_row_names,
      row_names_gp = row_names_gp,
      show_column_names = show_column_names,
      top_annotation = anno,
      row_dend_gp = row_dend_gp,
      column_dend_gp = column_dend_gp,
      row_title = row_title,
      column_title = column_title,
      ...
    )
  }

  Heatmap_flipped <-  function(anno, ...) {
    Heatmap(
      matrix = t(agg_mat),
      col = col_fun,
      name = name,
      heatmap_legend_param = heatmap_legend_param,
      row_dend_width = row_dend_width,
      column_dend_height = column_dend_height,
      column_dend_side = column_dend_side,
      show_row_names = show_row_names,
      row_names_gp = row_names_gp,
      show_column_names = show_column_names,
      right_annotation = anno,
      row_dend_gp = row_dend_gp,
      column_dend_gp = column_dend_gp,
      row_title = row_title,
      column_title = column_title,
      ...
    )
  }

  if (flip_axis) {
    Heatmap_use <- Heatmap_flipped
  } else {
    Heatmap_use <- Heatmap_normal
  }

  if (is.null(heatmap_highlights)) {
    return_heatmap <-
      grid.grabExpr(draw(Heatmap_use(anno = NULL, ...)), wrap = T)
  } else {
    return_heatmap <-
      grid.grabExpr(draw(Heatmap_use(anno = heatmap_gene_anno, ...)), wrap = T)
  }
  return(return_heatmap)

}
