#' @title Plot a UMAP Showing Cite-Seq Antibody Binding
#' @description Requires a cds with an alt experiment established.  Use bb_split_citeseq to generate this and to normalize binding data using the CLR method.  Returns a ggplot.
#' @param cds The cds with an "Antibody Capture" alt experiment to plot.
#' @param antibody The name of the antibody to plot. Equivalent to gene_short_name.  Accepts a character vector.
#' @param cell_size Size of points to plot, Default: 1
#' @param alpha Alpha for the plotted points, Default: 1
#' @param alt_dim_x Alternate/reference dimensions to plot by.
#' @param alt_dim_y Alternate/reference dimensions to plot by.
#' @param plot_title Optional title for the plot, Default: NULL
#' @param color_legend_title Optional title for the color scale., Default: NULL
#' @param order Whether or not to order cells by gene expression.  When ordered, non-expressing cells are plotted first, i.e. on the bottom. Default: TRUE
#' @param rescale Optional redefinition of the color scale, Default: NULL
#' @param ncol If specified, the number of columns for facet_wrap, Default: NULL
#' @return a ggplot
#' @seealso
#'  \code{\link[SingleCellExperiment]{reducedDims}}
#' @rdname bb_cite_umap
#' @export
#' @importFrom SingleCellExperiment swapAltExp reducedDims
#' @importFrom SummarizedExperiment assays
#' @importFrom cli cli_abort
bb_cite_umap <-
  function(cds,
           antibody,
           cell_size = 1,
           alpha = 1,
           alt_dim_x = NULL,
           alt_dim_y = NULL,
           plot_title = NULL,
           color_legend_title = NULL,
           order = TRUE,
           rescale = NULL,
           ncol = NULL) {
    cds_alt <- as(SingleCellExperiment::swapAltExp(cds, name = "Antibody Capture"), Class = "cell_data_set")
    data <- SummarizedExperiment::assays(cds_alt)$CLR_counts
    data <- t(data)
    data_tbl <- as_tibble(data) |>
      mutate(cell_id = rownames(data))

    # find the column we want to plot and rename it
    antibody_id <- bb_rowmeta(cds_alt) |>
      filter(gene_short_name %in% antibody) |>
      pull(id)
    if (length(antibody_id) == 0)
      cli::cli_abort("The requested antibody is not available in the data object.")
    data_tbl <- data_tbl |>
      select(cell_id, matches(antibody_id)) |>
      pivot_longer(cols = -cell_id) |>
      rename(antibody_id = name, binding = value) |>
      left_join(
        bb_rowmeta(cds_alt) |>
          filter(gene_short_name %in% antibody) |>
          select(antibody_id = feature_id, gene_short_name),
        by = "antibody_id"
      ) |>
      select(cell_id, antibody = gene_short_name, binding)

    dims <- SingleCellExperiment::reducedDims(cds)$UMAP
    colnames(dims) <- c("data_dim_1", "data_dim_2")

    dims_tbl <- as_tibble(dims) |>
      mutate(cell_id = rownames(dims))
    plot_tbl <- full_join(dims_tbl, data_tbl, by = "cell_id")

    # make the plot
    # plot <- ggplot(plot_tbl, mapping = aes(x = data_dim_1, y = data_dim_2, ))
    background_data <- plot_tbl |> filter(is.na(binding))
    foreground_data <- plot_tbl |> filter(!is.na(binding))

    # optionally order the cells to un-bury rare expressing cells
    if (order)
      foreground_data <- dplyr::arrange(foreground_data, !is.na(binding), binding)

    p <- ggplot() +
      geom_point(
        data = background_data,
        aes(x = data_dim_1, y = data_dim_2),
        color = "grey80",
        shape = 1,
        size = cell_size,
        stroke = 0.25
      ) +
      geom_point(
        data = foreground_data,
        aes(x = data_dim_1,
            y = data_dim_2,
            color = binding),
        shape = 16,
        size = cell_size,
        alpha = alpha
      )

    if (!is.null(rescale)) {
      p <- p +
        scale_color_viridis_c(option = "A", limits = rescale, na.value = "grey80")
    } else {
      p <- p +
        scale_color_viridis_c(option = "A")
    }
    p <- p +
      labs(
        x = ifelse(is.null(alt_dim_x), "UMAP 1", alt_dim_x),
        y = ifelse(is.null(alt_dim_y), "UMAP 2", alt_dim_y),
        color = color_legend_title,
        title = plot_title
      ) +
      facet_wrap(facets = vars(antibody), ncol = ncol) +
      theme(strip.background = element_blank()) +
      theme(plot.title = element_text(hjust = 0.5))
    return(p)
  }



#' @title Split Antibody Capture Data into Alt Experiment
#' @description If you have cite-seq data together with gene expression data, this function will move the cite seq data to a new separate experiment.  It will use Seurat to normalize these data using the CLR method and store them in a new assay.
#' @param cds the cell data set to split
#' @return a new CDS
#' @seealso
#'  \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @rdname bb_split_citeseq
#' @export
#' @importFrom SingleCellExperiment splitAltExps swapAltExp
#' @importFrom SummarizedExperiment assay
bb_split_citeseq <- function(cds) {
  check <- "Antibody Capture" %in% bb_rowmeta(cds)$data_type
  stopifnot("Your cds does not have any Antibody data." = check)
  # create the alternative experiment, making "Antibody Capture" the new reference
  cds <- SingleCellExperiment::splitAltExps(cds, rowData(cds)$data_type, ref = "Antibody Capture")

  # normalize using Seurat functions
  CLR_counts <- get_seurat_clr(cds)

  # set this as a new assay
  SummarizedExperiment::assay(cds, "CLR_counts", withDimnames = FALSE) <- CLR_counts

  # swap experiments back
  cds <- SingleCellExperiment::swapAltExp(cds, name = "Gene Expression")

  # convert from SingleCellExperiment class back to cell_data_set
  cds <- as(object = cds, Class = "cell_data_set")

  return(cds)

}



#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom SingleCellExperiment counts
#' @importFrom Seurat NormalizeData as.sparse
get_seurat_clr <- function(cds) {
  data <-
    SeuratObject::CreateSeuratObject(counts = SingleCellExperiment::counts(cds), assay = "ADT")
  data <-
    Seurat::NormalizeData(
      data,
      normalization.method = "CLR",
      margin = 2,
      assay = "ADT"
    )
  data <- data@assays$ADT@data
  data <- Seurat::as.sparse(data)
  return(data)
}
