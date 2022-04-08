#' Make a plot of gene expression in UMAP form
#'
#' @param cds A cell data set object
#' @param gene_or_genes A character vector of genes to plot
#' @param cell_size Size of cells to plot
#' @param alpha Alpha value for overlaid points
#' @param ncol Number of columns in facetted plot
#' @param plot_title Main title of the plot
#' @param color_legend_title Title for the color legend
#' @return A ggplot
#' @export
#' @import tidyverse monocle3
bb_gene_umap <-
  function (cds,
            gene_or_genes,
            cell_size = 1,
            alpha = 1 ,
            ncol = NULL,
            plot_title = NULL,
            color_legend_title = NULL
            ) {

    data <- plot_cells(cds = cds, genes = gene_or_genes)[["data"]]
    # data$feature_label <-
    #   factor(data$feature_label, levels = gene_or_genes)
    background_data <- data %>% filter(is.na(value))
    foreground_data <- data %>% filter(!is.na(value))
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
        aes(
          x = data_dim_1,
          y = data_dim_2,
          color = log10(value)
        ),
        shape = 16,
        size = cell_size,
        alpha = alpha
      ) +
      scale_color_viridis_c() +
      labs(
        x = "UMAP 1",
        y = "UMAP 2",
        color = color_legend_title,
        title = plot_title
      ) +
      facet_wrap(facets = vars(feature_label), ncol = ncol) +
      theme(strip.background = element_blank())+
      theme(plot.title = element_text(hjust = 0.5))
    return(p)
  }
