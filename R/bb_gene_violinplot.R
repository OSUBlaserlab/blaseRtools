#' Make a plot of gene expression in UMAP form
#'
#' @param cds A cell data set object
#' @param variable Stratification variable for x-axis
#' @param genes_to_plot Either a character vector of gene short names or a tbl/df where the first column is gene short name and the second is the gene grouping.
#' @param pseudocount Value to add to zero-cells
#' @param include_jitter Include jitter points
#' @param ytitle Title for y axis
#' @param plot_title Main title for the plot
#' @param rows Number of rows for facetting
#' @param show_x_label Option to show x label
#' @param legend_pos Position for label
#' @param comparison_list Optional list of comparisons for ggpubr
#' @param palette Color palette to use.  Viridis is default.
#' @param violin_alpha Alpha value for violin plot
#' @param jitter_alpha Alpha value for jitter plot
#' @param jitter_color Color for the jitter plot.  Defaults to black and ignored if jitter_match == TRUE
#' @param jitter_fill Fill for the jitter plot
#' @param jitter_size Size of the jitter points
#' @param facet_scales Scale option for facetting.  "Fixed" is default
#' @param order_genes If true, put genes in the same order as variable parameter
#' @param jitter_match If true, match jitter color to violin fill.
#' @param rasterize Whether to render the graphical layer as a raster image.  Default is FALSE.
#' @param raster_dpi If rasterize then this is the DPI used.  Default = 300.
#' @return A ggplot
#' @export
#' @import tidyverse monocle3 ggpubr
#' @importFrom ggrastr rasterise
bb_gene_violinplot <-
  function(cds,
           variable,
           genes_to_plot,
           pseudocount = 1,
           include_jitter = FALSE,
           ytitle = "Expression",
           plot_title = NULL,
           rows = 1,
           show_x_label = TRUE,
           legend_pos = "none",
           comparison_list = NULL,
           palette = NULL,
           violin_alpha = 1,
           jitter_alpha = 1,
           jitter_color = "black",
           jitter_fill = "transparent",
           jitter_size = 0.5,
           facet_scales = "fixed",
           order_genes = TRUE,
           jitter_match = FALSE,
           rasterize = FALSE,
           raster_dpi = 300
           ) {
    my_comparisons <-
      comparison_list#(list(c(comparator1,comparator2),c(comparator1,comparator3)...))
    if (length(dim(genes_to_plot)) > 1) {
      data_to_plot <-
        aggregate_gene_expression(cds = cds, gene_group_df = genes_to_plot) %>% as_tibble(rownames = "gene_group") %>% pivot_longer(
          cols = !gene_group,
          names_to = "barcode_sample",
          values_to = "expression"
        )
      data_to_plot <-
        colData(cds) %>% as_tibble(rownames = "barcode_sample") %>% left_join(data_to_plot) %>% mutate(expression = replace_na(expression, 0))
      p1 <-
        ggplot(data = data_to_plot, aes(x = !!as.name(variable),
                                        y = expression))
    } else {
      data_to_plot <-
        monocle3::plot_genes_violin(cds_subset = cds[rowData(cds)$gene_short_name %in% genes_to_plot,], group_cells_by = variable)[["data"]]

      if (order_genes) {
        data_to_plot <-
          data_to_plot %>% mutate(gene_short_name = factor(gene_short_name, levels = genes_to_plot))
        p1 <-
          ggplot(data = data_to_plot, aes(
            x = !!as.name(variable),
            y = log10(expression + pseudocount)
          ))
      }
    }
    if (include_jitter == TRUE) {
      if (jitter_match == TRUE){
       p1 <- p1 + geom_jitter(
         shape = 21,
         size = jitter_size,
         width = 0.2,
         alpha = jitter_alpha,
	 fill = jitter_fill,
         aes(color = !!as.name(variable)),
         show.legend = F
       )
      } else {
      p1 <-
        p1 + geom_jitter(
          shape = 21,
          size = jitter_size,
          color = jitter_color,
          fill = jitter_fill,
          alpha = jitter_alpha,
          width = 0.2
        )
      }
    }
    if (is.null(palette)) {
      p1 <- p1 +
        scale_color_viridis_d(alpha = jitter_alpha,
                             begin = 0.1,
                             end = 0.9)
    } else {
      p1 <- p1 +
        scale_color_manual(aesthetics = "color", values = alpha(palette, jitter_alpha), drop = TRUE)
    }
    p1 <- p1 +
      geom_violin(
        scale = "width",
        color = "black",
        trim = T,
        size = 0.5,
        aes(fill = !!as.name(variable)),
        draw_quantiles = 0.5
      ) +
      theme(legend.position = legend_pos) +
      theme(legend.direction = "horizontal") +
      theme(legend.justification = "center") +
      labs(
        x = "",
        y = ytitle,
        title = plot_title,
        fill = NULL
      )
    if (is.null(palette)) {
      p1 <- p1 +
        scale_fill_viridis_d(alpha = violin_alpha,
                             begin = 0.1,
                             end = 0.9)
    } else {
      p1 <- p1 +
        scale_fill_manual(values = alpha(palette, violin_alpha), drop = TRUE)
    }
    p1 <- p1 + theme(plot.title = element_text(hjust = 0.5))
    if (length(dim(genes_to_plot)) > 1) {
      p1 <- p1 +
        facet_wrap(~ gene_group, nrow = rows, scales = facet_scales) +
        theme(strip.background = element_rect(fill = "transparent"))
    } else {
      p1 <- p1 +
        facet_wrap(~ gene_short_name, nrow = rows, scales = facet_scales) +
        theme(strip.background = element_rect(fill = "transparent"))
    }
    if (show_x_label == F) {
      p1 <- p1 + theme(axis.text.x = element_blank())
    }
    # optionally rasterize the point layers
    if (rasterize) p1 <- ggrastr::rasterise(p1,
                                           dpi = raster_dpi)
    return(p1)
  }
