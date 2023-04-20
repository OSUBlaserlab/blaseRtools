#' @title Cluster Representation By Regression Per Sample
#' @description Use this function to determine the differential representation of cells in clusters. It uses a regression method to determine fold change between groups of biological samples.  It can only compare two sample groups, e.g. control vs experimental at this point.  See parameter descriptions for how to identify these properly.
#' @param obj The (possibly filtered) single cell object to operate on.  Can be either Seurat or monocle/CDS object.
#' @param sample_var The metadata column holding the biological sample information.
#' @param cluster_var The metadata column holding the clustering or other cell classification information.
#' @param comparison_var The metadata column holding the comparison group information.  There can be only two levels in this column.  Character data will be converted to factors.
#' @param comparison_levels A character vector identifying the order of the levels to compare. The first value will be shown with negative log2Fold Change and the second will be positive.  If NULL (default), R will pick for you.
#' @param color_pal Color palette for the comparison levels, Default: c("red3", "blue4")
#' @param sig_val Report PValue or FDR, Default:  "FDR"
#' @param return_val Value to return, Default: c("plot", "table)
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @source http://bioconductor.org/books/3.13/OSCA.multisample/differential-abundance.html
#' @rdname bb_cluster_representation2
#' @export
#' @importFrom dplyr count mutate case_when
#' @importFrom tidyr pivot_wider
#' @importFrom edgeR DGEList estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom tibble as_tibble
#' @importFrom forcats fct_reorder
#' @importFrom ggtext element_markdown
#' @import ggplot2
bb_cluster_representation2 <-
  function(obj,
           sample_var,
           cluster_var,
           comparison_var,
           comparison_levels = NULL,
           color_pal = c("red3", "blue4"),
           sig_val = c("FDR","PValue"),
           return_val = c("plot", "data")) {
    return_val <- match.arg(return_val)
    sig_val <- match.arg(sig_val)
    cli::cli_div(theme = list(span.emph = list(color = "orange")))
    cli::cli_alert_info("Reporting significance as {.emph {sig_val}}.")
    # get counts by cluster
    counts <- bb_cellmeta(obj) |>
      dplyr::count(!!sym(sample_var), !!sym(cluster_var)) |>
      tidyr::pivot_wider(
        names_from = sample_var,
        values_from = "n",
        values_fill = 0
      ) |>
      bb_tbl_to_matrix()

    dge_list <- edgeR::DGEList(
      counts = counts,
      samples = bb_cellmeta(obj) |>
        group_by(!!sym(sample_var),
                 !!sym(comparison_var)) |>
        summarise()
    )

    if (is.null(comparison_levels))
      comparison_levels <-
      unique(dge_list$samples[[comparison_var]])

    dge_list$samples[[comparison_var]] <-
      factor(dge_list$samples[[comparison_var]], levels = comparison_levels)
    design <-
      model.matrix( ~ dge_list$samples[[comparison_var]], dge_list$samples, )
    dge_list <-
      edgeR::estimateDisp(dge_list, design, trend = "none")
    fit <-
      edgeR::glmQLFit(dge_list,
                      design,
                      robust = TRUE,
                      abundance.trend = FALSE)
    res <- edgeR::glmQLFTest(fit, coef = ncol(design))
    res_top_tags <- edgeR::topTags(res, n = Inf)
    data <-
      tibble::as_tibble(res_top_tags@.Data[[1]], rownames = "cluster") |>
      dplyr::mutate(enriched = ifelse(logFC < 0, comparison_levels[1], comparison_levels[2])) |>
      dplyr::mutate(sig_col = !!sym(sig_val)) |>
      dplyr::mutate(
        sig = dplyr::case_when(
          sig_col < 0.05 & sig_col >= 0.01 ~ "*",
          sig_col < 0.01 & sig_col >= 0.001 ~ "**",
          sig_col < 0.001 & sig_col >= 0.0001 ~ "***",
          sig_col < 0.0001 ~ "****",
          sig_col >= 0.05 ~ ""
        )
      ) |>
      dplyr::select(-sig_col) |>
      dplyr::mutate(cluster = forcats::fct_reorder(cluster, logFC)) |>
      dplyr::mutate(texty = ifelse(logFC > 0,
                                   logFC,
                                   0))

    if (return_val == "data")
      return(data)

    ggplot(data, aes(
      x = cluster,
      y = logFC,
      fill = enriched,
      label = sig
    )) +
      geom_col(color = "black") +
      scale_fill_manual(values = color_pal) +
      labs(x = cluster_var,
           y = "Log<sub>2</sub>Fold Enrichment",
           fill = "Enriched Group") +
      theme(axis.title.y.left = ggtext::element_markdown()) +
      geom_text(aes(y = texty), color = "black", nudge_y = 0.05) +
      theme(legend.position = "right")

  }
