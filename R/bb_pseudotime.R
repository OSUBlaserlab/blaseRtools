#' @title Learn Graph and Calculate Pseudotime
#' @description This function determines a gene expression trajectory using `learn_graph` from monocle3 and then calculates pseudotime dimensions along this trajectory using `order_cells`.  So it is 2 functions wrapped into 1.  Usually we will not adjust the parameters for learn_graph with the possible exception of `close_loop` and `use_partition` which are also available in this function with the same defaults.  If you need to fine-tune the trajectory, use `monocle3::learn_graph` on the cds object first and then run this function to calculate pseudotime.  The graph learning will not be repeated on an object unless `force_graph` is set to `TRUE`.
#'
#' If you just want to look at the trajectory graph and not calculate pseudotime, change `calculate_pseudotime` to FALSE, or run `monocle3::learn_graph`.
#'
#' After the pseudotime values are calculated, they are handled differently than in monocle3.  In this function, they are copied from the hidden CDS slot and made an explicit cell metadata column.  Pseudotime needs a starting point or anchor.  There is no interactive option here as in monocle3.  To identify this starting point, you identify a cell metadata variable and provide it to `cluster_variable`.  This should identify a cohesive group of cells in UMAP space such as a leiden cluster, louvain cluster or partition.  Then provide a value corresponding to the cluster of interest to `cluster_value`.  The function will start pseudotime at the cell closest to the graph node in that cluster.  The pseudotime value column will be named automatically as a composite of the `cluster_variable` and `cluster_value` parameters.
#'
#' @param cds The cell data set object to calculate pseudotime upon.  Does not yet accept seurat objects.
#' @param calculate_pseudotime Logical, whether to calculate the pseudotime dimension.  If false, will only run learn_graph, Default: TRUE
#' @param cluster_variable The cell metadata column from which the pseudotime = 0 cell will be selected.
#' @param cluster_value The value of cluster_variable that identifies a cluster.  The cell closest to the root node closest to the center of this cluster will have pseudotime of 0.
#' @param use_partition Logical; If TRUE, learn_graph will construct trajectories within partitions.  If FALSE, it will connect partitions, Default: TRUE
#' @param close_loop Logical; Whether learn_graph will close looping trajectories, Default: TRUE
#' @param force_graph Logical; If TRUE, the function will recalculate the graph., Default: FALSE
#' @return A cell data set
#' @seealso
#'  \code{\link[cli]{cli_div}}, \code{\link[cli]{cli_alert}}
#'  \code{\link[monocle3]{learn_graph}}, \code{\link[monocle3]{order_cells}}, \code{\link[monocle3]{pseudotime}}
#'  \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @rdname bb_pseudotime
#' @export
#' @importFrom cli cli_div cli_alert_info
#' @importFrom monocle3 learn_graph order_cells pseudotime
#' @importFrom SummarizedExperiment colData
bb_pseudotime <-
  function(cds,
           calculate_pseudotime = TRUE,
           cluster_variable,
           cluster_value,
           use_partition = TRUE,
           close_loop = TRUE,
           force_graph = FALSE) {
    cli::cli_div(theme = list(span.emph = list(color = "orange")))
    # check and see if the principle graph needs to be calculated
    if (length(cds@principal_graph@listData) == 1) {
      cli::cli_alert_info("The graph has been learned.")
      if (force_graph) {
        cli::cli_alert_info("Re-learning graph")
        cds <-
          invisible(
            monocle3::learn_graph(
              cds = cds,
              use_partition = use_partition,
              close_loop = close_loop
            )
          )

      }
    } else {
      cli::cli_alert_info("Learning graph")
      cds <-
        invisible(
          monocle3::learn_graph(
            cds = cds,
            use_partition = use_partition,
            close_loop = close_loop
          )
        )
    }
    if (calculate_pseudotime) {
      pseudotime_colname <-
        paste0("pseudotime_", cluster_variable, "_", cluster_value)

      # identify root nodes
      cli::cli_alert_info(
        "Calculating pseudotime relative to the node nearest the center of {.emph {cluster_variable}} cluster {.emph {cluster_value}}."
      )
      cli::cli_alert_info(
        "Writing pseudotime dimensions to new cell metadata column {.emph {pseudotime_colname}}."
      )
      cds <-
        monocle3::order_cells(
          cds,
          root_pr_nodes = get_earliest_principal_node(
            cds,
            cluster_variable = cluster_variable,
            cluster_value = cluster_value
          )
        )
      SummarizedExperiment::colData(cds)[[pseudotime_colname]] <-
        monocle3::pseudotime(cds)

    }
    cds
  }

#' @importFrom SummarizedExperiment colData
#' @importFrom igraph V
#' @importFrom monocle3 principal_graph
get_earliest_principal_node <-
  function(cds, cluster_variable, cluster_value) {
    cell_ids <-
      which(SummarizedExperiment::colData(cds)[, cluster_variable] == cluster_value)

    closest_vertex <-
      cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds),])
    root_pr_nodes <-
      igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                          (which.max(table(closest_vertex[cell_ids, ]))))]

    root_pr_nodes
  }




#' Plots expression for one or more genes as a function of pseudotime
#'
#' @param cds Cell data set to plot.
#' @param gene_or_genes Gene or genes for which to plot pseudotime.
#' @param pseudotime_dim The column holding the pseudotime dimension to plot along.
#' @param min_expr the minimum (untransformed) expression level to plot.
#' @param cell_size the size (in points) of each cell used in the plot.
#' @param nrow the number of rows used when laying out the panels for each
#'   gene's expression.
#' @param ncol the number of columns used when laying out the panels for each
#'   gene's expression
#' @param panel_order vector of gene names indicating the order in which genes
#'   should be laid out (left-to-right, top-to-bottom). If
#'   \code{label_by_short_name = TRUE}, use gene_short_name values, otherwise
#'   use feature IDs.
#' @param color_cells_by the cell attribute (e.g. the column of colData(cds))
#'   to be used to color each cell.  Defaults to the value provided for pseudotime_dim.
#' @param trend_formula_df degrees of freedom for the model formula used to fit the expression
#'   trend over pseudotime.  The formulat takes the form of "~ splines::ns(pseudotime_dim, df = trend_formula_df)".
#' @param label_by_short_name label figure panels by gene_short_name (TRUE) or
#'   feature ID (FALSE).
#' @param vertical_jitter A value passed to ggplot to jitter the points in the
#'   vertical dimension. Prevents overplotting, and is particularly helpful for
#'   rounded transcript count data.
#' @param horizontal_jitter A value passed to ggplot to jitter the points in
#'   the horizontal dimension. Prevents overplotting, and is particularly
#'   helpful for rounded transcript count data.
#' @return a ggplot2 plot object
#' @export
#' @importFrom methods is
#' @importFrom cli cli_abort
#' @importFrom stringr str_detect
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom assertthat assert_that is.number is.count
#' @importFrom dplyr filter left_join
#' @importFrom monocle3 size_factors fit_models model_predictions
#' @importFrom SingleCellExperiment counts
#' @importFrom Matrix t
#' @importFrom reshape2 melt
#' @importFrom tibble as_tibble
#' @importFrom plyr ddply .
#' @importFrom ggplot2 ggplot aes_string geom_point position_jitter scale_color_viridis_c geom_line scale_y_log10 facet_wrap expand_limits ylab xlab labs theme element_blank
bb_plot_genes_in_pseudotime <- function(cds,
                                        gene_or_genes,
                                        pseudotime_dim,
                                        min_expr = NULL,
                                        cell_size = 0.75,
                                        nrow = NULL,
                                        ncol = 1,
                                        panel_order = NULL,
                                        color_cells_by = pseudotime_dim,
                                        trend_formula_df = 3,
                                        #"~ splines::ns(pseudotime, df=3)",
                                        label_by_short_name = TRUE,
                                        vertical_jitter = NULL,
                                        horizontal_jitter = NULL,
                                        legend_title = NULL) {
  if (!methods::is(cds, "cell_data_set"))
    cli::cli_abort("The object must be a cell_data_set for this function.")

  # check if we have a valid pseudotime dimension
  if (stringr::str_detect(pseudotime_dim, "pseudotime.*", negate = TRUE))
    cli::cli_abort("You must provide a valid pseudotime dimension to this function.")
  if (!is.numeric(SummarizedExperiment::colData(cds)[[pseudotime_dim]]))
    cli::cli_abort("You must provide a valid pseudotime dimension to identify the root node.")
  if (min(SummarizedExperiment::colData(cds)[[pseudotime_dim]]) != 0)
    cli::cli_abort("You must provide a valid pseudotime dimension to identify the root node.")
  if (!is.null(min_expr))
    assertthat::assert_that(assertthat::is.number(min_expr))
  assertthat::assert_that(assertthat::is.number(cell_size))
  if (!is.null(nrow))
    assertthat::assert_that(assertthat::is.count(nrow))
  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(is.logical(label_by_short_name))

  # subset the cds
  cds_subset <- filter_cds(cds,
                           genes = bb_rowmeta(cds) |>
                             dplyr::filter(gene_short_name %in% gene_or_genes))

  if (label_by_short_name) {
    assertthat::assert_that(
      "gene_short_name" %in% names(rowData(cds_subset)),
      msg = paste(
        "When label_by_short_name = TRUE,",
        "rowData must have a column of gene",
        "names called gene_short_name."
      )
    )
  }
  assertthat::assert_that(
    color_cells_by %in% c("cluster", "partition") |
      color_cells_by %in% names(colData(cds_subset)),
    msg = paste("color_cells_by must be a column in the",
                "colData table.")
  )

  if (!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(
        panel_order %in%
          SummarizedExperiment::rowData(cds_subset)$gene_short_name
      ))
    } else {
      assertthat::assert_that(all(
        panel_order %in%
          row.names(SummarizedExperiment::rowData(cds_subset))
      ))
    }
  }

  assertthat::assert_that(!is.null(monocle3::size_factors(cds_subset)))
  assertthat::assert_that(is.logical(label_by_short_name))
  assertthat::assert_that(
    color_cells_by %in% c("cluster", "partition") |
      color_cells_by %in% names(SummarizedExperiment::colData(cds_subset)),
    msg = paste("color_cells_by must be a column in the",
                "colData table.")
  )

  f_id <- NA
  Cell <- NA
  cds_subset = cds_subset[, is.finite(SummarizedExperiment::colData(cds_subset)[[pseudotime_dim]])]

  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <-
    Matrix::t(Matrix::t(cds_exprs) / monocle3::size_factors(cds_subset))
  cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))

  if (is.null(min_expr)) {
    min_expr <- 0
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_colData <- bb_cellmeta(cds_subset)
  cds_rowData <- bb_rowmeta(cds_subset)
  cds_exprs <-
    dplyr::left_join(tibble::as_tibble(cds_exprs),
                     cds_rowData,
                     by = c("f_id" = "feature_id"))
  cds_exprs <-
    dplyr::left_join(cds_exprs, cds_colData, by = c("Cell" = "cell_id"))

  cds_exprs$adjusted_expression <- cds_exprs$expression

  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <-
        cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$f_id <- as.character(cds_exprs$f_id)
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)


  new_data <- as.data.frame(colData(cds_subset))
  new_data$Size_Factor = 1
  trend_formula <-
    paste0("~ splines::ns(",
           pseudotime_dim,
           ", df=, ",
           trend_formula_df,
           ")")
  model_tbl = monocle3::fit_models(cds_subset, model_formula_str = trend_formula)

  model_expectation <- monocle3::model_predictions(model_tbl,
                                                  new_data = new_data)

  colnames(model_expectation) <- colnames(cds_subset)
  expectation <- plyr::ddply(cds_exprs, plyr::.(f_id, Cell),
                             function(x) {
                               data.frame("expectation" = model_expectation[x$f_id,
                                                                            x$Cell])
                             })
  cds_exprs <- merge(cds_exprs, expectation)

  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <-
    min_expr
  if (!is.null(panel_order)) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label,
                                      levels = panel_order)
  }
  q <-
    ggplot2::ggplot(ggplot2::aes_string(x = pseudotime_dim, y = "expression"),
                    data = cds_exprs) +
    ggplot2::geom_point(
      ggplot2::aes_string(color = color_cells_by),
      size = I(cell_size),
      position = ggplot2::position_jitter(horizontal_jitter,
                                          vertical_jitter)
    )

  if (methods::is(SummarizedExperiment::colData(cds_subset)[, color_cells_by], "numeric")) {
    q <- q + ggplot2::scale_color_viridis_c(option = "C")
  }
  q <-
    q + ggplot2::geom_line(ggplot2::aes_string(x = pseudotime_dim, y = "expectation"),
                           data = cds_exprs)

  q <-
    q + ggplot2::scale_y_log10() + ggplot2::facet_wrap( ~ feature_label,
                                                        nrow = nrow,
                                                        ncol = ncol,
                                                        scales = "free_y")
  if (min_expr < 1) {
    q <- q + ggplot2::expand_limits(y = c(min_expr, 1))
  }

  q <- q + ggplot2::ylab("Expression")

  q <- q + ggplot2::xlab("pseudotime")

  if (!is.null(legend_title))
    q <- q + ggplot2::labs(color = legend_title)

  q + ggplot2::theme(strip.background = ggplot2::element_blank())
}

