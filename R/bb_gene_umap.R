#' @title Make a Plot of Gene Expression in UMAP Form
#' @description Takes in a Seurat or cell_dat_set object, extracts UMAP dimensions and gene expression values.  For Seurat, default assay is "RNA"; can be changed to if necessary.  For cell_data_set, the assay parameter does nothing; the function extracts log and size-factor normalized counts which are similar but not identical to the Seurat "RNA" assay.  If a vector of genes is supplied to gene_or_genes, a faceted plot will be generated.  If a dataframe is supplied, an aggregated plot will be generated with a facet for each gene group.  The dataframe must be of 2 colums: the first containing feature ids and the second containing grouping information.  This is best generated using bb_rowmeta.
#' @param obj A Seurat or cell_data_set object.
#' @param gene_or_genes Individual gene or genes or aggregated genes to plot.  Supply a character string for a single gene, a vector for multiple genes or a dataframe for aggregated genes.  See description.
#' @param assay For Seurat objects only:  the gene expression assay to get expression data from, Default: 'RNA'
#' @param order Whether or not to order cells by gene expression.  When ordered, non-expressing cells are plotted first, i.e. on the bottom.  Caution:  when many cells are overplotted it may lead to a misleading presentation.  Generally bb_genebubbles is a better way to present, Default: TRUE
#' @param cell_size Size of the points, Default: 1
#' @param alpha Transparency of the points, Default: 1
#' @param ncol Specify the number of columns if faceting, Default: NULL
#' @param plot_title Optional title for the plot, Default: NULL
#' @param color_legend_title Option to change the color scale title, Default: 'Expression'
#' @param max_expr_val Maximum expression value to cap the color scale, Default: NULL
#' @param alt_dim_x Alternate/reference dimensions to plot by.
#' @param alt_dim_y Alternate/reference dimensions to plot by.
#' @param rasterize Whether to render the graphical layer as a raster image.  Default is FALSE.
#' @param raster_dpi If rasterize then this is the DPI used.  Default = 300.
#' @param cds Provided for backward compatibility.  If a value is supplied a warning will be emitted., Default: NULL
#' @return A ggplot
#' @seealso
#'  \code{\link[monocle3]{normalized_counts}}
#'  \code{\link[tibble]{as_tibble}}
#'  \code{\link[tidyr]{pivot_longer}}
#'  \code{\link[dplyr]{mutate-joins}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{arrange}}
#'  \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{aes}}, \code{\link[ggplot2]{geom_point}}, \code{\link[ggplot2]{scale_colour_viridis_d}}, \code{\link[ggplot2]{labs}}, \code{\link[ggplot2]{facet_wrap}}, \code{\link[ggplot2]{vars}}, \code{\link[ggplot2]{theme}}, \code{\link[ggplot2]{margin}}
#' @rdname bb_gene_umap
#' @export
#' @importFrom monocle3 normalized_counts
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join mutate select arrange
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_c scale_fill_viridis_c labs facet_wrap vars theme element_blank element_text
#' @importFrom ggrastr rasterise
bb_gene_umap <-
  function (obj,
            gene_or_genes,
            assay = "RNA",
            order = TRUE,
            cell_size = 1,
            alpha = 1 ,
            ncol = NULL,
            plot_title = NULL,
            color_legend_title = "Expression",
            max_expr_val = NULL,
            alt_dim_x = NULL,
            alt_dim_y = NULL,
            rasterize = FALSE,
            raster_dpi = 300,
            cds = NULL) {
    cds_warn(cds)
    obj_stop(obj)

    # use alternate dimensions if desired, otherwise use UMAP
    dim_x <- ifelse(is.null(alt_dim_x), "UMAP_1", alt_dim_x)
    dim_y <- ifelse(is.null(alt_dim_y), "UMAP_2", alt_dim_y)

    # check to see if you are dealing with single or aggregated genes
    aggregated <-
      ifelse(length(dim(gene_or_genes)) == 2, TRUE, FALSE)

    # next get the data you will need for each gene in the list

    if (!aggregated) {
      # convert gene_or_genes into feature_id
      ids <- get_gene_ids(obj, gene_or_genes)

      if ("Seurat" %in% class(obj)) {
        dat <- obj[[assay]]@data
        dat <- dat[rownames(dat) %in% ids,]
        if (class(dat) == "numeric") {
          dat <- as.matrix(dat)
          colnames(dat) <- ids
        } else {
          dat <- t(as.matrix(dat))
        }
        if (assay == "RNA") {
          dat[dat == 0] <- NA
        }
      } else {
        dat <- monocle3::normalized_counts(obj)
        dat <- dat[rownames(dat) %in% ids,]
        if (class(dat) == "numeric") {
          dat <- as.matrix(dat)
          colnames(dat) <- ids
        } else {
          dat <- t(as.matrix(dat))
        }
        dat[dat == 0] <- NA
      }
    } else {
      dat <-
        bb_aggregate(obj,
                     gene_group_df = gene_or_genes,
                     scale_agg_values = F
                     )
      dat <- t(as.matrix(dat))
    }
    dat <- tibble::as_tibble(dat, rownames = "cell_id") |>
      tidyr::pivot_longer(-cell_id)

    # always replace ensemble ids with gene short names
    if (!aggregated) {
      dat <-
        dplyr::left_join(dat, bb_rowmeta(obj), by = c("name" = "feature_id")) |>
        dplyr::mutate(name = ifelse(is.na(gene_short_name), name, gene_short_name)) |>
        dplyr::select(cell_id, name, value)
    }

    # get the plot data
    if ("cell_data_set" %in% class(obj)) {
      plot_data <- get_cds_umap_dims(obj)
    } else {
      plot_data <- get_seurat_umap_dims(obj, assay)
    }

    # join the expression data back onto the plot data
    plot_data <- dplyr::left_join(plot_data, dat, by = "cell_id")

    # join the cell metadata back on in case you want to facet or something else
    plot_data <- dplyr::left_join(plot_data, bb_cellmeta(obj), by = "cell_id")

    # optionally order the cells to un-bury rare expressing cells
    if (order)
      plot_data <- dplyr::arrange(plot_data, !is.na(value), value)

    # optionally cap the color scale

    if (!is.null(max_expr_val)) {
      plot_data <- plot_data |>
        mutate(value = ifelse(value > max_expr_val, max_expr_val, value))
    }

    p <- ggplot2::ggplot(plot_data,
                mapping = ggplot2::aes(
                  x = !!sym(dim_x),
                  y = !!sym(dim_y),
                  color = value,
                  fill = value
                )) +
      ggplot2::geom_point(shape = 16,
                 size = cell_size,
                 stroke = 0.25) +
      ggplot2::scale_color_viridis_c(na.value = "grey80") +
      ggplot2::scale_fill_viridis_c(na.value = "transparent", guide = "none") +
      ggplot2::labs(
        x = ifelse(is.null(alt_dim_x), "UMAP 1", alt_dim_x),
        y = ifelse(is.null(alt_dim_y), "UMAP 2", alt_dim_y),
        color = color_legend_title,
        title = plot_title
      ) +

      ggplot2::facet_wrap(facets = ggplot2::vars(name), ncol = ncol) +
      ggplot2::theme(strip.background = ggplot2::element_blank()) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    # optionally rasterize the point layers
    if (rasterize) p <- ggrastr::rasterise(p,
                                           layers = "Point",
                                           dpi = raster_dpi)
    return(p)
  }


#' @title Aggregate Single Cell Gene Expression
#' @description Generates a matrix of counts aggregated by gene and/or cell group.
#' @param obj A Seurat or cell data set object
#' @param assay Gene expression assay to use for aggregation; currently only applies to Seurat objects, Default: 'RNA'
#' @param gene_group_df A 2-column dataframe with gene names or ids and gene groupings, Default: NULL
#' @param cell_group_df A 2-coumn dataframe with cell ids and gene groupings, Default: NULL
#' @param norm_method Gene normalization method, Default: c("log", "binary", "size_only")
#' @param pseudocount Pseudocount, Default: 1
#' @param scale_agg_values Whether to scale the aggregated values, Default: TRUE
#' @param max_agg_value If scaling, make this the maximum aggregated value, Default: 3
#' @param min_agg_value If scaling, make this the minimum aggregated value, Default: -3
#' @param binary_min Minimum value below which a cell is considered not to express a feature, Default: 0
#' @param exclude.na Exclude NA?, Default: TRUE
#' @return A dense or sparse matrix.
#' @details The best way to group genes or cells is by using bb_*meta and then select cell_id or feature_id plus one metadata column with your group labels.
#' @rdname bb_aggregate
#' @seealso
#'  \code{\link[cli]{cli_div}}, \code{\link[cli]{cli_alert}}
#'  \code{\link[monocle3]{normalized_counts}}, \code{\link[monocle3]{my.aggregate.Matrix}}
#'  \code{\link[Matrix]{character(0)}}
#' @export
#' @importFrom cli cli_div cli_alert_info
#' @import monocle3
#' @importFrom Matrix t
bb_aggregate <-
  function (obj,
            assay = "RNA",
            experiment_type = "Gene Expression",
            gene_group_df = NULL,
            cell_group_df = NULL,
            norm_method = c("log",
                            "binary",
                            "size_only"),
            pseudocount = 1,
            scale_agg_values = TRUE,
            max_agg_value = 3,
            min_agg_value = -3,
            binary_min = 0,
            exclude.na = TRUE) {
    norm_method <- match.arg(norm_method)
    obj_stop(obj)
    if ("cell_data_set" %in% class(obj)) {
      # check to be sure experiment_type is available
      all_exps <- c(
        SingleCellExperiment::mainExpName(obj),
        SingleCellExperiment::altExpNames(obj)
      )
    if (!is.null(SingleCellExperiment::mainExpName(obj))) {
      if (experiment_type %notin% all_exps)
        cli::cli_abort("The requested experiment name is not available.")
      if (experiment_type != "Gene Expression") {
        obj <-
          as(SingleCellExperiment::swapAltExp(obj, name = experiment_type),
             Class = "cell_data_set")
      }

    }

    }
    # set the cli alert aesthetics
    cli::cli_div(theme = list(span.emph = list(color = "orange")))

    # check to be sure the grouping was properly assigned
    if (is.null(gene_group_df) && is.null(cell_group_df))
      stop("Error: one of either gene_group_df or cell_group_df must not be NULL")

    # get the normalized coutns from the object using method appropriate to object class
    if ("cell_data_set" %in% class(obj)) {
      if (experiment_type == "Gene Expression") {
        cli::cli_alert_info("Getting normalized {.emph gene expression} counts from the monocle cell_data_set.")
        agg_mat <- monocle3::normalized_counts(obj,
                                             norm_method = norm_method,
                                             pseudocount = pseudocount)
      } else if (experiment_type == "Antibody Capture") {
        cli::cli_alert_info("Getting CLR-normalized {.emph antibody capture} counts from the monocle cell_data_set.")
        agg_mat <- SummarizedExperiment::assay(obj, i = "CLR_counts")
      } else
        cli::cli_abort("The requested experiment assay is not available.")
    } else if ("Seurat" %in% class(obj)) {
      cli::cli_alert_info("Getting normalized counts from the {.emph Seurat} {.emph {assay}} assay")
      agg_mat <- obj[[assay]]@data
    }
    if (norm_method == "binary")
      agg_mat <- agg_mat > binary_min

    if (!is.null(gene_group_df)) {
      gene_aggregator <- colnames(gene_group_df)[2]
      cli::cli_alert_info("Aggregating features by {.emph {gene_aggregator}}")
      gene_group_df <- as.data.frame(gene_group_df)
      gene_group_df <- gene_group_df[gene_group_df[, 1] %in%
                                       bb_rowmeta(obj)$gene_short_name |
                                       gene_group_df[, 1] %in%
                                       bb_rowmeta(obj)$feature_id, , drop = FALSE]
      short_name_mask <-
        gene_group_df[[1]] %in% bb_rowmeta(obj)$gene_short_name
      if (any(short_name_mask)) {
        geneids <- as.character(gene_group_df[[1]])
        geneids[short_name_mask] <-
          bb_rowmeta(obj)$feature_id[match(geneids[short_name_mask],
                                           bb_rowmeta(obj)$gene_short_name)]
        gene_group_df[[1]] <- geneids
      }
      agg_mat = as.matrix(monocle3:::my.aggregate.Matrix(agg_mat[gene_group_df[,
                                                                               1],], as.factor(gene_group_df[, 2]), fun = "sum"))
      if (scale_agg_values) {
        agg_mat <- t(scale(t(agg_mat)))
        agg_mat[agg_mat < min_agg_value] <- min_agg_value
        agg_mat[agg_mat > max_agg_value] <- max_agg_value
      }
    }
    if (!is.null(cell_group_df)) {
      cell_aggregator <- colnames(cell_group_df)[2]
      cli::cli_alert_info("Aggregating cells by {.emph {cell_aggregator}}")
      cell_group_df <- as.data.frame(cell_group_df)
      cell_group_df <- cell_group_df[cell_group_df[, 1] %in%
                                       bb_cellmeta(obj)$cell_id, , drop = FALSE]
      agg_mat <- agg_mat[, cell_group_df[, 1]]
      agg_mat <-
        monocle3:::my.aggregate.Matrix(Matrix::t(agg_mat), as.factor(cell_group_df[,
                                                                                   2]), fun = "mean")
      agg_mat <- Matrix::t(agg_mat)
    }
    if (exclude.na) {
      agg_mat <- agg_mat[rownames(agg_mat) != "NA", colnames(agg_mat) !=
                           "NA", drop = FALSE]
    }
    return(agg_mat)
  }


#' @importFrom dplyr mutate case_when filter
get_gene_ids <- function(obj, gen, et = "Gene Expression") {
  rowmet <- bb_rowmeta(obj, experiment_type = et)
  if (!"gene_short_name" %in% colnames(rowmet)) {
    cli::cli_alert_warning("Copying feature_id to a new column, gene_short_name.")
    cli::cli_alert_warning("You may want to change your object to have a proper gene_short_name metadata column.")
    rowmet <- rowmet |> dplyr::mutate(gene_short_name = feature_id)
  }
  rowmet |>
    dplyr::mutate(
      ids =
        dplyr::case_when(
          gene_short_name %in% gen ~ feature_id,
          feature_id %in% gen ~ feature_id,
          TRUE ~ "other"
        )
    ) |>
    dplyr::filter(ids != "other") |>
    pull(ids)
}
