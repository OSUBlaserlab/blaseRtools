#' A function to generate a UMAP with colors mapped to colData variables
#'
#' @param obj A Seurat or cell data set object
#' @param var The variable to map colors to.  Special exceptions are "density", "local_n" and "log_local_n" which calculate the 2 d kernel density estimate or binned cell counts and maps to color scale.
#' @param assay The gene expression assay to draw reduced dimensions from.  Default is "RNA".  Does not do anything with cell_data_set objects.
#' @param value_to_highlight Option to highlight a single value
#' @param foreground_alpha Alpha value for foreground points
#' @param legend_pos Legend position
#' @param cell_size Cell point size
#' @param alt_stroke_color Alternative color for the data point stroke
#' @param legend_title Title for the legend
#' @param plot_title Main title for the plot
#' @param palette Color palette to use.  "Rcolorbrewer", "Viridis" are builtin options.  Otherwise provide manual values.
#' @param alt_dim_x Alternate/reference dimensions to plot by.
#' @param alt_dim_y Alternate/reference dimensions to plot by.
#' @param overwrite_labels Whether to overwrite the variable value labels
#' @param group_label_size Size of the overwritten labels
#' @param alt_label_col Alternate column to label cells by
#' @param shape Shape for data points
#' @param nbin Number of bins if using var %in% c("density". "local_n", "log_local_n")
#' @param facet_by Variable or variables to facet by.
#' @param sample_equally Whether or not you should downsample to the same number of cells in each plot.  Default is FALSE or no.
#' @param rasterize Whether to render the graphical layer as a raster image.  Default is FALSE.
#' @param raster_dpi If rasterize then this is the DPI used.  Default = 300.
#' @param ... Additional params for facetting.
#' @param man_text_df A data frame in the form of text_x = numeric_vector, text_y = numeric_vector, label = character_vector for manually placing text labels.
#' @param show_trajectory_graph Whether to render the principal graph for the
#'   trajectory. Requires that learn_graph() has been called on cds.
#' @param trajectory_graph_color The color to be used for plotting the
#'   trajectory graph.
#' @param trajectory_graph_segment_size The size of the line segments used for
#'   plotting the trajectory graph.
#' @param graph_label_size How large to make the branch, root, and leaf labels.
#' @param label_root_node Logical; whether to label the root node for the selected pseudotime trajectory.  The function will requires that a valid pseudotime column be identified, usually as the value of the "var" argument in the form of "pseudotime_cluster_value".  If you wish to use var to color the cells in some other way, the pseudotime_dim argument needs to be supplied with the correct pseudotime dimension to pick the root node from.
#' @param pseudotime_dim An alternative column to pick the pseudoetime root node from, if not supplied to var.
#' @param label_principal_points Logical indicating whether to label roots,
#'   leaves, and branch points with principal point names. This is useful for
#'   order_cells and choose_graph_segments in non-interactive mode.
#' @param cds Provided for backward compatibility with prior versions.  If a value is supplied, a warning will be emitted and the value will be transferred to the obj argument, Default: NULL
#' @return a ggplot
#' @rdname bb_var_umap
#' @export
#' @importFrom dplyr left_join group_by summarise slice_sample ungroup n mutate filter pull
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_point aes aes_string scale_color_viridis_c scale_fill_viridis_c scale_fill_viridis_d scale_color_viridis_d scale_color_brewer scale_fill_brewer scale_color_manual scale_fill_manual scale_color_discrete scale_fill_discrete guides guide_legend theme element_text geom_text facet_wrap element_blank facet_grid geom_label
#' @importFrom purrr map_dfr
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @importFrom ggrastr rasterise
bb_var_umap <- function(obj,
                        var,
                        assay = "RNA",
                        value_to_highlight = NULL,
                        foreground_alpha = 1,
                        legend_pos = "right",
                        cell_size = 0.5,
                        alt_stroke_color = NULL,
                        legend_title = NULL,
                        plot_title = NULL,
                        palette = NULL,
                        alt_dim_x = NULL,
                        alt_dim_y = NULL,
                        overwrite_labels = FALSE,
                        group_label_size = 3,
                        alt_label_col = NULL,
                        shape = 21,
                        nbin = 100,
                        facet_by = NULL,
                        sample_equally = FALSE,
                        rasterize = FALSE,
                        raster_dpi = 300,
                        show_trajectory_graph = FALSE,
                        trajectory_graph_color = "grey28",
                        trajectory_graph_segment_size = 0.75,
                        label_root_node = FALSE,
                        pseudotime_dim = var,
                        label_principal_points = FALSE,
                        graph_label_size = 2,
                        cds = NULL,
                        test = TRUE,
                        ...,
                        man_text_df = NULL,
                        text_geom = "text") {
  cds_warn(cds)
  obj_stop(obj)

  if (text_geom == "label") {
    text_geom <- geom_label
    text_geom_repel <- geom_label_repel
  } else {
    text_geom <- geom_text
    text_geom_repel <- geom_text_repel
  }

  if ("cell_data_set" %in% class(obj)) {
    dims <- get_cds_umap_dims(obj)
  } else if ("Seurat" %in% class(obj)) {
    dims <- get_seurat_umap_dims(obj, assay)
  } else {
    stop("You must use this function with a cell_data_set or Seurat object")
  }

  cellmeta <- bb_cellmeta(obj)
  plot_data <- dplyr::left_join(dims, cellmeta, by = "cell_id")

  if (sample_equally) {
    if (is.null(facet_by)) {
      return("You need to provide a faceting variable.")
    } else {
      if (length(facet_by) == 1) {
        sample_sizes <-
          plot_data |>
          dplyr::group_by(!!sym(facet_by)) |>
          dplyr::summarise(ncells = n())
        plot_data <-
          plot_data |>
          dplyr::group_by(!!sym(facet_by)) |>
          dplyr::slice_sample(n = min(sample_sizes$ncells)) |>
          dplyr::ungroup()


      } else if (length(facet_by) == 2) {
        # create a composite variable to group by
        plot_data$cross_product <-
          paste0(plot_data[, facet_by[1]], "_", plot_data[, facet_by[2]])
        sample_sizes <-
          plot_data |>
          dplyr::group_by(cross_product) |>
          dplyr::summarise(ncells = n())
        plot_data <-
          plot_data |>
          dplyr::group_by(cross_product) |>
          dplyr::slice_sample(n = min(sample_sizes$ncells)) |>
          dplyr::ungroup()
      } else if (length(facet_by) > 2) {
        return("You can only facet by 2 variables")
      }
    }
  }

  if (var %in% c("density", "local_n", "log_local_n")) {
    data_long <- plot_data
  } else {
    data_long <-
      plot_data |> tidyr::pivot_longer(cols = (!!sym(var)), names_to = "var")
  }
  dim_x <- ifelse(is.null(alt_dim_x), "UMAP_1", alt_dim_x)
  dim_y <- ifelse(is.null(alt_dim_y), "UMAP_2", alt_dim_y)

  # generate text data frame for variable labels if you are going to use them
  if (is.null(alt_label_col)) {
    if (var %in% c("density", "local_n", "log_local_n")) {
      text_df <- data_long
    } else {
      text_df <- data_long |> dplyr::group_by(value)
    }

  } else {
    text_df <-
      plot_data |> tidyr::pivot_longer(cols = !!sym(alt_label_col),
                                       names_to = "var") |> dplyr::group_by(value)
  }
  if (overwrite_labels == T && is.null(man_text_df)) {
    median_coord_df <-
      text_df |>
      dplyr::summarise(
        fraction_of_group = dplyr::n(),
        text_x = median(!!sym(dim_x)),
        text_y = median(!!sym(dim_y))
      )
    text_df <- dplyr::left_join(text_df, median_coord_df) |>
      dplyr::mutate(label = value)
    text_df <-
      text_df |> dplyr::group_by(label, text_x, text_y) |> dplyr::summarise()
  }

  if (!is.null(man_text_df))
    text_df <- man_text_df

  # make the main plot
  plot <- ggplot2::ggplot()
  if (!is.null(value_to_highlight)) {
    data_background <-
      data_long |> dplyr::filter(value %notin% value_to_highlight)
    data_long <-
      data_long |> dplyr::filter(value %in% value_to_highlight)
    plot <- plot +
      ggplot2::geom_point(
        data = data_background,
        ggplot2::aes(x = !!sym(dim_x),
                     y = !!sym(dim_y)),
        stroke = 0.25,
        shape = 1,
        size = cell_size,
        color = "grey80"
      )
  }

  #option to color dots by local density
  if (var == "density") {
    if (!is.null(facet_by)) {
      if (length(facet_by) == 1) {
        data_long <-
          purrr::map_dfr(
            .x = data_long |> dplyr::group_by(!!sym(facet_by)) |> dplyr::summarise() |> dplyr::pull() |> as.character(),
            .f = function(x, density_data = data_long) {
              density_data <-
                density_data |> dplyr::filter(!!sym(facet_by) == x) |> as.data.frame()
              density_data$density <-
                get_density(xgd = density_data[, dim_x],
                            ygd = density_data[, dim_y],
                            ngd = nbin)
              return(density_data)
            }
          )
      } else if (length(facet_by) == 2) {
        data_long$cross_product <-
          paste0(data_long[, facet_by[1]], "_", data_long[, facet_by[2]])
        data_long <-
          purrr::map_dfr(
            .x = data_long$cross_product |> unique(),
            .f = function(x, density_data = data_long) {
              density_data <-
                density_data |> dplyr::filter(cross_product == x) |> as.data.frame()
              density_data$density <-
                get_density(xgd = density_data[, dim_x],
                            ygd = density_data[, dim_y],
                            ngd = nbin)
              return(density_data)
            }
          )
      }
    } else {
      data_long <-
        data_long |> dplyr::mutate(density = get_density(x = data[, dim_x], y = data[, dim_y], n = nbin))
    }

    plot <- ggplot2::ggplot(data_long) +
      ggplot2::geom_point(
        ggplot2::aes_string(
          x = dim_x,
          y = dim_y,
          color = "density",
          fill = "density"
        ),
        size = cell_size,
        stroke = 0.25,
        shape = shape,
        alpha = foreground_alpha
      ) +
      ggplot2::scale_color_viridis_c(option = "inferno",
                                     begin = 0.1,
                                     end = 0.9) +
      ggplot2::scale_fill_viridis_c(
        option = "inferno",
        begin = 0.1,
        end = 0.9,
        guide = "none"
      )
  } else if (var %in% c("local_n", "log_local_n")) {
    if (!is.null(facet_by)) {
      if (length(facet_by) == 1) {
        data_long <-
          purrr::map_dfr(
            .x = data_long |> dplyr::pull(!!sym(facet_by)) |> unique(),
            .f = function(x, data = data_long) {
              data <- data |> dplyr::filter(!!sym(facet_by) == x)
              data <-
                data |> dplyr::mutate(local_n = get_hexcount(
                  data = data,
                  x = dim_x,
                  y = dim_y,
                  n = nbin
                ))
              return(data)
            }
          )
      } else if (length(facet_by) == 2) {
        data_long$cross_product <-
          paste0(data_long[, facet_by[1]], "_", data_long[, facet_by[2]])
        data_long <-
          purrr::map_dfr(
            .x = data_long$cross_product |> unique(),
            .f = function(x, data = data_long) {
              data <- data |> dplyr::filter(cross_product == x)
              data <-
                data |> dplyr::mutate(local_n = get_hexcount(
                  data = data,
                  x = dim_x,
                  y = dim_y,
                  n = nbin
                ))
              return(data)
            }
          )
      }
    } else {
      data_long <-
        data_long |> dplyr::mutate(local_n = get_hexcount(
          data = data_long,
          x = dim_x,
          y = dim_y,
          n = nbin
        ))
    }
    data_long <-
      data_long |> dplyr::mutate(log_local_n = log10(local_n))
    plot <- ggplot2::ggplot(data_long) +
      ggplot2::geom_point(
        ggplot2::aes_string(
          x = dim_x,
          y = dim_y,
          color = var,
          fill = var
        ),
        size = cell_size,
        stroke = 0.25,
        shape = shape,
        alpha = foreground_alpha
      ) +
      ggplot2::scale_color_viridis_c(option = "inferno",
                                     begin = 0.1,
                                     end = 0.9) +
      ggplot2::scale_fill_viridis_c(
        option = "inferno",
        begin = 0.1,
        end = 0.9,
        guide = "none"
      )
  } else {
    plot <- plot +
      ggplot2::geom_point(
        data = data_long,
        ggplot2::aes(
          x = !!sym(dim_x),
          y = !!sym(dim_y),
          fill = value,
          color = value
        ),
        stroke = 0.25,
        shape = shape,
        alpha = foreground_alpha,
        size = cell_size
      )
    if (class(data_long$value) == "numeric") {
      plot <- plot +
        ggplot2::scale_fill_viridis_c(guide = "colorbar", na.value = "transparent") +
        ggplot2::scale_color_viridis_c(guide = "none", na.value = "grey80")
    } else if (length(palette) == 1 && palette == "viridis") {
      plot <- plot +
        ggplot2::scale_fill_viridis_d(begin = 0.1, end = 0.9) +
        ggplot2::scale_color_viridis_d(begin = 0.1,
                                       end = 0.9,
                                       guide = "none")
    } else if (length(palette) == 1 && palette == "rcolorbrewer") {
      plot <-
        plot + ggplot2::scale_color_brewer(palette = "Paired", guide = "none") +
        ggplot2::scale_fill_brewer(palette = "Paired")
    } else if (!is.null(palette)) {
      plot <-
        plot + ggplot2::scale_color_manual(values = palette, guide = "none") +
        ggplot2::scale_fill_manual(values = palette)
    } else {
      plot <- plot + ggplot2::scale_color_discrete(guide = "none") +
        ggplot2::scale_fill_discrete()
    }
    if (class(data_long$value) != "numeric") {
      plot <-
        plot + ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(
          size = 2,
          alpha = 1,
          color = "transparent"
        )))
    }
  }

  plot <- plot + labs(
    x = ifelse(is.null(alt_dim_x), "UMAP 1", alt_dim_x),
    y = ifelse(is.null(alt_dim_y), "UMAP 2", alt_dim_y),
    title = plot_title,
    fill = legend_title
  ) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  plot <- plot +
    ggplot2::theme(legend.position = legend_pos)#+coord_fixed()

  #option to overwrite labels
  if (overwrite_labels == T && is.null(man_text_df)) {
    plot <- plot +
      ggplot2::theme(legend.position = "none") +
      text_geom_repel(
        data = text_df,
        mapping = ggplot2::aes_string(x = "text_x", y = "text_y", label = "label"),
        size = group_label_size,
        min.segment.length = 1
      )
  } else if (!is.null(man_text_df)) {
    plot <- plot +
      ggplot2::theme(legend.position = "none") +
      text_geom(
        data = text_df,
        mapping = ggplot2::aes_string(x = "text_x", y = "text_y", label = "label"),
        size = group_label_size
      )
  }
  if (!is.null(facet_by)) {
    if (length(facet_by) == 1) {
      plot <- plot +
        ggplot2::facet_wrap(facets = facet_by, ...) +
        ggplot2::theme(strip.background = ggplot2::element_blank())
    } else if (length(facet_by) == 2) {
      plot <- plot +
        ggplot2::facet_grid(...) +
        ggplot2::theme(strip.background = ggplot2::element_blank())
    } else {
      return("Too many dimensions to facet by.")
    }
  }

  if (show_trajectory_graph)
    plot <- make_trajectory_graph(
      g = plot,
      cds = obj,
      trajectory_graph_segment_size = trajectory_graph_segment_size,
      trajectory_graph_color = trajectory_graph_color
    )

  if (label_principal_points)
    plot <-
    make_pp_labels(
      g = plot,
      cds = obj,
      graph_label_size = graph_label_size,
      trajectory_graph_segment_size = trajectory_graph_segment_size
    )

  if (label_root_node) {
    # check if we have a valid pseudotime dimension
    if (stringr::str_detect(pseudotime_dim, "pseudotime.*", negate = TRUE))
      cli::cli_abort("You must provide a valid pseudotime dimension to identify the root node.")
    if (!is.numeric(SummarizedExperiment::colData(obj)[[pseudotime_dim]]))
      cli::cli_abort("You must provide a valid pseudotime dimension to identify the root node.")
    if (min(SummarizedExperiment::colData(obj)[[pseudotime_dim]]) != 0)
      cli::cli_abort("You must provide a valid pseudotime dimension to identify the root node.")
    plot <- make_root_node_labels(g = plot,
                                  cds = obj,
                                  pseudotime_dim = pseudotime_dim,
                                  trajectory_graph_segment_size = trajectory_graph_segment_size,
                                  graph_label_size = graph_label_size)


  }


  # if (label_roots)
  #   plot <-
  #   make_root_labels(
  #     g = plot,
  #     cds = obj,
  #     trajectory_graph_segment_size = trajectory_graph_segment_size,
  #     graph_label_size = graph_label_size
  #   )

  # if (label_leaves)
  #   plot <-
  #   make_leaf_labels(
  #     g = plot,
  #     cds = obj,
  #     trajectory_graph_segment_size = trajectory_graph_segment_size,
  #     graph_label_size = graph_label_size
  #   )

  # if (label_branch_points)
  #   plot <-
  #   make_branch_labels(
  #     g = plot,
  #     cds = obj,
  #     graph_label_size = graph_label_size,
  #     trajectory_graph_segment_size = trajectory_graph_segment_size
  #   )

  # optionally rasterize the point layers
  if (rasterize)
    plot <- ggrastr::rasterise(plot,
                               layers = "Point",
                               dpi = raster_dpi)
  # finally trim the plot
  plot <- plot +
    xlim(range(obj@int_colData@listData$reducedDims@listData$UMAP[, 1])) +
    ylim(range(obj@int_colData@listData$reducedDims@listData$UMAP[, 2]))

  return(plot)
}

#' @importFrom MASS kde2d
get_density <- function(xgd, ygd, ngd) {
  dens <- MASS::kde2d(x = xgd, y = ygd, n = ngd)
  ix <- findInterval(xgd, dens$x)
  iy <- findInterval(ygd, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#' @importFrom hexbin hexbin
#' @importFrom dplyr pull mutate left_join
#' @importFrom tibble tibble
get_hexcount <- function(data, x, y, n) {
  hexdata <- hexbin::hexbin(
    x = data |> dplyr::pull(!!sym(x)),
    y = data |> dplyr::pull(!!sym(y)),
    IDs = T,
    xbins = n
  )
  hexdata_cells_counts <-
    tibble::tibble(hexcell = as.character(hexdata@cell),
                   local_n = hexdata@count)

  data <-
    data |>
    dplyr::mutate(hexcell = as.character(hexdata@cID)) |>
    dplyr::left_join(hexdata_cells_counts)
  return(data$local_n)

}

#' @importFrom SingleCellExperiment reducedDims
#' @importFrom tibble tibble
get_cds_umap_dims <- function(obj) {
  res <-
    tibble::tibble(
      cell_id = rownames(SingleCellExperiment::reducedDims(obj)$UMAP),
      UMAP_1 = SingleCellExperiment::reducedDims(obj)$UMAP[, 1],
      UMAP_2 = SingleCellExperiment::reducedDims(obj)$UMAP[, 2]
    )
  return(res)
}

#' @importFrom SeuratObject DefaultAssay
#' @importFrom tibble tibble
get_seurat_umap_dims <- function(obj, assay) {
  SeuratObject::DefaultAssay(obj) <- assay
  mat <- obj[["umap"]]@cell.embeddings
  res <- tibble::tibble(cell_id = rownames(mat),
                        UMAP_1 = mat[, 1],
                        UMAP_2 = mat[, 2])
  return(res)
}

`%notin%` <- Negate(`%in%`)

get_ica_space_df <- function(cds) {
  mat <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
  tibble::tibble(
    prin_graph_dim_1 = mat[, 1],
    prin_graph_dim_2 = mat[, 2],
    sample_name = rownames(mat),
    sample_state = rownames(mat)
  )
}

get_edge_df <- function(cds) {
  ica_space_df <- get_ica_space_df(cds)
  dp_mst <- cds@principal_graph[["UMAP"]]
  edge_df <- dp_mst %>%
    igraph::as_data_frame() %>%
    dplyr::select(source = "from", target = "to") %>%
    dplyr::left_join(
      ica_space_df %>%
        dplyr::select(
          source = "sample_name",
          source_prin_graph_dim_1 = "prin_graph_dim_1",
          source_prin_graph_dim_2 = "prin_graph_dim_2"
        ),
      by = "source"
    ) %>%
    dplyr::left_join(
      ica_space_df %>%
        dplyr::select(
          target = "sample_name",
          target_prin_graph_dim_1 = "prin_graph_dim_1",
          target_prin_graph_dim_2 = "prin_graph_dim_2"
        ),
      by = "target"
    )
  edge_df

}


make_trajectory_graph <-
  function(g,
           cds,
           trajectory_graph_segment_size,
           trajectory_graph_color) {
    edge_df <- get_edge_df(cds)
    g + geom_segment(
      aes_string(
        x = "source_prin_graph_dim_1",
        y = "source_prin_graph_dim_2",
        xend = "target_prin_graph_dim_1",
        yend = "target_prin_graph_dim_2"
      ),
      size = trajectory_graph_segment_size,
      color = I(trajectory_graph_color),
      linetype = "solid",
      na.rm = TRUE,
      data = edge_df
    )


  }

branch_nodes <- function(cds) {
  g <- principal_graph(cds)[["UMAP"]]
  branch_points <- which(igraph::degree(g) > 2)
  branch_points <-
    branch_points[branch_points %in% root_nodes(cds) == FALSE]
  return(branch_points)
}

root_nodes <- function(cds) {
  g <- principal_graph(cds)[["UMAP"]]
  root_pr_nodes <-
    which(names(igraph::V(g)) %in% cds@principal_graph_aux[["UMAP"]]$root_pr_nodes)
  names(root_pr_nodes) <-
    cds@principal_graph_aux[["UMAP"]]$root_pr_nodes
  return(root_pr_nodes)
}

leaf_nodes <- function (cds) {
  g <- principal_graph(cds)[["UMAP"]]
  leaves <- which(igraph::degree(g) == 1)
  leaves <- leaves[leaves %in% root_nodes(cds) == FALSE]
  return(leaves)
}

make_pp_labels <-
  function(g,
           cds,
           graph_label_size,
           trajectory_graph_segment_size) {
    ica_space_df <- get_ica_space_df(cds)
    mst_branch_nodes <- branch_nodes(cds)
    mst_leaf_nodes <- leaf_nodes(cds)
    mst_root_nodes <- root_nodes(cds)
    pps <- c(mst_branch_nodes, mst_leaf_nodes, mst_root_nodes)
    princ_point_df <- ica_space_df %>%
      dplyr::slice(match(names(pps), sample_name))

    g +
      geom_point(
        aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"),
        shape = 21,
        stroke = I(trajectory_graph_segment_size),
        color = "white",
        fill = "black",
        size = I(graph_label_size * 1.5),
        na.rm = TRUE,
        princ_point_df
      ) +
      ggrepel::geom_text_repel(
        aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2",
                   label = "sample_name"),
        size = I(graph_label_size * 1.5),
        color = "Black",
        na.rm = TRUE,
        princ_point_df
      )

  }










# make_branch_labels <- function(g,
#                                cds,
#                                graph_label_size,
#                                trajectory_graph_segment_size) {
#   ica_space_df <- get_ica_space_df(cds)
#   mst_branch_nodes <- branch_nodes(cds)
#   branch_point_df <- ica_space_df %>%
#     dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
#     dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
#
#   g +
#     geom_point(
#       aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"),
#       shape = 21,
#       stroke = I(trajectory_graph_segment_size),
#       color = "white",
#       fill = "black",
#       size = I(graph_label_size * 1.5),
#       na.rm = TRUE,
#       branch_point_df
#     ) +
#     geom_text(
#       aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2",
#                  label = "branch_point_idx"),
#       size = I(graph_label_size),
#       color = "white",
#       na.rm = TRUE,
#       branch_point_df
#     )
#
# }
# make_leaf_labels <- function(g,
#                              cds,
#                              trajectory_graph_segment_size,
#                              graph_label_size) {
#   mst_leaf_nodes <- leaf_nodes(cds)
#   ica_space_df <- get_ica_space_df(cds)
#   leaf_df <- ica_space_df %>%
#     dplyr::slice(match(names(mst_leaf_nodes), sample_name)) %>%
#     dplyr::mutate(leaf_idx = seq_len(dplyr::n()))
#
#   g +
#     geom_point(
#       aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"),
#       shape = 21,
#       stroke = I(trajectory_graph_segment_size),
#       color = "black",
#       fill = "lightgray",
#       size = I(graph_label_size * 1.5),
#       na.rm = TRUE,
#       leaf_df
#     ) +
#     geom_text(
#       aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2",
#                  label = "leaf_idx"),
#       size = I(graph_label_size),
#       color = "black",
#       na.rm = TRUE,
#       leaf_df
#     )
#
# }

# make_root_labels <- function(g,
#                              cds,
#                              trajectory_graph_segment_size,
#                              graph_label_size) {
#   mst_root_nodes <- root_nodes(cds)
#   ica_space_df <- get_ica_space_df(cds)
#   root_df <- ica_space_df %>%
#     dplyr::slice(match(names(mst_root_nodes), sample_name)) %>%
#     dplyr::mutate(root_idx = seq_len(dplyr::n()))
#
#   g +
#     geom_point(
#       aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"),
#       shape = 21,
#       stroke = I(trajectory_graph_segment_size),
#       color = "black",
#       fill = "white",
#       size = I(graph_label_size * 1.5),
#       na.rm = TRUE,
#       root_df
#     ) +
#     geom_text(
#       aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2",
#                  label = "root_idx"),
#       size = I(graph_label_size),
#       color = "black",
#       na.rm = TRUE,
#       root_df
#     )
# }

make_root_node_labels <- function(g,
                                  cds,
                                  pseudotime_dim,
                                  trajectory_graph_segment_size,
                                  graph_label_size) {
  zero_cell <- bb_cellmeta(cds) |>
    dplyr::filter(!!sym(pseudotime_dim) == 0) |>
    dplyr::pull(cell_id)
  if (length(zero_cell) > 1)
    cli::cli_abort("Couldn't identify a single cell to label as root node.")
  closest_pp <-
    cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex[zero_cell, ] |> unname()
  closest_pp <- paste0("Y_", closest_pp)
  ica_space_df <- get_ica_space_df(cds)
  root_df <- ica_space_df |>
    dplyr::filter(sample_name == closest_pp) |>
    dplyr::mutate(root_idx = "1")

  g +
    geom_point(
      aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"),
      shape = 21,
      stroke = I(trajectory_graph_segment_size),
      color = "black",
      fill = "white",
      size = I(graph_label_size * 1.5),
      na.rm = TRUE,
      root_df
    ) +
    geom_text(
      aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2",
                 label = "root_idx"),
      size = I(graph_label_size),
      color = "black",
      na.rm = TRUE,
      root_df
    )
}
