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
#' @param ... Additional params for facetting.
#' @param man_text_df A data frame in the form of text_x = numeric_vector, text_y = numeric_vector, label = character_vector for manually placing text labels.
#' @param cds Provided for backward compatibility with prior versions.  If a value is supplied, a warning will be emitted and the value will be transferred to the obj argument, Default: NULL
#' @return a ggplot
#' @rdname bb_var_umap
#' @export
#' @importFrom dplyr left_join group_by summarise slice_sample ungroup n mutate filter pull
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_point aes aes_string scale_color_viridis_c scale_fill_viridis_c scale_fill_viridis_d scale_color_viridis_d scale_color_brewer scale_fill_brewer scale_color_manual scale_fill_manual scale_color_discrete scale_fill_discrete guides guide_legend theme element_text geom_text facet_wrap element_blank facet_grid
#' @importFrom purrr map_dfr
#' @importFrom ggrepel geom_text_repel
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
                        cds = NULL,
                        ...,
                        man_text_df = NULL) {


  cds_warn(cds)
  obj_stop(obj)
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
    data_long <- data_long |> dplyr::filter(value %in% value_to_highlight)
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
    data_long <- data_long |> dplyr::mutate(log_local_n = log10(local_n))
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
      plot <- plot + ggplot2::scale_color_brewer(palette = "Paired", guide = "none") +
        ggplot2::scale_fill_brewer(palette = "Paired")
    } else if (!is.null(palette)) {
      plot <- plot + ggplot2::scale_color_manual(values = palette, guide = "none") +
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
      ggrepel::geom_text_repel(
        data = text_df,
        mapping = ggplot2::aes_string(x = "text_x", y = "text_y", label = "label"),
        size = group_label_size,
        min.segment.length = 1
      )
  } else if (!is.null(man_text_df)) {
    plot <- plot +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_text(
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
  res <- tibble::tibble(cell_id = rownames(SingleCellExperiment::reducedDims(obj)$UMAP),
                        UMAP_1 = SingleCellExperiment::reducedDims(obj)$UMAP[,1],
                        UMAP_2 = SingleCellExperiment::reducedDims(obj)$UMAP[,2])
  return(res)
}

#' @importFrom SeuratObject DefaultAssay
#' @importFrom tibble tibble
get_seurat_umap_dims <- function(obj, assay) {
  SeuratObject::DefaultAssay(obj) <- assay
  mat <- obj[["umap"]]@cell.embeddings
  res <- tibble::tibble(cell_id = rownames(mat),
                        UMAP_1 = mat[,1],
                        UMAP_2 = mat[,2])
  return(res)
}

`%notin%` <- negate(`%in%`)
