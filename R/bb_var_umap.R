#' A function to generate a UMAP with colors mapped to colData variables
#'
#' @param cds A cell data set object
#' @param var The variable to map colors to.  Special exceptions are "density", "local_n" and "log_local_n" which calculate the 2 d kernel density estimate or binned cell counts and maps to color scale.
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
#' @param man_text_df A data frame in the form of text_x = <numeric vector>, text_y = <numeric vector>, label = <character vector> for manually placing text labels.
#' @return A ggplot
#' @export
#' @import tidyverse monocle3 ggrepel MASS hexbin
bb_var_umap <- function(cds,
                        var,
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
                        ...,
                        man_text_df = NULL) {
  plot_data <- plot_cells(cds)[["data"]]
  
  if (sample_equally) {
    if (is.null(facet_by)) {
      return("You need to provide a faceting variable.")
    } else {
      if (length(facet_by) == 1) {
        sample_sizes <-
          plot_data %>%
          group_by(!!sym(facet_by)) %>%
          summarise(ncells = n())
        plot_data <-
          plot_data %>%
          group_by(!!sym(facet_by)) %>%
          slice_sample(n = min(sample_sizes$ncells)) %>%
          ungroup()
        
        
      } else if (length(facet_by) == 2) {
        # create a composite variable to group by
        plot_data$cross_product <-
          paste0(plot_data[, facet_by[1]], "_", plot_data[, facet_by[2]])
        sample_sizes <-
          plot_data %>%
          group_by(cross_product) %>%
          summarise(ncells = n())
        plot_data <-
          plot_data %>%
          group_by(cross_product) %>%
          slice_sample(n = min(sample_sizes$ncells)) %>%
          ungroup()
      } else if (length(facet_by) > 2) {
        return("You can only facet by 2 variables")
      }
    }
  }
  
  if (var %in% c("density", "local_n", "log_local_n")) {
    data_long <- plot_data
  } else {
    data_long <-
      plot_data %>% pivot_longer(cols = (!!sym(var)), names_to = "var")
  }
  
  dim_x <- ifelse(is.null(alt_dim_x), "data_dim_1", alt_dim_x)
  dim_y <- ifelse(is.null(alt_dim_y), "data_dim_2", alt_dim_y)
  
  # generate text data frame for variable labels if you are going to use them
  if (is.null(alt_label_col)) {
    if (var %in% c("density", "local_n", "log_local_n")) {
      text_df <- data_long
    } else {
      text_df <- data_long %>% group_by(value)
    }
    
  } else {
    text_df <-
      plot_data %>% pivot_longer(cols = !!sym(alt_label_col),
                                 names_to = "var") %>% group_by(value)
  }
  if (overwrite_labels == T && is.null(man_text_df)) {
    median_coord_df <-
      text_df %>%
      summarise(
        fraction_of_group = n(),
        text_x = median(!!sym(dim_x)),
        text_y = median(!!sym(dim_y))
      )
    text_df <- left_join(text_df, median_coord_df) %>%
      mutate(label = value)
    text_df <-
      text_df %>% group_by(label, text_x, text_y) %>% summarise()
  }
  
  if (!is.null(man_text_df))
    text_df <- man_text_df
  
  # make the main plot
  plot <- ggplot()
  if (!is.null(value_to_highlight)) {
    data_background <-
      data_long %>% filter(value %notin% value_to_highlight)
    data_long <- data_long %>% filter(value %in% value_to_highlight)
    plot <- plot +
      geom_point(
        data = data_background,
        aes(x = !!sym(dim_x),
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
          map_dfr(
            .x = data_long %>% group_by(!!sym(facet_by)) %>% summarise() %>% pull() %>% as.character(),
            .f = function(x, density_data = data_long) {
              density_data <-
                density_data %>% filter(!!sym(facet_by) == x) %>% as.data.frame()
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
          map_dfr(
            .x = data_long$cross_product %>% unique(),
            .f = function(x, density_data = data_long) {
              density_data <-
                density_data %>% filter(cross_product == x) %>% as.data.frame()
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
        data_long %>% mutate(density = get_density(x = data[, dim_x], y = data[, dim_y], n = nbin))
    }
    
    plot <- ggplot(data_long) +
      geom_point(
        aes_string(
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
      scale_color_viridis_c(option = "inferno",
                            begin = 0.1,
                            end = 0.9) +
      scale_fill_viridis_c(
        option = "inferno",
        begin = 0.1,
        end = 0.9,
        guide = F
      )
  } else if (var %in% c("local_n", "log_local_n")) {
    if (!is.null(facet_by)) {
      if (length(facet_by) == 1) {
        data_long <-
          map_dfr(
            .x = data_long %>% pull(!!sym(facet_by)) %>% unique(),
            .f = function(x, data = data_long) {
              data <- data %>% filter(!!sym(facet_by) == x)
              data <-
                data %>% mutate(local_n = get_hexcount(
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
          map_dfr(
            .x = data_long$cross_product %>% unique(),
            .f = function(x, data = data_long) {
              data <- data %>% filter(cross_product == x)
              data <-
                data %>% mutate(local_n = get_hexcount(
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
        data_long %>% mutate(local_n = get_hexcount(
          data = data_long,
          x = dim_x,
          y = dim_y,
          n = nbin
        ))
    }
    data_long <- data_long %>% mutate(log_local_n = log10(local_n))
    plot <- ggplot(data_long) +
      geom_point(
        aes_string(
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
      scale_color_viridis_c(option = "inferno",
                            begin = 0.1,
                            end = 0.9) +
      scale_fill_viridis_c(
        option = "inferno",
        begin = 0.1,
        end = 0.9,
        guide = F
      )
  } else {
    plot <- plot +
      geom_point(
        data = data_long,
        aes(
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
        scale_fill_viridis_c(guide = "colorbar", na.value = "transparent") +
        scale_color_viridis_c(guide = F, na.value = "grey80")
    } else if (length(palette) == 1 && palette == "viridis") {
      plot <- plot +
        scale_fill_viridis_d(begin = 0.1, end = 0.9) +
        scale_color_viridis_d(begin = 0.1,
                              end = 0.9,
                              guide = F)
    } else if (length(palette) == 1 && palette == "rcolorbrewer") {
      plot <- plot + scale_color_brewer(palette = "Paired", guide = F) +
        scale_fill_brewer(palette = "Paired")
    } else if (!is.null(palette)) {
      plot <- plot + scale_color_manual(values = palette, guide = F) +
        scale_fill_manual(values = palette)
    } else {
      plot <- plot + scale_color_discrete(guide = F) +
        scale_fill_discrete()
    }
    if (class(data_long$value) != "numeric") {
      plot <-
        plot + guides(fill = guide_legend(override.aes = list(
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
    theme(plot.title = element_text(hjust = 0.5))
  plot <- plot +
    theme(legend.position = legend_pos)#+coord_fixed()
  
  #option to overwrite labels
  if (overwrite_labels == T && is.null(man_text_df)) {
    plot <- plot +
      theme(legend.position = "none") +
      ggrepel::geom_text_repel(
        data = text_df,
        mapping = aes_string(x = "text_x", y = "text_y", label = "label"),
        size = group_label_size,
        min.segment.length = 1
      )
  } else if (!is.null(man_text_df)) {
    plot <- plot +
      theme(legend.position = "none") +
      geom_text(
        data = text_df,
        mapping = aes_string(x = "text_x", y = "text_y", label = "label"),
        size = group_label_size
      )
  }
  if (!is.null(facet_by)) {
    if (length(facet_by) == 1) {
      plot <- plot +
        facet_wrap(facets = facet_by, ...) +
        theme(strip.background = element_blank())
    } else if (length(facet_by) == 2) {
      plot <- plot +
        facet_grid(...) +
        theme(strip.background = element_blank())
    } else {
      return("Too many dimensions to facet by.")
    }
  }
  
  return(plot)
}

get_density <- function(xgd, ygd, ngd) {
  dens <- MASS::kde2d(x = xgd, y = ygd, n = ngd)
  ix <- findInterval(xgd, dens$x)
  iy <- findInterval(ygd, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

get_hexcount <- function(data, x, y, n) {
  hexdata <- hexbin::hexbin(
    x = data %>% pull(!!sym(x)),
    y = data %>% pull(!!sym(y)),
    IDs = T,
    xbins = n
  )
  hexdata_cells_counts <-
    tibble(hexcell = as.character(hexdata@cell),
           local_n = hexdata@count)
  
  data <-
    data %>%
    mutate(hexcell = as.character(hexdata@cID)) %>%
    left_join(hexdata_cells_counts)
  return(data$local_n)
  
}
