
#' @title Plot a Heatmap Row Dendrogram
#' @description Takes in a SummarizedHeatmap object and returns a ggplot of the rowDendro slot.  This can be positioned with the side parameter.  Default is to position it on the left.  If flipping the heatmap so that the rows run vertically, you will need to change the side argument to top or bottom.
#' @param obj a Summarized Heatmap
#' @param side Orientation/side of the heatmap to plot, Default: c("left", "right", "top", "bottom")
#' @param linewidth Weight of the lines for the dendrogram, Default: 0.5
#' @return A ggplot
#' @rdname bb_plot_heatmap_rowDendro
#' @export
#' @importFrom cli cli_abort
#' @importFrom ggdendro dendro_data segment
#' @import ggplot2
#' @import patchwork
bb_plot_heatmap_rowDendro <-
  function(obj,
           side = c("left", "right", "top", "bottom"),
           linewidth = 0.5) {
    if (!"SummarizedHeatmap" %in% class(obj)) {
      cli::cli_abort("This function must be run on a SummarizedHeatmap object")
    }
    side <- match.arg(side)
    dg <- rowDendro(obj)
    ddata <- ggdendro::dendro_data(dg, type = "rectangle")

    if (side %in% c("left", "right")) {
      dend <- ggplot(ggdendro::segment(ddata)) +
        geom_segment(aes(
          x = x,
          y = y,
          xend = xend,
          yend = yend
        ), linewidth = linewidth) +
        coord_flip()
      if (side == "left") {
        dend <- dend +
          scale_y_reverse(expand = expansion(add = c(0.1, 0))) +
          scale_x_continuous(expand = expansion(add = c(0.5, 0.5))) +
          theme_dendro()
      } else {
        dend <- dend +
          scale_x_continuous(expand = expansion(add = c(0.5, 0.5))) +
          scale_y_continuous(expand = expansion(add = c(0.1, 0.1))) +
          theme_dendro()
      }
    } else {
      dend <- ggplot(segment(ddata)) +
        geom_segment(aes(
          x = y,
          y = x,
          xend = yend,
          yend = xend
        ), linewidth = linewidth) +
        coord_flip()
      if (side == "bottom") {
        dend <- dend +
          scale_x_reverse(expand = expansion(add = c(0.1, 0.1))) +
          scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
          theme_dendro()
      } else {
        dend <- dend +
          scale_x_continuous(expand = expansion(add = c(0.1, 0.1))) +
          scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
          theme_dendro()

      }


    }

    dend
  }

#' @title Plot a Heatmap Column Dendrogram
#' @description Takes in a SummarizedHeatmap object and returns a ggplot of the rowDendro slot.  This can be positioned with the side parameter.  Default is to position it on the left.  If flipping the heatmap so that the rows run vertically, you will need to change the side argument to top or bottom.
#' @param obj a Summarized Heatmap
#' @param side Orientation/side of the heatmap to put the dendrogram, Default: c("top", "bottom", "left", "right")
#' @param linewidth Weight of the dendrogram plot, Default: 0.5
#' @return a ggplot
#' @rdname bb_plot_heatmap_colDendro
#' @export
#' @importFrom cli cli_abort
#' @importFrom ggdendro dendro_data segment
#' @import ggplot2
#' @import patchwork
bb_plot_heatmap_colDendro <-
  function(obj,
           side = c("top", "bottom", "left", "right"),
           linewidth = 0.5) {
    if (!"SummarizedHeatmap" %in% class(obj)) {
      cli::cli_abort("This function must be run on a SummarizedHeatmap object")
    }
    side <- match.arg(side)
    dg <- colDendro(obj)
    ddata <- ggdendro::dendro_data(dg, type = "rectangle")

    if (side %in% c("top", "bottom")) {
      dend <- ggplot(segment(ddata)) +
        geom_segment(aes(
          x = x,
          y = y,
          xend = xend,
          yend = yend
        ), linewidth = linewidth)
      if (side == "bottom") {
        dend <- dend +
          scale_y_reverse(add = c(0.1, 0.1)) +
          scale_x_continuous(expand = expansion(add = c(0.5, 0.5))) +
          theme_dendro()
      } else {
        dend <- dend +
          scale_y_continuous(expand = expansion(add = c(0.1, 0.1))) +
          scale_x_continuous(expand = expansion(add = c(0.5, 0.5))) +
          theme_dendro()
      }


    } else {
      dend <- ggplot(ggdendro::segment(ddata)) +
        geom_segment(aes(
          x = y,
          y = x,
          xend = yend,
          yend = xend
        ), linewidth = linewidth)
      if (side == "left") {
        dend <- dend +
          scale_x_reverse(expand = expansion(add = c(0.1, 0.1))) +
          scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
          theme_dendro()
      } else {
        dend <- dend +
          scale_x_continuous(expand = expansion(add = c(0.1, 0.1))) +
          scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
          theme_dendro()

      }


    }
    dend
  }

#' @title Plot the Body of  Heatmap
#' @description Takes in a Summarized Heatmap object and returns a ggplot of the matrix data.  Currently only supports clustering.
#' @param obj A SummarizedHeatmap
#' @param tile_color Outline of the color tiles, Default: 'white'
#' @param high Color for high values, applied to scale_fill_gradient_2, Default: 'red3'
#' @param mid Color for mid values, applied to scale_fill_gradient_2, Default: 'white'
#' @param low Color for low values, applied to scale_fill_gradient_2, Default: 'blue4'
#' @param flip Whether to transpose the matrix, i.e. plot the rows as columns and columns as rows, Default: FALSE
#' @return a ggplot
#' @rdname bb_plot_heatmap_main
#' @export
#' @importFrom cli cli_abort
#' @importFrom ggdendro dendro_data
#' @importFrom SummarizedExperiment assay
#' @import ggplot2
#' @import patchwork
bb_plot_heatmap_main <-
  function(obj,
           tile_color = "white",
           high = "red3",
           mid = "white",
           low = "blue4",
           flip = FALSE) {
    if (!"SummarizedHeatmap" %in% class(obj)) {
      cli::cli_abort("This function must be run on a SummarizedHeatmap object")
    }
    dg <- rowDendro(obj)
    ddata <- ggdendro::dendro_data(dg, type = "rectangle")
    dh <- colDendro(obj)
    hdata <- ggdendro::dendro_data(dh, type = "rectangle")
    hm_dat <- SummarizedExperiment::assay(obj) |>
      as_tibble(rownames = "rownames") |>
      pivot_longer(-rownames) |>
      mutate(rownames = factor(rownames, levels = ddata$labels$label)) |>
      mutate(name = factor(name, levels = hdata$labels$label))

    if (!flip) {
      hm <- ggplot(hm_dat, aes(x = name, y = rownames, fill = value)) +
        geom_tile(color = tile_color) +
        scale_fill_gradient2(high = high,
                             mid = mid,
                             low = low) +
        scale_y_discrete(position = "right", expand = expansion(add = 0)) +
        scale_x_discrete(expand = expansion(add = 0)) +
        labs(x = NULL, y = NULL, fill = "Expression") +
        panel_border(color = "black") +
        theme_hm_main()

    } else {
      hm <- ggplot(hm_dat, aes(y = name, x = rownames, fill = value)) +
        geom_tile(color = tile_color) +
        scale_fill_gradient2(high = high,
                             mid = mid,
                             low = low) +
        scale_y_discrete(position = "right", expand = expansion(add = 0)) +
        scale_x_discrete(expand = expansion(add = 0)) +
        labs(x = NULL, y = NULL, fill = "Expression") +
        panel_border(color = "black") +
        theme_hm_main()

    }

    hm

  }

#' @title Plot a Column Highligh
#' @description Use geom_text_repel to selectively highlight some column names.  Useful when there are too many to highlight to be able to use the axis directly.
#' @param obj A summarized heatmap
#' @param highlights A vector of columns to highlight, Default: character(0)
#' @param side Side on which to put the higlight, Default: c("top", "bottom", "right", "left")
#' @param ... Other arguments to pass to geom_text_repel
#' @return a ggplot
#' @rdname bb_plot_heatmap_colHighlight
#' @export
#' @importFrom cli cli_abort
#' @importFrom ggdendro dendro_data
#' @import ggplot2
#' @import patchwork
#' @import ggrepel
bb_plot_heatmap_colHighlight <-
  function(obj,
           highlights = character(0),
           side = c("top", "bottom", "right", "left"),
           ...) {
    if (!"SummarizedHeatmap" %in% class(obj)) {
      cli::cli_abort("This function must be run on a SummarizedHeatmap object")
    }
    side <- match.arg(side)
    dg <- rowDendro(obj)
    ddata <- ggdendro::dendro_data(dg, type = "rectangle")
    dh <- colDendro(obj)
    hdata <- ggdendro::dendro_data(dh, type = "rectangle")

    hm_dat <- assay(obj) |>
      as_tibble(rownames = "rownames") |>
      pivot_longer(-rownames) |>
      mutate(rownames = factor(rownames, levels = ddata$labels$label)) |>
      mutate(name = factor(name, levels = hdata$labels$label))


    if (side %in% c("top", "bottom")) {
      label_dat <- hm_dat |>
        group_by(name) |>
        summarize() |>
        mutate(y = 0) |>
        mutate(label = ifelse(name %in% highlights, as.character(name), NA))
      th <- ggplot(label_dat, aes(x = name, y = y, label = label))

      if (side == "bottom") {
        th <- th +
          geom_text_repel(
            ...,
            angle = 90,
            force        = 1,
            box.padding = 0.25,
            nudge_y      = -0.01,
            direction    = "x",
            hjust        = 1,
            segment.size = 0.5,
            min.segment.length = 0,
            segment.curvature = -0.01,
            aes(
              segment.square  = TRUE,
              segment.inflect = TRUE,
            )
          ) +
          labs(x = NULL, y = NULL) +
          theme(plot.margin = margin(0, 0, 0, 0)) +
          scale_y_continuous(limits = c(NA, 0)) +
          theme_nothing(font_size = theme_get()$text$size)

      } else {
        th <- th +
          geom_text_repel(
            ...,
            angle = 90,
            force        = 1,
            box.padding = 0.25,
            nudge_y      = 0.01,
            direction    = "x",
            hjust        = 1,
            segment.size = 0.5,
            min.segment.length = 0,
            segment.curvature = -0.01,
            aes(
              segment.square  = TRUE,
              segment.inflect = TRUE,
            )
          ) +

          labs(x = NULL, y = NULL) +
          theme(plot.margin = margin(0, 0, 0, 0)) +
          scale_y_reverse(limits = c(0, NA)) +
          theme_nothing(font_size = theme_get()$text$size)
      }


    } else {
      label_dat <- hm_dat |>
        group_by(name) |>
        summarize() |>
        mutate(y = 0) |>
        mutate(label = ifelse(name %in% highlights, as.character(name), NA))
      th <- ggplot(label_dat, aes(y = name, x = y, label = label))
      if (side == "left") {
        th <- th +
          geom_text_repel(
            ...,
            angle = 0,
            force        = 1,
            box.padding = 0.25,
            nudge_x      = -0.01,
            direction    = "y",
            hjust        = 1,
            segment.size = 0.5,
            min.segment.length = 0,
            segment.curvature = -0.01,
            aes(
              segment.square  = TRUE,
              segment.inflect = TRUE,
            )
          ) +
          labs(x = NULL, y = NULL) +
          theme(plot.margin = margin(0, 0, 0, 0)) +
          scale_x_continuous(limits = c(NA, 0)) +
          theme_nothing(font_size = theme_get()$text$size)

      } else {
        th <- th +
          geom_text_repel(
            ...,
            angle = 0,
            force        = 1,
            box.padding = 0.25,
            nudge_x      = 0.01,
            direction    = "y",
            hjust        = 1,
            segment.size = 0.5,
            min.segment.length = 0,
            segment.curvature = -0.01,
            aes(
              segment.square  = TRUE,
              segment.inflect = TRUE,
            )
          ) +
          labs(x = NULL, y = NULL) +
          theme(plot.margin = margin(0, 0, 0, 0)) +
          scale_x_reverse(limits = c(0, NA)) +
          theme_nothing(font_size = theme_get()$text$size)
      }

    }

    th

  }

#' @title Plot a Row Highlight
#' @description Plots a selected row name from a SummarizedHeatmap object.  Useful when there are too many to highlight so you can't alter the plot axis.
#' @param obj A summarized Heatmap object
#' @param highlights A vector of rows to highlight, Default: character(0)
#' @param side Side on which to plot the highlight, Default: c("top", "bottom", "right", "left")
#' @param ... Other arguments to pass to geom_text_repel
#' @return a ggplot
#' @rdname bb_plot_heatmap_rowHighlight
#' @export
#' @importFrom cli cli_abort
#' @importFrom ggdendro dendro_data
#' @import ggplot2
#' @import patchwork
#' @import ggrepel
bb_plot_heatmap_rowHighlight <-
  function(obj,
           highlights = character(0),
           side = c("top", "bottom", "right", "left"),
           ...) {
    if (!"SummarizedHeatmap" %in% class(obj)) {
      cli::cli_abort("This function must be run on a SummarizedHeatmap object")
    }
    side <- match.arg(side)
    dg <- rowDendro(obj)
    ddata <- ggdendro::dendro_data(dg, type = "rectangle")
    dh <- colDendro(obj)
    hdata <- ggdendro::dendro_data(dh, type = "rectangle")

    hm_dat <- assay(obj) |>
      as_tibble(rownames = "rownames") |>
      pivot_longer(-rownames) |>
      mutate(rownames = factor(rownames, levels = ddata$labels$label)) |>
      mutate(name = factor(name, levels = hdata$labels$label))


    if (side %in% c("top", "bottom")) {
      label_dat <- hm_dat |>
        group_by(rownames) |>
        summarize() |>
        mutate(y = 0) |>
        mutate(label = ifelse(rownames %in% highlights, as.character(rownames), NA))

      th <-
        ggplot(label_dat, aes(x = rownames, y = y, label = label))

      if (side == "bottom") {
        th <- th +
          geom_text_repel(
            ...,
            angle = 90,
            force        = 1,
            box.padding = 0.25,
            nudge_y      = -0.01,
            direction    = "x",
            hjust        = 1,
            segment.size = 0.5,
            min.segment.length = 0,
            segment.curvature = -0.01,
            aes(
              segment.square  = TRUE,
              segment.inflect = TRUE,
            )
          ) +
          labs(x = NULL, y = NULL) +
          theme(plot.margin = margin(0, 0, 0, 0)) +
          scale_y_continuous(limits = c(NA, 0)) +
          theme_nothing(font_size = theme_get()$text$size)

      } else {
        th <- th +
          geom_text_repel(
            ...,
            angle = 90,
            force        = 1,
            box.padding = 0.25,
            nudge_y      = 0.01,
            direction    = "x",
            hjust        = 1,
            segment.size = 0.5,
            min.segment.length = 0,
            segment.curvature = -0.01,
            aes(
              segment.square  = TRUE,
              segment.inflect = TRUE,
            )
          ) +

          labs(x = NULL, y = NULL) +
          theme(plot.margin = margin(0, 0, 0, 0)) +
          scale_y_reverse(limits = c(0, NA)) +
          theme_nothing(font_size = theme_get()$text$size)
      }


    } else {
      label_dat <- hm_dat |>
        group_by(rownames) |>
        summarize() |>
        mutate(y = 0) |>
        mutate(label = ifelse(rownames %in% highlights, as.character(rownames), NA))
      th <-
        ggplot(label_dat, aes(y = rownames, x = y, label = label))
      if (side == "left") {
        th <- th +
          geom_text_repel(
            ...,
            angle = 0,
            force        = 1,
            box.padding = 0.25,
            nudge_x      = -0.01,
            direction    = "y",
            hjust        = 1,
            segment.size = 0.5,
            min.segment.length = 0,
            segment.curvature = -0.01,
            aes(
              segment.square  = TRUE,
              segment.inflect = TRUE,
            )
          ) +
          labs(x = NULL, y = NULL) +
          theme(plot.margin = margin(0, 0, 0, 0)) +
          scale_x_continuous(limits = c(NA, 0)) +
          theme_nothing(font_size = theme_get()$text$size)

      } else {
        th <- th +
          geom_text_repel(
            ...,
            angle = 0,
            force        = 1,
            box.padding = 0.25,
            nudge_x      = 0.01,
            direction    = "y",
            hjust        = 1,
            segment.size = 0.5,
            min.segment.length = 0,
            segment.curvature = -0.01,
            aes(
              segment.square  = TRUE,
              segment.inflect = TRUE,
            )
          ) +
          labs(x = NULL, y = NULL) +
          theme(plot.margin = margin(0, 0, 0, 0)) +
          scale_x_reverse(limits = c(0, NA)) +
          theme_nothing(font_size = theme_get()$text$size)
      }

    }

    th

  }

#' @title Plot SummarizedHeatmap colData
#' @description Will generate a ggplot from the colData of a SummarizedHeatmap object.  Typically this will be placed at the top or bottom of a plot.  If flipped, use the side argument to put the colData on the left or right.
#' @param obj a SummarizedHeatmap object
#' @param tile_color Outline color for the tiles, Default: 'white'
#' @param vars Variables to plot.  Supply a named vector to change the axis text and legend titles for each variable, Default: colnames(colData(obj))
#' @param side Side on which to plot, Default: c("top", "bottom", "right", "left")
#' @param manual_pal a color palette, preferably a named vector corresponding to values of colData, Default: NULL
#' @return a ggplot
#' @rdname bb_plot_heatmap_colData
#' @export
#' @importFrom cli cli_abort cli_alert_info
#' @importFrom ggdendro dendro_data
#' @importFrom dplyr filter
#' @importFrom patchwork wrap_plots
#' @import ggplot2
#' @import patchwork
bb_plot_heatmap_colData <-
  function(obj,
           tile_color = "white",
           vars = colnames(colData(obj)),
           side = c("top", "bottom", "right", "left"),
           manual_pal = NULL) {
    if (!"SummarizedHeatmap" %in% class(obj)) {
      cli::cli_abort("This function must be run on a SummarizedHeatmap object")
    }
    side <- match.arg(side)
    dh <- colDendro(obj)
    hdata <- ggdendro::dendro_data(dh, type = "rectangle")

    dat <- colData(obj) |>
      as_tibble(rownames = "colData_rownames") |>
      pivot_longer(-colData_rownames) |>
      dplyr::filter(name %in% vars) |>
      mutate(colData_rownames = factor(colData_rownames, levels = hdata$labels$label))

    if (is.null(names(vars))) {
      cli_div(theme = list(span.emph = list(color = "orange")))
      cli::cli_alert_info(
        "To change the appearance of the annotation labels and legend titles,\n supply a named vector to the variable {.emph vars}."
      )
      names(vars) <- vars
    }

    if (side %in% c("top", "bottom")) {
      plotlist <- map2(.x = vars,
                       .y = names(vars),
                       .f = \(x, y, d = dat, pal = manual_pal) {
                         d <- dplyr::filter(d, name == x) |>
                           mutate(label = y)
                         p <- ggplot(data = d,
                                     aes(x = colData_rownames,
                                         fill = value,
                                         y = label)) +
                           geom_tile(color = tile_color) +
                           theme_hm_colData_top() +
                           scale_y_discrete(position = "right", expand = expansion(add = 0)) +
                           scale_x_discrete(expand = expansion(add = 0)) +
                           labs(fill = y, x = NULL, y = NULL)
                         if (!is.null(pal))
                           p <- p + scale_fill_manual(values = pal)
                         p

                       })

      p <-
        patchwork::wrap_plots(plotlist, ncol = 1, axes = "collect") &
        theme(plot.margin = margin(0, 0, 0, 0))
      p

    } else {
      plotlist <- map2(.x = vars,
                       .y = names(vars),
                       .f = \(x, y, d = dat, pal = manual_pal) {
                         d <- dplyr::filter(d, name == x) |>
                           mutate(label = y)
                         p <- ggplot(data = d,
                                     aes(y = colData_rownames,
                                         fill = value,
                                         x = label)) +
                           geom_tile(color = tile_color) +
                           theme_hm_colData_left() +
                           scale_x_discrete(position = "top",
                                            expand = expansion(add = 0)) +
                           scale_y_discrete(expand = expansion(add = 0)) +
                           labs(fill = y, x = NULL, y = NULL)
                         if (!is.null(pal))
                           p <- p + scale_fill_manual(values = pal)
                         p

                       })

      p <-
        patchwork::wrap_plots(plotlist, nrow = 1, axes = "collect") &
        theme(plot.margin = margin(0, 0, 0, 0))
      # free(p, type = "space")
      p
    }
  }








#' @title Plot a SummarizedHeatmap rowData
#' @description Use this function to create a plot annotating SummarizedHeatmap rowDAta
#' @param obj A SummarizedHeatmap objectr
#' @param tile_color Color for the tile outlines, Default: 'white'
#' @param vars rowData variables to plot.  Supply a named vector to change the names shown on the axis and legend, Default: colnames(rowData(obj))
#' @param side Side on which to plotj, Default: c("right", "left", "top", "bottom")
#' @param manual_pal Color palette for filling the tiles, preferably a named vector, Default: NULL
#' @return a ggplot
#' @rdname bb_plot_heatmap_rowData
#' @export
#' @importFrom cli cli_abort cli_alert_info
#' @importFrom ggdendro dendro_data
#' @importFrom dplyr filter
#' @importFrom patchwork wrap_plots
#' @import ggplot2
#' @import patchwork
bb_plot_heatmap_rowData <-
  function(obj,
           tile_color = "white",
           vars = colnames(rowData(obj)),
           side = c("right", "left", "top", "bottom"),
           manual_pal = NULL) {
    if (!"SummarizedHeatmap" %in% class(obj)) {
      cli::cli_abort("This function must be run on a SummarizedHeatmap object")
    }
    side <- match.arg(side)
    dg <- rowDendro(obj)
    ddata <- ggdendro::dendro_data(dg, type = "rectangle")

    dat <- rowData(obj) |>
      as_tibble(rownames = "rowData_rownames") |>
      pivot_longer(-rowData_rownames) |>
      dplyr::filter(name %in% vars) |>
      mutate(rowData_rownames = factor(rowData_rownames, levels = ddata$labels$label))

    if (is.null(names(vars))) {
      cli_div(theme = list(span.emph = list(color = "orange")))
      cli::cli_alert_info(
        "To change the appearance of the annotation labels and legend titles,\n supply a named vector to the variable {.emph vars}."
      )
      names(vars) <- vars

    }

    if (side %in% c("left", "right")) {
      plotlist <- map2(.x = vars,
                       .y = names(vars),
                       .f = \(x, y, d = dat, pal = manual_pal) {
                         d <- dplyr::filter(d, name == x) |>
                           mutate(label = y)

                         p <- ggplot(data = d,
                                     aes(y = rowData_rownames, fill = value, x = label)) +
                           geom_tile(color = tile_color) +
                           scale_x_discrete(position = "top", expand =
                                              expansion(add = 0)) +
                           theme_hm_rowData_left() +
                           labs(fill = y, x = NULL, y = NULL)
                         if (!is.null(pal))
                           p <- p + scale_fill_manual(values = pal)
                         p

                       })

      p <-
        patchwork::wrap_plots(plotlist, nrow = 1, axes = "collect") &
        theme(plot.margin = margin(0, 0, 0, 0))
      # free(p, type = "space")
      p
    } else {
      plotlist <- map2(.x = vars,
                       .y = names(vars),
                       .f = \(x, y, d = dat, pal = manual_pal) {
                         d <- dplyr::filter(d, name == x) |>
                           mutate(label = y)

                         p <- ggplot(data = d,
                                     aes(x = rowData_rownames,
                                         fill = value,
                                         y = label)) +
                           geom_tile(color = tile_color) +
                           scale_y_discrete(position = "right",
                                            expand =
                                              expansion(add = 0)) +
                           theme_hm_rowData_top() +

                           labs(fill = y, x = NULL, y = NULL)
                         if (!is.null(pal))
                           p <- p + scale_fill_manual(values = pal)
                         p

                       })

      p <-
        patchwork::wrap_plots(plotlist, ncol = 1, axes = "collect") &
        theme(plot.margin = margin(0, 0, 0, 0))
      p
    }
  }
