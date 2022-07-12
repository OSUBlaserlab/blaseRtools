# geom split violin

GeomSplitViolin <- ggplot2::ggproto(
  "GeomSplitViolin",
  ggplot2::GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <-
      transform(
        data,
        xminv = x - violinwidth * (x - xmin),
        xmaxv = x + violinwidth * (xmax - x)
      )
    grp <- data[1, "group"]
    newdata <-
      plyr::arrange(transform(data, x = if (grp %% 2 == 1)
        xminv
        else
          xmaxv), if (grp %% 2 == 1)
            y
        else-y)
    newdata <-
      rbind(newdata[1,], newdata, newdata[nrow(newdata),], newdata[1,])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <-
      round(newdata[1, "x"])

    if (length(draw_quantiles) > 0 &
        !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                1))
      quantiles <-
        ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <-
        data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <-
        rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <-
        ggplot2::GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin",
                       grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    }
    else {
      ggplot2:::ggname("geom_split_violin",
                       ggplot2::GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

#' @title geom_split_violin
#' @description create a split violin geom
#' @param mapping Set of aesthetitc mappings, Default: NULL
#' @param data Data to display, Default: NULL
#' @param stat Statistic to plot, Default: 'ydensity'
#' @param position Position, Default: 'identity'
#' @param ...
#' @param draw_quantiles Quantiles to draw, Default: NULL
#' @param trim Trim ends?, Default: TRUE
#' @param scale Analogous to violin scale, Default: 'area'
#' @param na.rm Default: FALSE
#' @param show.legend Default: NA
#' @param inherit.aes Default: TRUE
#' @rdname geom_split_violin
#' @export
geom_split_violin <-
  function(mapping = NULL,
           data = NULL,
           stat = "ydensity",
           position = "identity",
           ...,
           draw_quantiles = NULL,
           trim = TRUE,
           scale = "area",
           na.rm = FALSE,
           show.legend = NA,
           inherit.aes = TRUE) {
    layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomSplitViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        trim = trim,
        scale = scale,
        draw_quantiles = draw_quantiles,
        na.rm = na.rm,
        ...
      )
    )
  }
