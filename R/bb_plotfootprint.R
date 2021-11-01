#' Plot motif footprinting results
#'
#' @param object A Seurat object
#' @param features A vector of features to plot
#' @param alt_main_title Alternative title for the main plot.  Accepts markdown.
#' @param alt_color_title Alternative title for the color scale
#' @param legend_pos Position to place the legend
#' @param assay Name of assay to use
#' @param group.by A grouping variable
#' @param idents Set of identities to include in the plot
#' @param show.expected Plot the expected Tn5 integration frequency below the
#' main footprint plot
#' @param normalization Method to normalize for Tn5 DNA sequence bias. Options
#' are "subtract", "divide", or NULL to perform no bias correction.
#' @param label TRUE/FALSE value to control whether groups are labeled.
#' @param repel Repel labels from each other
#' @param label.top Number of groups to label based on highest accessibility
#' in motif flanking region.
#' @param label.idents Vector of identities to label. If supplied,
#' @param fontsize Theme font size
#' @param linesize Size to draw the footprint lines
#' \code{label.top} will be ignored.
#' @export
#' @concept visualization
#' @concept footprinting
#' @importFrom Seurat DefaultAssay
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap xlab ylab theme_classic
#' theme element_blank geom_label guides guide_legend
#' @importFrom dplyr group_by summarize top_n
#' @import cowplot Signac Seurat
bb_PlotFootprint <- function(object,
                            features,
                            alt_main_title = NULL,
                            alt_color_title = NULL,
                            legend_pos = "right",
                            assay = NULL,
                            group.by = NULL,
                            idents = NULL,
                            label = TRUE,
                            repel = TRUE,
                            show.expected = TRUE,
                            normalization = "subtract",
                            label.top = 3,
                            label.idents = NULL,
                            fontsize = 14,
                            linesize = 0.2) {
  assay <-
    Signac:::SetIfNull(x = assay, y = DefaultAssay(object = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay.")
  }
  plot.data <- Signac::GetFootprintData(
    object = object,
    features = features,
    assay = assay,
    group.by = group.by,
    idents = idents
  )
  motif.sizes <- Signac:::GetMotifSize(object = object,
                                       features = features,
                                       assay = assay)
  obs <- plot.data[plot.data$class == "Observed",]
  expect <- plot.data[plot.data$class == "Expected",]

  # flanks are motif edge to 50 bp each side
  # add flank information (T/F)
  base <- ceiling(motif.sizes / 2)
  obs$flanks <- sapply(
    X = seq_len(length.out = nrow(x = obs)),
    FUN = function(x) {
      pos <- abs(obs[x, "position"])
      size <- base[[obs[x, "feature"]]]
      return((pos > size) & (pos < (size + 50)))
    }
  )

  if (!is.null(normalization)) {
    # need to group by position and motif
    correction.vec <- expect$norm.value
    names(correction.vec) <- paste(expect$position, expect$feature)
    if (normalization == "subtract") {
      obs$norm.value <- obs$norm.value - correction.vec[paste(obs$position, obs$feature)]
    } else if (normalization == "divide") {
      obs$norm.value <- obs$norm.value / correction.vec[paste(obs$position, obs$feature)]
    } else {
      stop("Unknown normalization method requested")
    }
  }

  # find flanking accessibility for each group and each feature
  flanks <- obs[obs$flanks,]
  flanks <- group_by(.data = flanks, feature, group)
  flankmeans <- summarize(.data = flanks, mn = mean(x = norm.value))

  # find top n groups for each feature
  topmean <- top_n(x = flankmeans, n = label.top, wt = mn)

  # find the top for each feature to determine axis limits
  ymax <- top_n(x = flankmeans, n = 1, wt = mn)
  ymin <- top_n(x = flankmeans, n = 1, wt = -mn)

  # make df for labels
  label.df <- data.frame()
  sub <- obs[obs$position == 75,]
  for (i in seq_along(along.with = features)) {
    if (is.null(x = label.idents)) {
      # determine which idents to label based on flanking accessibility
      groups.use <-
        topmean[topmean$feature == features[[i]],]$group
    } else {
      # supplied list of idents to label
      groups.use <- label.idents
    }
    df.sub <- sub[(sub$feature == features[[i]]) &
                    (sub$group %in% groups.use),]
    label.df <- rbind(label.df, df.sub)
  }
  obs$label <- NA
  label.df$label <- label.df$group
  obs <- rbind(obs, label.df)

  # plot each feature separately rather than using facet
  # easier to manage the "expected" track
  df <- obs[obs$feature == features[[i]],]
  min.use <-
    ifelse(test = normalization == "subtract",
           yes = -0.5,
           no = 0.5)
  axis.min <-
    min(min.use, ymin[ymin$feature == features[[i]],]$mn)
  axis.max <- ymax[ymax$feature == features[[i]],]$mn + 0.5

  p <- ggplot(
    data = df,
    mapping = aes(
      x = position,
      y = norm.value,
      color = group,
      label = label
    )
  )
  p <- p +
    geom_line(size = linesize) +
    xlab("Distance from motif") +
    ylab(label = "Tn5 insertion\nenrichment") +
    theme_cowplot(font_size = fontsize) +
    labs(subtitle = features[[i]]) +
    ylim(c(axis.min, axis.max)) +
    guides(color = guide_legend(override.aes = list(size = 1)))
  if (label) {
    if (repel) {
      p <- p + geom_label_repel(box.padding = 0.5, show.legend = FALSE)
    } else {
      p <- p + geom_label(show.legend = FALSE)
    }
  }
  if (!is.null(alt_main_title)) {
    p <- p + labs(subtitle = alt_main_title)
  }
  if (!is.null(alt_color_title)) {
    p <- p + labs(color = alt_color_title)
  }
  p <- p +
    theme(plot.subtitle = element_markdown(hjust = 0.5)) +
    theme(legend.position = legend_pos)
  plots <- p
  if (show.expected) {
    df <- expect[expect$feature == features[[i]],]
    p1 <- ggplot(data = df,
                 mapping = aes(x = position, y = norm.value)) +
      geom_line(size = linesize) +
      xlab("Distance from motif") +
      ylab(label = "Expected\nTn5 enrichment") +
      theme_cowplot(font_size = fontsize)

    # remove x-axis labels from top plot
    p <- p + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank()
    )
    plots <-
      plot_grid(
        p,
        p1,
        nrow = 2,
        rel_heights = c(3, 1) ,
        align = "v",
        axis = "lr"
      )
  }
  return(plots)
}
