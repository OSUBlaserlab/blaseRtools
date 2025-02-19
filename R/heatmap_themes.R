#' @import cowplot
#' @import ggplot2
theme_hm_colData_top <- function(font_size = theme_get()$text$size) {
  theme_cowplot(font_size = font_size) +
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.text.y.right = element_text(),
      legend.position = "right"
    )

}

#' @import cowplot
#' @import ggplot2
theme_hm_colData_left <-
  function(font_size = theme_get()$text$size) {
    theme_cowplot(font_size = font_size) +
      theme(
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "right"
      )

  }

#' @import cowplot
#' @import ggplot2
theme_hm_rowData_left <-
  function(font_size = theme_get()$text$size) {
    theme_cowplot(font_size = font_size) +
      theme(
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        axis.text.x.top = element_text(angle = 90, vjust = 0.5),
        legend.position = "right"

      )

  }

#' @import cowplot
#' @import ggplot2
theme_hm_rowData_top <-
  function(font_size = theme_get()$text$size) {
    theme_cowplot(font_size = font_size) +
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        axis.text.y.right = element_text(),
        legend.position = "right"

      )

  }

#' @import cowplot
#' @import ggplot2
theme_hm_main <- function(font_size = theme_get()$text$size) {
  theme_cowplot(font_size = font_size) +
    theme(axis.line = element_blank(),
          plot.margin = margin(2, 0, 0, 2))

}

#' @import cowplot
#' @import ggplot2
theme_dendro <- function(font_size = theme_get()$text$size) {
  theme_nothing(font_size = font_size) +
    theme(plot.margin = margin(0, 0, 0, 0))

}
