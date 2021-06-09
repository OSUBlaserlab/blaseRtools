#' Make a scatter plot of GO term associations 
#' 
#' @param simMatrix  
#' @param reducedTerms
#' @param size
#' @param addLabel
#' @param labelSize 
#' @return A ggplot 
#' @export
#' @import tidyverse ggrepel
#' @importFrom stats cmdscale as.dist
bb_goscatter <-
  function (simMatrix,
            reducedTerms,
            size = "score",
            addLabel = TRUE,
            labelSize = 4) {
    x <- stats::cmdscale(
			 as.matrix(
				   stats::as.dist(1 - simMatrix)), 
		  eig = TRUE,
                  k = 2)
    df <-
      cbind(as.data.frame(x$points), reducedTerms[match(rownames(x$points),
                                                        reducedTerms$go), c("term", "parent", "parentTerm", "size", "score")])
    p <-
      ggplot(df, aes(x = V1, y = V2, color = parentTerm)) +
      geom_point(aes(size = !!sym(size)), alpha = 0.5) +
      scale_color_discrete(guide = FALSE) +
      scale_size_continuous(guide = FALSE, range = c(0, 10)) +
      geom_text_repel(
        aes(label = parentTerm),
        segment.size = 0.25,
        data = subset(df, parent == rownames(df)),
        box.padding = grid::unit(0.5, "lines"),
        size = labelSize,
        color = "black",
        max.overlaps = 100,
        force = 2,
        seed = 1234,
        segment.curvature = -0.1,
        segment.square = TRUE,
        segment.inflect = TRUE,
        min.segment.length = 0
      ) +
      labs(x = "PCoA 1", y = "PCoA 2")
   
   return(p)
  }
