#' A function to reduce go terms by semantic similarity
#' 
#' @param  x A list go term enrichment results produced by bb_goenrichment.
#' @param  reduce_threshold The degree of term reduction. 0 to 1.  Higher is more reduction. 
#' @param  go_db The database to query.  Choose from c("org.Hs.eg.db", "org.Dr.eg.db", "org.Mm.eg.db", ...). 
#' @return A list of items for downstream plotting
#' @export
#' @import rrvgo
bb_gosummary <- function(x, 
			 reduce_threshold = 0.8,
			 go_db = c("org.Hs.eg.db", "org.Dr.eg.db", "org.Mm.eg.db")) {
    simMatrix <-
      calculateSimMatrix(x = x[[3]]$GO.ID,
                         ont = "BP",
                         orgdb = go_db)
    scores <- setNames(-log10(ifelse(is.na(
      as.numeric(x[[3]]$classicFisher)),
      1e-30,
      as.numeric(x[[3]]$classicFisher)
    )),
    x[[3]]$GO.ID)
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold = reduce_threshold,
                                    orgdb = go_db)
    returnlist <- list(simMatrix, scores, reducedTerms)
    names(returnlist) <- c("simMatrix", "scores", "reducedTerms")
    return(returnlist)
}
