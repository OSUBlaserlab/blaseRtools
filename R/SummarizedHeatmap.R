setOldClass("dendrogram", prototype = "dendrogram")

# class definition --------------------
#' @rdname SummarizedHeatmap
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom stats dendrogram
.SummarizedHeatmap <- setClass(
  "SummarizedHeatmap",
  slots = representation(
    colDendro = "dendrogram",
    rowDendro = "dendrogram"
  ),
  contains = "SummarizedExperiment"
)

# constructor ----------------------
#' @title An S4 Class for Holding Heatmap Data
#' @description
#' This is an S4 Class that is part of a solution to optimize plotting heatmaps.  This class is derived from the SummarizedExperiment class.  It inherits structure and methods and adds some structures.
#'
#' Like the SummarizedExperiment Class it is built around a matrix.  Like the SingleCellExperiment and cell_data_set classes which are also derived from SummarizedExperiment, SummarizedHeatmap holds metadata about the columns and rows of the matrix.  This enables plotting useful annotation information with a set of plotting functions (bb_plot_heatmap...)
#'
#' Use the SummarizedHeatmap constructor to make an instance of the class from a matrix.  Use colData and rowData to get or set these values.  Internal validity checks will ensure the columns and rows match.
#'
#' New to this object are colDendro and rowDendro slots.  These hold hierarchical clustering information used for ordering the heatmap plot and plotting the dendrogrms.  These are generated automatically when the object is created.
#' @param mat A matrix to build the object from.
#' @param ... other arguments to pass into SummarizedExperiment
#' @return A SummarizedHeatmap object
#' @examples
#'
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#' mat <- matrix(rnorm(100), ncol=5)
#' colnames(mat) <- letters[1:5]
#' rownames(mat) <- letters[6:25]
#' test_sh <- SummarizedHeatmap(mat)
#' colData(test_sh)$sample_type <- c("vowel", "consonant", "consonant", "consonant", "vowel")
#' colData(test_sh)$sample_type2 <- c("vowel2", "consonant2", "consonant2", "consonant2", "vowel2")
#' isVowel <- function(char) char %in% c('a', 'e', 'i', 'o', 'u')
#' rowData(test_sh)$feature_type <- ifelse(isVowel(letters[6:25]), "vowel", "consonant")
#' rowData(test_sh)$feature_type2 <- paste0(rowData(test_sh)$feature_type, "2")

#'  }
#' }
#' @seealso
#'  \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}, \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'  \code{\link[S4Vectors]{DataFrame-class}}, \code{\link[S4Vectors]{S4VectorsOverview}}
#' @rdname SummarizedHeatmap
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
SummarizedHeatmap <- function(mat, ...) {
  se <-
    SummarizedExperiment::SummarizedExperiment(
      assays = list(matrix = mat),
      rowData = S4Vectors::DataFrame(row.names = rownames(mat)),
      colData = S4Vectors::DataFrame(row.names = colnames(mat)),
      ...
    )
  cd <- as.dendrogram(hclust(dist(t(mat)), "ave"))
  rd <- as.dendrogram(hclust(dist(mat), "ave"))
  obj <-
    .SummarizedHeatmap(se,
                       colDendro = cd,
                       rowDendro = rd)
}


# validity ------------------------------------
S4Vectors::setValidity2("SummarizedHeatmap", function(object) {
  msg <- NULL

  if (SummarizedExperiment::assayNames(object)[1] != "matrix") {
    msg <- c(msg, "'matrix' must be first assay")
  }

  if (is.null(msg)) {
    TRUE
  } else
    msg
})

# getters ----------------------------

#' @export
setGeneric("colDendro", function(x, ...)
  standardGeneric("colDendro"))

#' @export
#' @importFrom SummarizedExperiment assay
setMethod("colDendro", "SummarizedHeatmap", function(x) {
  x@colDendro

})


#' @export
setGeneric("rowDendro", function(x, ...)
  standardGeneric("rowDendro"))

#' @export
#' @importFrom SummarizedExperiment assay
setMethod("rowDendro", "SummarizedHeatmap", function(x) {
  x@rowDendro

})

#' @export
#' @importMethodsFrom SummarizedExperiment rowData
setMethod("colData", "SummarizedHeatmap", function(x, ...) {
  out <- callNextMethod()
  out
  # as_tibble(out)
})

#' @export
#' @importMethodsFrom SummarizedExperiment rowData
setMethod("rowData", "SummarizedHeatmap", function(x, ...) {
  out <- callNextMethod()
  out
  # as_tibble(out)
})

# setters ----------------------------

#' @export
setGeneric("rowData<-", function(x, ..., value)
  standardGeneric("rowData<-"))

setReplaceMethod("rowData", "SummarizedHeatmap", function(x, value) {
  x@elementMetadata <- value
  validObject(x)
  x
})


#' @export
setGeneric("colData<-", function(x, ..., value)
  standardGeneric("colData<-"))

setReplaceMethod("colData", "SummarizedHeatmap", function(x, value) {
  x@colData <- value
  validObject(x)
  x
})


