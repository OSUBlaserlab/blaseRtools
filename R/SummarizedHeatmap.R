setOldClass("dendrogram", prototype = "dendrogram")
setClassUnion("dendrogramOrNULL", c("dendrogram", "NULL"))

# class definition --------------------
#' @rdname SummarizedHeatmap
#' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.SummarizedHeatmap <- setClass(
  "SummarizedHeatmap",
  slots = list(
    colDendro = "dendrogramOrNULL",
    rowDendro = "dendrogramOrNULL",
    colOrder = "character",
    rowOrder = "character"
  ),
  prototype = list(
    colDendro = NULL,
    rowDendro = NULL,
    colOrder = NULL,
    rowOrder = NULL
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
#'
#' In order to manually set the order of the columns or rows, supply values to the rowOrder or colOrder parameters.  This will prevent creation of dendrograms for the respective colums or rows.
#'
#' @param mat A matrix to build the object from.
#' @param colOrder A character string corresponding to matrix column names.
#' @param rowOrder A character string corresponding to matrix row names.
#' @param cluster_method Clusterihng algorithm.  See stats::hclust.
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
SummarizedHeatmap <- function(
    mat,
    colOrder = NULL,
    rowOrder = NULL,
    cluster_method = "ave",
    ...) {
  se <-
    SummarizedExperiment::SummarizedExperiment(
      assays = list(matrix = mat),
      rowData = S4Vectors::DataFrame(row.names = rownames(mat)),
      colData = S4Vectors::DataFrame(row.names = colnames(mat)),
      ...
    )

  cd <- as.dendrogram(hclust(dist(t(mat)), method = cluster_method))
  rd <- as.dendrogram(hclust(dist(mat), method = cluster_method))

  if (!is.null(colOrder)) {
    co <- colOrder
    cd <- NULL
  } else {
    ddata <- ggdendro::dendro_data(cd, type = "rectangle")
    co <- ddata$labels$label
  }

  if (!is.null(rowOrder)) {
    ro <- rowOrder
    rd <- NULL
  } else {
    ddata <- ggdendro::dendro_data(rd, type = "rectangle")
    ro <- ddata$labels$label
  }

  obj <-
    .SummarizedHeatmap(se,
      colDendro = cd,
      rowDendro = rd,
      colOrder = co,
      rowOrder = ro
    )
}


# validity ------------------------------------
S4Vectors::setValidity2("SummarizedHeatmap", function(object) {
  msg <- NULL

  if (SummarizedExperiment::assayNames(object)[1] != "matrix") {
    msg <- c(msg, "'matrix' must be first assay")
  }


  if (!is.null(rowDendro(object))) {
    dg <- rowDendro(object)
    ddata <- ggdendro::dendro_data(dg, type = "rectangle")
    if (any(rowOrder(object) != ddata$labels$label)) {
      msg <- c(msg, "The rowOrder slot does not match the row dendrogram")
    }

  }

  if (!is.null(colDendro(object))) {
    dh <- colDendro(object)
    hdata <- ggdendro::dendro_data(dh, type = "rectangle")
    if (any(colOrder(object) != hdata$labels$label)) {
      msg <- c(msg, "The colOrder slot does not match the row dendrogram")
    }

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

#' @export
setGeneric("colOrder", function(x, ...)
  standardGeneric("colOrder"))

#' @export
#' @importFrom SummarizedExperiment assay
setMethod("colOrder", "SummarizedHeatmap", function(x) {
  x@colOrder

})


#' @export
setGeneric("rowOrder", function(x, ...)
  standardGeneric("rowOrder"))

#' @export
#' @importFrom SummarizedExperiment assay
setMethod("rowOrder", "SummarizedHeatmap", function(x) {
  x@rowOrder

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




