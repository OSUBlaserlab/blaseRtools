#' @title Infer Cell-Cell Communication
#' @description Use this function to identify ligand/receptor pairs expressed by cell clusters in human, mouse or zebrafish single cell data.  A CellChat object is generated which can be used to visualize these connections using bb_cellchat_heatmap or other tools from package CellChat.
#' @param cds The cell data set object.  It should usually be pre-filtered to conatin a single biological sample.
#' @param group_var The cell metadata column identifying cell groups for cell-cell communication inference.
#' @param n_cores Number of cores for the analysis, Default: 12
#' @param species Species for the assay, Default: c("human", "mouse", "zebrafish")
#' @param min_cells Cell clusters smaller than this value will be ignored., Default: 10
#' @param prob_type Methods for computing the average gene expression per cell group. By default = "triMean", producing fewer but stronger interactions; When setting ‘type = "truncatedMean"', a value should be assigned to ’trim', producing more interactions, Default: c("triMean", "truncatedMean", "median")
#' @param prob_trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed if using truncatedMean, Default: NULL
#' @param project Whether or not to smooth gene expression, Default: TRUE
#' @param pop_size_arg Whether consider the proportion of cells in each group across all sequenced cells. Set population.size = FALSE if analyzing sorting-enriched single cells, to remove the potential artifact of population size. Set population.size = TRUE if analyzing unsorted single-cell transcriptomes, with the reason that abundant cell populations tend to send collectively stronger signals than the rare cell populations., Default: TRUE
#' @return A CellChat object
#' @details see github::sqjin/CellChat
#' @seealso
#'  \code{\link[monocle3]{normalized_counts}}
#'  \code{\link[dplyr]{mutate-joins}},\code{\link[dplyr]{pull}}
#'  \code{\link[tibble]{tibble}}
#'  \code{\link[CellChat]{createCellChat}},\code{\link[CellChat]{CellChatDB.human}},\code{\link[CellChat]{CellChatDB.mouse}},\code{\link[CellChat]{character(0)}},\code{\link[CellChat]{subsetData}},\code{\link[CellChat]{identifyOverExpressedGenes}},\code{\link[CellChat]{identifyOverExpressedInteractions}},\code{\link[CellChat]{projectData}},\code{\link[CellChat]{computeCommunProb}},\code{\link[CellChat]{filterCommunication}},\code{\link[CellChat]{computeCommunProbPathway}},\code{\link[CellChat]{aggregateNet}}
#'  \code{\link[SingleCellExperiment]{character(0)}}
#'  \code{\link[future]{plan}}
#' @rdname bb_cellchat
#' @export
#' @importFrom monocle3 normalized_counts
#' @importFrom dplyr left_join pull
#' @importFrom tibble tibble
#' @importFrom CellChat createCellChat CellChatDB.human CellChatDB.mouse CellChatDB.zebrafish subsetData identifyOverExpressedGenes identifyOverExpressedInteractions projectData computeCommunProb filterCommunication computeCommunProbPathway aggregateNet
#' @importFrom SingleCellExperiment colData
#' @importFrom future plan
bb_cellchat <-
  function(cds,
           group_var,
           n_cores = 12,
           species = c("human", "mouse", "zebrafish"),
           min_cells = 10,
           prob_type = c("triMean", "truncatedMean", "median"),
           prob_trim = NULL,
           project = TRUE,
           pop_size_arg = TRUE) {
    species <- match.arg(species)
    prob_type <- match.arg(prob_type)

    # get the expression data
    cellchat_data <- monocle3::normalized_counts(cds)
    row.names(cellchat_data) <-
      dplyr::left_join(tibble::tibble(feature_id = row.names(cellchat_data)),
                       bb_rowmeta(cds_main)) |>
      dplyr::pull(gene_short_name)

    # set up the cellchat object
    cellchat <- CellChat::createCellChat(
      object = cellchat_data,
      meta = droplevels(as.data.frame(SingleCellExperiment::colData(cds))),
      group.by = group_var
    )

    # Set the ligand-receptor interaction database
    if (species == "human") {
      CellChatDB <-
        CellChat::CellChatDB.human
    } else if (species == "mouse") {
      CellChatDB <-
        CellChat::CellChatDB.mouse
    } else if (species == "zebrafish") {
      CellChatDB <-
        CellChat::CellChatDB.zebrafish
    } else {
      stop("The species must be one of human, mouse or zebrafish.")
    }

    # use all CellChatDB for cell-cell communication analysis
    CellChatDB.use <- CellChatDB

    # set the used database in the object
    cellchat@DB <- CellChatDB.use

    # subset the expression data of signaling genes for saving computation cost
    cellchat <-
      CellChat::subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multiprocess", workers = n_cores) # do parallel

    cellchat <- CellChat::identifyOverExpressedGenes(cellchat)

    cellchat <-
      CellChat::identifyOverExpressedInteractions(cellchat)

    # project gene expression data onto PPI network (optional)
    if (project) {
      if (species == "human") {
        cellchat <- CellChat::projectData(cellchat, PPI.human)
      } else if (species == "mouse") {
        cellchat <- CellChat::projectData(cellchat, PPI.human)
      } else {
        stop("Can only project if the species is human or mouse.")
      }

    }


    # Compute the communication probability and infer cellular communication network
    cellchat <-
      CellChat::computeCommunProb(cellchat, population.size = pop_size_arg, type = prob_type)

    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <-
      CellChat::filterCommunication(cellchat, min.cells = min_cells)


    # Infer the cell-cell communication at a signaling pathway level
    cellchat <- CellChat::computeCommunProbPathway(cellchat)

    # Calculate the aggregated cell-cell communication network
    cellchat <- CellChat::aggregateNet(cellchat)

    future::plan("sequential")
    return(cellchat)
  }


#' @title Make a Complex Heatmap From a CellChat Object With Sensible Defaults
#' @description This will generate heatmap from a CellChat object using ComplexHeatmap::Heatmap.  Options are provided to filter for sender and receiver cells, to generate simple marginal annotations and for aesthetic control.
#' @param object The CellChat object to plot
#' @param source_filter Optional filter for source cell clusters from the object metadata., Default: NULL
#' @param target_filter Optional filter for target cell clusters from the object metadata., Default: NULL
#' @param colors Color scale endpoints, Default: c("transparent", "red3")
#' @param rowanno Options for simple row annotation; must be one of c(NULL, "Annotation", "Pathway")
#' @param rowanno_colors Optional colors to replace the poor color selections from Complex heatmap.  Must be supplied as a named list with one element each for "Annotation" and "Pathway".  Not required if not showing these annotations.  The list should be of the form:  list(Annotation = c("name1" = "color value1", "name2" = "color_value2")), Default: NULL
#' @param colanno Options for simple column annotation; must be one of c(NULL, "Source", "Target")
#' @param colanno_colors See rowanno_colors, Default: NULL
#' @param pval_filter Filter for significance of associations.  CellChat returns pvalues of 0, 0.01, and 0.05; this function will filter and retain values less than or equal to the provided value. Default: 0.05
#' @param heatmap_name Name for the main color scale of the heatmap, Default: 'InteractionScore'
#' @param heatmap_show_row_dend Show row dendrograms? Default: TRUE
#' @param heatmap_row_dend_width Width of row dendrograms Default: unit(5, "mm")
#' @param heatmap_show_column_dend Show column dendrograms?' Default: TRUE
#' @param heatmap_column_dend_height Height of column dendrograms. Default: unit(5, "mm")
#' @param heatmap_row_names_gp Row name graphical params, Default: gpar(fontsize = 10)
#' @param heatmap_column_names_gp Column name graphical params, Default: gpar(fontsize = 10)
#' @param heatmap_column_names_rot Column name rotation, Default: 90
#' @param heatmap_column_title Column title, Default: NULL
#' @param heatmap_column_title_gp Column title graphical params, Default: gpar(fontsize = 12, fontface = "bold")
#' @param col_anno_name_gp Column annotation name graphical params, Default: gpar(fonmtsize = 10, fontface = "bold")
#' @param row_anno_name_gp Row annotation name graphical params, Default: gpar(fontsize = 10, fontface = "bold")
#' @return a heatmap as a grid object; plot using cowplot::plot_grid
#' @details see github::sqjin/CellChat
#' @seealso
#'  \code{\link[CellChat]{subsetCommunication}}
#'  \code{\link[tibble]{as_tibble}},\code{\link[tibble]{c("tibble", "tibble")}},\code{\link[tibble]{rownames}}
#'  \code{\link[dplyr]{filter}},\code{\link[dplyr]{mutate}},\code{\link[dplyr]{select}},\code{\link[dplyr]{mutate-joins}},\code{\link[dplyr]{group_by}},\code{\link[dplyr]{summarise}}
#'  \code{\link[tidyr]{pivot_wider}}
#'  \code{\link[ComplexHeatmap]{rowAnnotation}},\code{\link[ComplexHeatmap]{columnAnnotation}},\code{\link[ComplexHeatmap]{draw-dispatch}},\code{\link[ComplexHeatmap]{Heatmap}}
#'  \code{\link[circlize]{colorRamp2}}
#'  \code{\link[grid]{grid.grab}}
#' @rdname bb_cellchat_heatmap
#' @export
#' @importFrom CellChat subsetCommunication
#' @importFrom tibble as_tibble tibble column_to_rownames
#' @importFrom dplyr filter mutate select left_join group_by summarise
#' @importFrom tidyr pivot_wider
#' @importFrom ComplexHeatmap rowAnnotation columnAnnotation draw Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom grid grid.grabExpr
bb_cellchat_heatmap <- function(object,
                                source_filter = NULL,
                                target_filter = NULL,
                                colors = c("transparent", "red3"),
                                rowanno = c(NULL, "Annotation", "Pathway"),
                                rowanno_colors = NULL,
                                colanno = c(NULL, "Source", "Target"),
                                colanno_colors = NULL,
                                pval_filter = 0.05,
                                heatmap_name = "Interaction\nScore",
                                heatmap_show_row_dend = TRUE,
                                heatmap_row_dend_width = unit(5, "mm"),
                                heatmap_show_column_dend = TRUE,
                                heatmap_column_dend_height = unit(5, "mm"),
                                heatmap_row_names_gp = gpar(fontsize = 10),
                                heatmap_column_names_gp = gpar(fontsize = 10),
                                heatmap_column_names_rot = 90,
                                heatmap_column_title = NULL,
                                heatmap_column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                col_anno_name_gp = gpar(fontsize = 10, fontface = "bold"),
                                row_anno_name_gp = gpar(fontsize = 10, fontface = "bold")) {
  stopifnot("The object must be of class CellChat." = "CellChat" %in% class(object))
  stopifnot("You must select valid colors." = all(areColors(colors)))
  stopifnot("You can only select 2 colors." = length(colors) == 2)
  stopifnot(
    "You can only annotate rows by c(NULL, 'Annotation', 'Pathway')" = rowanno %in% c(NULL, "Annotation", "Pathway")
  )
  stopifnot(
    "You can only annotate columns by c(NULL, 'Source', 'Target')" = colanno %in% c(NULL, "Source", "Target")
  )

  mtx <- CellChat::subsetCommunication(object = object) |>
    tibble::as_tibble() |>
    dplyr::filter(pval <= pval_filter)
  if (!is.null(source_filter))
    mtx <- mtx |> dplyr::filter(source %in% source_filter)
  if (!is.null(target_filter))
    mtx <- mtx |> dplyr::filter(target %in% target_filter)
  mtx <- mtx |>
    dplyr::mutate(source_target = paste0(source, " -> ", target)) |>
    dplyr::select(source_target, interaction_name_2, prob) |>
    tidyr::pivot_wider(
      names_from = interaction_name_2,
      values_from = prob,
      values_fill = 0
    ) |>
    bb_tbl_to_matrix() |>
    t()

  if (!is.null(rowanno)) {
    row_anno_df <- tibble::tibble(rownames = rownames(mtx)) |>
      dplyr::left_join(
        CellChat::subsetCommunication(object) |>
          dplyr::group_by(rownames = interaction_name_2, Annotation = annotation, Pathway = pathway_name) |>
          dplyr::summarise()
      )
    if (rowanno == "Annotation") row_anno_df <- dplyr::select(row_anno_df, c(rownames, Annotation))
    if (rowanno == "Pathway") row_anno_df <- dplyr::select(row_anno_df, c(rownames, Pathway))

    row_anno_df <- row_anno_df |>
      tibble::column_to_rownames("rownames") |>
      as.data.frame()


    rowanno <-
      ComplexHeatmap::rowAnnotation(df = row_anno_df,
                                    col = rowanno_colors,
                                    annotation_name_gp = row_anno_name_gp)

  }

  if (!is.null(colanno)) {
    col_anno_df <- tibble::tibble(colnames = colnames(mtx)) |>
      dplyr::left_join(
        CellChat::subsetCommunication(object) |>
          dplyr::mutate(source_target = paste0(source, " -> ", target)) |>
          dplyr::group_by(colnames = source_target, Source = source, Target = target) |>
          dplyr::summarise()
      )

    if (colanno == "source") col_anno_df <- dplyr::select(col_anno_df, c(colnames, Source))
    if (colanno == "target") col_anno_df <- dplyr::select(col_anno_df, c(colnames, Target))

    col_anno_df <- col_anno_df |>
      tibble::column_to_rownames("colnames") |>
      as.data.frame()

    colanno <-
      ComplexHeatmap::columnAnnotation(df = col_anno_df,
                                       col = colanno_colors,
                                       annotation_name_gp = col_anno_name_gp)



  }

  col_fun <-
   circlize::colorRamp2(breaks = c(min(mtx), max(mtx)),
                               colors = colors)
  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      ComplexHeatmap::Heatmap(
        mtx,
        name = heatmap_name,
        col = col_fun,
        left_annotation = rowanno,
        top_annotation = colanno,
        show_row_dend = heatmap_show_row_dend,
        row_dend_width = heatmap_row_dend_width,
        show_column_dend = heatmap_show_column_dend,
        column_dend_height = heatmap_column_dend_height,
        row_names_gp = heatmap_row_names_gp,
        column_names_gp = heatmap_column_names_gp,
        column_names_rot = heatmap_column_names_rot,
        column_title = heatmap_column_title,
        column_title_gp = heatmap_column_title_gp

      ),
      heatmap_legend_side = "right",
      annotation_legend_side = "bottom",
      legend_grouping = "original"
    ),
    wrap = TRUE
  )



}



areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  })
}
