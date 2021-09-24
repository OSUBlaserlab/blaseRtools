#' Make a dotplot of gene expression by cell population
#'
#' @param cds A cell data set object
#' @param markers A character vector of genes to plot
#' @param group_cells_by A cds colData column.  Use "multifactorial" to pick 2 categorical variables to put on X axis and to facet by.  See ordering below.
#' @param norm_method How to normalize gene expression.  Default = "log"
#' @param lower_threshold Lower cutoff for gene expression
#' @param max.size The maximum size of the dotplot
#' @param group_ordering Defaults to "biclustering" method from pheatmap.  Optionally will take a vector of group values to set the axis order explicitly.  If using group_cells_by = "multifactorial" you will need a df to define facet and axis levels.  See example.
#' @param gene_ordering Optional vector of gene names to order the plot.
#' @param pseudocount Add to zero expressors.  Default = 1
#' @param scale_max Expression scale max
#' @param scale_min Expression scale min
#' @param colorscale_name Label for the color scale
#' @param sizescale_name Label for the size scale
#' @param ... Additional parameters to pass to facet_wrap.
#' @return A ggplot
#' @export
#' @import tidyverse monocle3
#' @importFrom stats as.dist cor
bb_gene_dotplot <- function(cds,
                            markers,
                            group_cells_by,
                            reduction_method = "UMAP",
                            norm_method = c("log", "size_only"),
                            lower_threshold = 0,
                            max.size = 10,
                            group_ordering = "bicluster",
                            gene_ordering = NULL,
                            pseudocount = 1,
                            scale_max = 3,
                            scale_min = -3,
                            colorscale_name = NULL,
                            sizescale_name = NULL,
                            ...)
{
  norm_method = match.arg(norm_method)
  gene_ids = as.data.frame(fData(cds)) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(rowname %in% markers | gene_short_name %in%
                    markers) %>%
    dplyr::pull(rowname)
  major_axis <- 2
  minor_axis <- 1
  exprs_mat <- t(as.matrix(exprs(cds)[gene_ids,]))
  exprs_mat <- reshape2::melt(exprs_mat)
  colnames(exprs_mat) <- c("Cell", "Gene", "Expression")
  exprs_mat$Gene <- as.character(exprs_mat$Gene)

  if (group_cells_by == "multifactorial") {
    multivar <- paste0(unique(group_ordering$variable)[1],"_AND_",unique(group_ordering$variable)[2])
    multivar_val <- paste0(unique(group_ordering$value))
    colData(cds)[,multivar] <- paste0(colData(cds)[,unique(group_ordering$variable)[1]],"_AND_",colData(cds)[,unique(group_ordering$variable)[2]])
    cell_group <- colData(cds)[,multivar]
  } else {
    cell_group <- colData(cds)[, group_cells_by]
  }
  names(cell_group) = colnames(cds)
  exprs_mat$Group <- cell_group[exprs_mat$Cell]
  exprs_mat = exprs_mat %>% dplyr::filter(is.na(Group) == FALSE)

  ExpVal <- exprs_mat %>% dplyr::group_by(Group, Gene) %>%
    dplyr::summarize(
      mean = mean(log(Expression + pseudocount)),
      percentage = sum(Expression > lower_threshold) / length(Expression)
    )
  ExpVal$mean <- ifelse(ExpVal$mean < scale_min, scale_min,
                        ExpVal$mean)
  ExpVal$mean <- ifelse(ExpVal$mean > scale_max, scale_max,
                        ExpVal$mean)
  ExpVal$Gene <- fData(cds)[ExpVal$Gene, "gene_short_name"]
    res <-
      reshape2::dcast(ExpVal[, 1:4], Group ~ Gene, value.var = colnames(ExpVal)[2 +
                                                                                  major_axis])
    group_id <- res[, 1]
    res <- res[,-1]
    row.names(res) <- group_id
    row_dist <- stats::as.dist((1 - stats::cor(t(res))) / 2)
    row_dist[is.na(row_dist)] <- 1
    col_dist <- stats::as.dist((1 - stats::cor(res)) / 2)
    col_dist[is.na(col_dist)] <- 1
    ph <-
      pheatmap::pheatmap(
        res,
        useRaster = T,
        cluster_cols = TRUE,
        cluster_rows = TRUE,
        show_rownames = F,
        show_colnames = F,
        clustering_distance_cols = col_dist,
        clustering_distance_rows = row_dist,
        clustering_method = "ward.D2",
        silent = TRUE,
        filename = NA
      )
    ExpVal$Gene <-
      factor(ExpVal$Gene, levels = colnames(res)[ph$tree_col$order])
    ExpVal$Group <-
      factor(ExpVal$Group, levels = row.names(res)[ph$tree_row$order])
  if (group_ordering != "bicluster" && group_cells_by != "multifactorial") {
    ExpVal$Group <-
      factor(ExpVal$Group, levels = group_ordering)
  }
  if (!is.null(gene_ordering)) {
    ExpVal$Gene <- factor(ExpVal$Gene, levels = gene_ordering)
  }

  if (group_cells_by != "multifactorial") {
    g <-
      ggplot(ExpVal, aes(y = Gene, x = Group)) +
      geom_point(aes(colour = mean,
                     size = percentage)) +
      viridis::scale_color_viridis(name = ifelse(is.null(colorscale_name),
                                                 "log(mean + 0.1)",
                                                 colorscale_name)) +
      scale_size(
        name = ifelse(is.null(sizescale_name), "proportion", sizescale_name),
        range = c(0, max.size))
    return(g)
  }
  if (group_cells_by == "multifactorial") {
    facet_choice <- group_ordering %>% filter(aesthetic == "facet") %>% pull(value) %>% paste(collapse = "|")
    axis_choice <- group_ordering %>% filter(aesthetic == "axis") %>% pull(value) %>% paste(collapse = "|")
    ExpVal <-
      ExpVal %>%
      mutate(facet = str_extract(Group, pattern = facet_choice)) %>%
      mutate(facet = factor(facet, levels = group_ordering %>% filter(aesthetic == "facet") %>% arrange(level) %>% pull(value))) %>%
      mutate(axis = str_extract(Group, pattern = axis_choice)) %>%
      mutate(axis = factor(axis, levels = group_ordering %>% filter(aesthetic == "axis") %>% arrange(level) %>% pull(value)))
      g <-
      ggplot(ExpVal, mapping = aes(y = Gene, x = axis)) +
      geom_point(aes(colour = mean,
                     size = percentage)) +
      viridis::scale_color_viridis(name = ifelse(is.null(colorscale_name),
                                                 "log(mean + 0.1)",
                                                 colorscale_name
                                                                )) +
      scale_size(
        name = ifelse(is.null(sizescale_name), "proportion", sizescale_name),
        range = c(0, max.size)
      ) +
      facet_wrap(facets = vars(facet),...) +
      theme(strip.background = element_blank())
     return(g)
  }

}
