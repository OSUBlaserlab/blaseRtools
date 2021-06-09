#' @title GO Term Enrichment
#' @description A function to find enriched go terms from a query list of gene names relative to a reference list of gene names.
#' @param query A vector of gene names
#' @param reference The background gene list.  Usually will be rowData(cds_main).
#' @param group_pval P value to determine enrichment.  Default: 0.01.
#' @param go_db PARAM_DESCRIPTION, Default: c("org.Hs.eg.db", "org.Dr.eg.db", "org.Mm.eg.db")
#' @return A list of items including the enrichment results.
#' @import tidyverse topGO
#' @rdname bb_goenrichment
bb_goenrichment <- function(query,
                            reference,
                            group_pval = 0.01,
                            go_db = c("org.Hs.eg.db", "org.Dr.eg.db", "org.Mm.eg.db")) {
  genes <- query
  genes_named <- reference %>%
    as_tibble() %>%
    mutate(selected = ifelse(gene_short_name %in% genes, 1, 0)) %>%
    pull(selected)
  names(genes_named) <- reference %>%
    as_tibble() %>%
    pull(gene_short_name)

  sampleGOdata <- new(
    "topGOdata",
    description = "Simple session",
    ontology = "BP",
    allGenes = genes_named,
    geneSel = selector,
    nodeSize = 10,
    annot = annFUN.org,
    mapping = go_db,
    ID = "symbol"
  )

  resultFisher <-
    runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

  res_table <- GenTable(
    sampleGOdata,
    classicFisher = resultFisher,
    orderBy = "classicFisher",
    ranksOf = "classicFisher",
    topNodes = 100
  ) %>%
    as_tibble(rownames = "Rank")

  resultFisher_tbl <-
    tibble(goterm = names(resultFisher@score),
           pval = resultFisher@score)
  return_list <- list(sampleGOdata, resultFisher, res_table)
  names(return_list) <- c("sampleGOdata", "resultFisher", "res_table")

  return(return_list)
}


selector <- function(theScore) {
  return (theScore == 1)
}
