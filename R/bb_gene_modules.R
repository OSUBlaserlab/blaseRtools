#' A function to generate gene modules and add them to the CDS rowData
#'
#' Based on Monocle3's gene module functions.  Implemented with default values.
#'
#' @param cds A cell data set object
#' @param n_cores Number of processor cores to use for the analysis
#' @return A modified cell data set object with a new rowData column holding gene module data.
#' @export
#' @import tidyverse monocle3
bb_gene_modules <- function(cds, 
			    n_cores = 8) {
  graph_test_res <-
    graph_test(
      cds = cds,
      neighbor_graph = "knn",
      cores = n_cores,
      verbose = TRUE
    )
  
  mod_deg_ids <- row.names(subset(graph_test_res, q_value < 0.05))
  
  gene_module_df <-
    find_gene_modules(cds[mod_deg_ids, ], cores = 39)
  
  rowData(cds)$module <-
    left_join(as_tibble(rowData(cds)),
              gene_module_df) %>%
    pull(module)
  
  rowData(cds)$module_labeled <-
    left_join(as_tibble(rowData(cds)) %>% select(-module),
              gene_module_df) %>%
    mutate(module_labeled = ifelse(is.na(module),"No Module", paste0("Module ",module))) %>%
    mutate(module_labeled = factor(module_labeled, levels = str_sort(unique(module_labeled),numeric = TRUE))) %>%
    pull(module_labeled)
  
  rowData(cds)$supermodule <-
    left_join(as_tibble(rowData(cds)) %>% select(-c(module, module_labeled)),
              gene_module_df) %>%
    pull(supermodule)
  
  rowData(cds)$supermodule_labeled <-
    left_join(as_tibble(rowData(cds)) %>% select(-c(module, module_labeled, supermodule)),
              gene_module_df) %>%
    mutate(supermodule_labeled = ifelse(is.na(supermodule), "No Supermodule", paste0("Supermodule ",supermodule))) %>%
    mutate(supermodule_labeled = factor(supermodule_labeled, levels = str_sort(unique(supermodule_labeled),numeric = TRUE))) %>%
    pull(supermodule_labeled) %>%
    fct_shift()
  
  return(cds)
}
