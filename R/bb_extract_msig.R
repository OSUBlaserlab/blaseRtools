#' Extract MSIGDB Gene Sets
#'
#' @description Use this to extract gene sets from the MSIGDB.  Most gene sets are known by "STANDARD_NAME".  You can filter the gene set list by supplying a named filter list to the bb_extract_msig function.  The name of each list element should be one of the metadata column names and the list element contents should be the values to filter for.  Filtering works in an additive way, meaning if you supply a filter list with two elements it will extract gene sets passing filters 1 AND 2.
#' @param filter_list A named list to filter the MSIGDB by.  Defaults to NULL which will return the whole MSIGDB
#' @param return_form Select from a list of gene ids or a list of gene names by gene set.  This is a useful format for the fgsea package.  Alternatively "tibble" can be select and all filtered gene sets will be bound into a long-form (tidy) tibble.
#' @return Gene set as a list or tibble.
#' @export
#' @import blaseRdata tidyverse
bb_extract_msig <- function(filter_list = NULL,
                            return_form = c("id_list",
                                            "name_list",
                                            "tibble")) {
    meta <- msigdb_geneset_metadata
    if (!is.null(filter_list)) {
      stopifnot(
        "Your filter must be either NULL or a named list.\n  Names must correspond to msigdb_geneset_metadata column names." = all(
          names(filter_list) %in% colnames(msigdb_geneset_metadata)
        )
      )
      for (i in 1:length(filter_list)) {
        meta <-
          meta %>%
          filter(!!sym(names(filter_list[i])) %in% filter_list[[i]])
      }
      names_to_get <- meta$STANDARD_NAME
    } else {
      names_to_get <- meta$STANDARD_NAME
    }

    res <- msigdb_genesets[names_to_get]

    if (return_form == "tibble") {
      res <- bind_rows(res, .id = "gene_set")
    } else if (return_form == "id_list") {
      temp <- map(res, ~ .x %>%
                    pull(id)) %>%
        set_names(nm = names(res))
      res <- temp
    } else {
      temp <- map(res, ~ .x %>%
                    filter(!is.na(gene_short_name)) %>%
                    pull(gene_short_name)) %>%
        set_names(nm = names(res))
      res <- temp
    }
    return(res)

  }
