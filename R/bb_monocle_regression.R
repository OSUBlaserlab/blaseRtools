#' A function to perform regression on single cell data.
#'
#' @param cds A cell data set object.
#' @param gene_or_genes Genes to regress by.
#' @param stratification_variable Optional colData column to subset the cds by internal to the function.
#' @param stratification_value Optional value to stratify by.
#' @param form The regression formula in the form of "~var1+var2+..."
#' @param linking_function For the generalized linear model.
#' @return A tibble containing the regression results.
#' @export
#' @import tidyverse monocle3
bb_monocle_regression <-
  function(cds,
           gene_or_genes,
           stratification_variable = NULL,
           stratification_value = NULL,
           form,
           linking_function = "negbinomial") {
    if (!is.null(stratification_variable)) {
      cds <-
        cds[rowData(cds)$gene_short_name %in% gene_or_genes,
            colData(cds)[[stratification_variable]] == stratification_value]
    } else {
      cds <- cds[rowData(cds)$gene_short_name %in% gene_or_genes,]
    }
    gene_fits <-
      fit_models(
        cds = cds,
        model_formula_str = form,
        expression_family = linking_function,
        cores = 1
      )
    fit_coefs <- coefficient_table(gene_fits)

    if (!is.null(stratification_variable)) {
      fit_coefs_return <-
        fit_coefs %>%
        filter(term != "(Intercept)") %>%
        mutate(stratification = stratification_value, formula = form) %>%
        select(id, gene_short_name, stratification, formula, term, estimate, p_value, q_value)
        # select(stratification, formula, term, estimate, p_value, q_value)
    } else {
      fit_coefs_return <-
        fit_coefs %>%
        filter(term != "(Intercept)") %>%
        mutate(stratification = "no stratification", formula = form) %>%
        select(id, gene_short_name, stratification, formula, term, estimate, p_value, q_value)

    }
    return(fit_coefs_return)
  }
