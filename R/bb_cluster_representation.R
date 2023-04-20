#' @title Cluster Representation By Fisher Exact Test Per Cell
#' @description Use this function to determine the differential representation of cells in clusters.  It will determine fold change in a single experimental class over a single control or reference class.  This value is normalized to the number of cells captured in all clusters from the class.  Significance is determined using Fisher's exact test.  This test may overestimate significance in large data sets.  In this case, bb_cluster_representation2 may be more robust.
#' @param cds A cell data set object
#' @param cluster_var The CDS cell metadata column holding cluster data.  There can be any number of clusters in this column.
#' @param class_var The CDS cell metadata column holding sample class data.  There can be only 2 classes in this column.  You may need to subset or reclass the samples to achieve this.
#' @param experimental_class The value from the class column indicating the experimental group.
#' @param control_class The value from the class column indicating the control or reference class.
#' @param return_value Option to return either a plot or a data table for plotting in a separate step.  Must be either "plot" or "table".
#' @return A ggplot or a table of data for plotting
#' @export
#' @import tidyverse monocle3
#' @importFrom rstatix fisher_test
#' @importFrom cli::cli_alert_info
bb_cluster_representation <- function(cds,
                                      cluster_var,
                                      class_var,
                                      experimental_class,
                                      control_class,
                                      pseudocount = 1,
                                      return_value = c("table", "plot")) {
  cli::cli_alert_info("This function calculates differential representation on a per-cell basis and uses the Fisher exact test to calculate significance.")
  cli::cli_alert_info("In large datasets this may inflate the significance to an undesirable level.")
  cli::cli_alert_info("Alternatively consider using bb_cluster_representation2 which calculates these values using a regression method and may be more robust.")
  stopifnot("You must select a table or a plot to return." = return_value %in% c("table", "plot"))
  stopifnot("You can only return a table OR a plot." = length(return_value) == 1)
  res <- bb_cellmeta(cds) %>%
    group_by(!!sym(cluster_var), !!sym(class_var)) %>%
    summarise(n = n())
  all <- res |>
    ungroup() |>
    tidyr::expand(!!sym(cluster_var), !!sym(class_var))
  res <- right_join(res, all) |>
    mutate(n = replace_na(n, 0)) |>
    mutate(n = n + pseudocount)
  res <-
    left_join(res,
              bb_cellmeta(cds) %>%
                group_by(!!sym(class_var)) %>%
                summarise(class_total = n())) %>%
    mutate(overall_total = nrow(bb_cellmeta(cds))) %>%
    mutate(normalized_count = n * overall_total / class_total / 2) %>%
    select(!!sym(cluster_var),!!sym(class_var), normalized_count) %>%
    pivot_wider(names_from = !!(class_var), values_from = normalized_count) %>%
    mutate(fold_change_over_control = !!sym(experimental_class) / !!sym(control_class)) %>%
    mutate(log2fold_change_over_control = log2(!!sym(experimental_class) /
                                                 !!sym(control_class))) %>%
    mutate(enriched = ifelse(
      fold_change_over_control > 1,
      experimental_class,
      control_class
    ))
  # return(res)
  fisher <- map_dfr(
    .x = bb_cellmeta(cds) %>% pull(!!sym(cluster_var)) %>% unique() %>% as.character(),
    .f = function(x,
                  cds_int = cds,
                  cluster_var_int = cluster_var,
                  class_var_int = class_var,
                  pscount = pseudocount) {
      fish <- bb_cellmeta(cds_int) %>%
        mutate(cluster_of_interest = ifelse(!!as.name(cluster_var_int) %notin% x,
                                            paste0("not ", x),
                                            x)) %>%
        count(!!sym(class_var_int), cluster_of_interest) %>%
        mutate(n = n + pscount) |>
        pivot_wider(names_from = !!sym(class_var_int),
                    values_from = n, values_fill = pseudocount) %>%
        bb_tbl_to_matrix() |>
        rstatix::fisher_test() %>%
        mutate(!!cluster_var_int := x) %>%
        relocate(!!sym(cluster_var_int)) %>%
        dplyr::rename(fisher_exact_p = p)
      return(fish)
    }
  )
  res_final <- left_join(res, fisher) %>%
    mutate(texty = ifelse(
      log2(fold_change_over_control) > 0,
      log2(fold_change_over_control),
      0
    ))
  if (return_value == "table") {
    return(res_final)
  } else {
    plot <- ggplot(
      data = res_final,
      aes_string(x = cluster_var, y = "log2fold_change_over_control", fill = "enriched")
    ) +
      geom_bar(stat = "identity", color = "black") +
      labs(x = NULL, fill = "Enriched:", y = "l2FC") +
      geom_text(
        mapping = aes(y = texty, label = p.signif),
        size = 3,
        show.legend = F,
        vjust = -1.5

      )
    return(plot)
  }
}
