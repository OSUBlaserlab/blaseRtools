#' A function to run qc tests on cds objects. 
#' 
#' @param cds A cell data set object to run qc functions on
#' @param cds_name The name of the cds   
#' @param genome The species to use for identifying mitochondrial genes.
#' @return A list of qc objects
#' @export
#' @import tidyverse scater
#' @examples
bb_qc <- function(cds, 
		  cds_name, 
		  genome = c("human", 
			     "mouse", 
			     "zfish")) {
  if (genome == "human") {
    mito_pattern <-"^MT-"
  }
  if (genome == "mouse") {
    mito_pattern <- "^mt-"
  }
  if (genome == "zfish") {
    mito_pattern <- "^mt-"
  }
  cds <-
    scater::addPerCellQC(cds,subsets=list(Mito=grep(mito_pattern, rowData(cds)$gene_short_name)))
  cds_tbl <- as_tibble(colData(cds)) %>%
    mutate(
      qc.detected = isOutlier(
        detected,
        log = TRUE,
        type = "lower",
        nmads = 2
      ),
      qc.mito = isOutlier(
        subsets_Mito_percent,
        type = "higher",
        nmads = 2
      ),
      qc.any = qc.detected | qc.mito,
      pre_post = "pre"
    )
  qc_detected_thresh <-
    attr(
      isOutlier(
        cds_tbl$detected,
        log = TRUE,
        type = "lower",
        nmads = 2
      ),
      "thresholds"
    )
  qc_mito_thresh <-
    attr(isOutlier(cds_tbl$subsets_Mito_percent, type = "higher", nmads = 2),
         "thresholds")
  cds_tbl_filtered_any <-
    cds_tbl %>% filter(qc.any == FALSE) %>% mutate(pre_post = "post")
  cds_tbl_filtered_detected <-
    cds_tbl %>% filter(qc.detected == FALSE) %>% mutate(pre_post = "post")
  cds_tbl_filtered_mito <-
    cds_tbl %>% filter(qc.mito == FALSE) %>% mutate(pre_post = "post")
  cds_to_plot_any <- bind_rows(cds_tbl, cds_tbl_filtered_any)
  cds_to_plot_detected <-
    bind_rows(cds_tbl, cds_tbl_filtered_detected)
  cds_to_plot_mito <- bind_rows(cds_tbl, cds_tbl_filtered_mito)
  plot_any <-
    ggplot(cds_to_plot_any, aes(x = fct_rev(pre_post), y = log10(detected))) +
    geom_violin() + 
    labs(y = "log10(detected features)", title = paste0("qc.any:  ",cds_name), x = "thresholding") +
    theme(plot.title = element_text(hjust = 0.5))
  plot_mito <-
    ggplot(cds_to_plot_mito, aes(x = fct_rev(pre_post), y = subsets_Mito_percent)) +
    geom_violin() +
    geom_hline(yintercept = qc_mito_thresh[[2]]) +
    scale_y_continuous(breaks = c(qc_mito_thresh[[2]], seq(
      0, max(cds_to_plot_mito$subsets_Mito_percent), by = 20
    ))) +
    labs(y = "Percent Mitochondrial", title = paste0("qc.mito:  ",cds_name), x = "thresholding") +
    theme(plot.title = element_text(hjust = 0.5))
  plot_detected <-
    ggplot(cds_to_plot_detected, aes(
      x = fct_rev(pre_post),
      y = log10(detected)
    )) +
    geom_violin() +
    geom_hline(yintercept = log10(qc_detected_thresh[[1]])) +
    scale_y_continuous(breaks = c(log10(qc_detected_thresh[[1]]), seq(0, max(
      log10(cds_to_plot_detected$detected)
    ), by = 1))) +
    labs(y = "Log10(detected features)", title = paste0("qc.detected:  ",cds_name), x = "thresholding") +
    theme(plot.title = element_text(hjust = 0.5))
  cds_return <- cds_tbl %>% select(barcode, qc.any)
  return_list <-
    list(cds_return,
         plot_any,
         plot_mito,
         plot_detected,
         qc_detected_thresh,
         qc_mito_thresh)
  return(return_list)
}
