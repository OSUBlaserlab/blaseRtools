#' A function to run qc tests on cds objects.
#'
#' @param cds A cell data set object to run qc functions on
#' @param cds_name The name of the cds
#' @param genome The species to use for identifying mitochondrial genes. Choose from "human", "mouse", "zfish", "human_mouse" for pdx.
#' @param max_mito Manual cutoff for mitochondrial percentage.  May be more strict, i.e. lower, than the automated cutoff but not less strict,  Default: NULL
#' @param min_log_detected Manual cutoff for log detected features.  May be more strict, i.e. higher, than the automated cutoff not not less strict, Default: NULL
#' @return A list of qc objects
#' @export
#' @import tidyverse scater
bb_qc <-
  function (cds,
            cds_name,
            genome = c("human", "mouse", "zfish"),
            max_mito = NULL,
            min_log_detected = NULL) {
    if (genome == "human") {
      mito_pattern <- "^MT-"
    }
    if (genome == "mouse") {
      mito_pattern <- "^mt-"
    }
    if (genome == "zfish") {
      mito_pattern <- "^mt-"
    }
    if (genome == "human_mouse") {
      mito_pattern <- "^MT-|^mt-"
    }
    cds <-
      scater::addPerCellQC(cds, subsets = list(Mito = grep(
        mito_pattern,
        rowData(cds)$gene_short_name
      )))
    cds_tbl <-
      tibble::as_tibble(colData(cds)) %>%
      dplyr::mutate(
        qc.detected = scater::isOutlier(
          detected,
          log = TRUE,
          type = "lower",
          nmads = 2
        ),
        qc.mito = scater::isOutlier(subsets_Mito_percent,
                            type = "higher", nmads = 2),
        qc.any = qc.detected | qc.mito,
        pre_post = "pre"
      )
    if (!is.null(max_mito)) {
      cds_tbl <- cds_tbl |>
        mutate(max_mito_col = subsets_Mito_percent > max_mito) |>
        mutate(qc.mito = qc.mito | max_mito_col) |>
        select(-max_mito_col)
    }
    if (!is.null(min_log_detected)) {
      cds_tbl <- cds_tbl |>
        mutate(min_detected_col = detected < min_log_detected) |>
        mutate(qc.detected = qc.detected | min_log_detected) |>
        select(-min_log_detected)
    }

    cds_tbl <- cds_tbl |>
      mutate(qc.any = qc.detected | qc.mito,
             pre_post = "pre")



    qc_detected_thresh <-
      attr(scater::isOutlier(
        cds_tbl$detected,
        log = TRUE,
        type = "lower",
        nmads = 2
      ),
      "thresholds")
    qc_mito_thresh <- attr(scater::isOutlier(
      cds_tbl$subsets_Mito_percent,
      type = "higher",
      nmads = 2
    ),
    "thresholds")
    cds_tbl_filtered_any <- cds_tbl %>%
      dplyr::filter(qc.any == FALSE) %>%
      dplyr::mutate(pre_post = "post")
    cds_tbl_filtered_detected <- cds_tbl %>%
      dplyr::filter(qc.detected == FALSE) %>%
      dplyr::mutate(pre_post = "post")
    cds_tbl_filtered_mito <- cds_tbl %>%
      dplyr::filter(qc.mito == FALSE) %>%
      dplyr::mutate(pre_post = "post")
    cds_to_plot_any <- dplyr::bind_rows(cds_tbl, cds_tbl_filtered_any)
    cds_to_plot_detected <-
      dplyr::bind_rows(cds_tbl, cds_tbl_filtered_detected)
    cds_to_plot_mito <- dplyr::bind_rows(cds_tbl, cds_tbl_filtered_mito)
    plot_any <- ggplot2::ggplot(cds_to_plot_any,
                                ggplot2::aes(x = forcats::fct_rev(pre_post),
                                            y = log10(detected))) +
      ggplot2::geom_violin() +
      ggplot2::labs(y = "log10(detected features)",
                    title = paste0("qc.any:  ", cds_name),
                    x = "thresholding") +
      ggplot2::theme(plot.title = element_text(hjust = 0.5))
    plot_mito <- ggplot2::ggplot(cds_to_plot_mito, ggplot2::aes(x = forcats::fct_rev(pre_post),
                                              y = subsets_Mito_percent)) +
      ggplot2::geom_violin() +
      ggplot2::geom_hline(yintercept = qc_mito_thresh[[2]]) +
      ggplot2::scale_y_continuous(breaks = c(qc_mito_thresh[[2]], seq(
        0,
        max(cds_to_plot_mito$subsets_Mito_percent), by = 20
      ))) +
      ggplot2::labs(y = "Percent Mitochondrial",
           title = paste0("qc.mito:  ",
                          cds_name),
           x = "thresholding") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    plot_detected <-
      ggplot2::ggplot(cds_to_plot_detected,
                      ggplot2::aes(x = forcats::fct_rev(pre_post),
                                       y = log10(detected))) +
      ggplot2::geom_violin() +
      ggplot2::geom_hline(yintercept = log10(qc_detected_thresh[[1]])) +
      ggplot2::scale_y_continuous(breaks = c(log10(qc_detected_thresh[[1]]),
                                    seq(0, max(
                                      log10(cds_to_plot_detected$detected)
                                    ),
                                    by = 1))) + ggplot2::labs(
                                      y = "Log10(detected features)",
                                      title = paste0("qc.detected:  ", cds_name),
                                      x = "thresholding"
                                    ) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    cds_return <- cds_tbl %>%
      dplyr::select(barcode, qc.any)
    return_list <-
      list(
        cds_return,
        plot_any,
        plot_mito,
        plot_detected,
        qc_detected_thresh,
        qc_mito_thresh
      )
    return(return_list)
  }



