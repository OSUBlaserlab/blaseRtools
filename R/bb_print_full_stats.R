#' Print out a stats report
#'
#' @param data A Tibble in tidy data format.  Must contain or be filtered to contain only 2 levels in "classification_variable" for comparisons.
#' @param classification_variable Column containing the class variable
#' @param numeric_variable The column containing the numeric values to summarize and compare
#' @param test_type Must be one of "Student", "Welch", and "Wilcox"
#' @param output Output file; if null prints to screen.
#' @return A text file
#' @export
#' @import tidyverse
bb_print_full_stats <- function(data,
                                classification_variable,
                                numeric_variable,
                                test_type = c("Student", "Welch", "Wilcox"),
                                output = NULL) {
  if (test_type == "Student") {
    if (!is.null(output))
      sink(file = output)
    data %>%
      group_by(!!as.name(classification_variable)) %>%
      rstatix::get_summary_stats(!!as.name(numeric_variable), type = "common") %>%
      huxtable::as_hux() %>%
      huxtable::set_all_borders(TRUE) %>%
      huxtable::set_caption("Summary Stats") %>%
      huxtable::print_screen()
    cat("\n\n")
    data %>%
      rstatix::t_test(as.formula(
        str_glue("{numeric_variable} ~ {classification_variable}")
      ), var.equal = TRUE) %>%
      huxtable::as_hux() %>%
      huxtable::set_all_borders(TRUE) %>%
      huxtable::set_caption("Student's T Test") %>%
      huxtable::print_screen()
    if (!is.null(output))
      sink()

  }
  if (test_type == "Welch") {
    if (!is.null(output))
      sink(file = output)
    data %>%
      group_by(!!as.name(classification_variable)) %>%
      rstatix::get_summary_stats(!!as.name(numeric_variable), type = "common") %>%
      huxtable::as_hux() %>%
      huxtable::set_all_borders(TRUE) %>%
      huxtable::set_caption("Summary Stats") %>%
      huxtable::print_screen()
    cat("\n\n")
    data %>%
      rstatix::t_test(as.formula(
        str_glue("{numeric_variable} ~ {classification_variable}")
      ), var.equal = FALSE) %>%
      huxtable::as_hux() %>%
      huxtable::set_all_borders(TRUE) %>%
      huxtable::set_caption("Welch's T Test") %>%
      huxtable::print_screen()
    if (!is.null(output))
      sink()
  }
  if (test_type == "Wilcox") {
    if (!is.null(output))
      sink(file = output)
    data %>%
      group_by(!!as.name(classification_variable)) %>%
      rstatix::get_summary_stats(!!as.name(numeric_variable), type = "common") %>%
      huxtable::as_hux() %>%
      huxtable::set_all_borders(TRUE) %>%
      huxtable::set_caption("Summary Stats") %>%
      huxtable::print_screen()
    cat("\n\n")
    data %>%
      rstatix::wilcox_test(as.formula(
        str_glue("{numeric_variable} ~ {classification_variable}")
      )) %>%
      huxtable::as_hux() %>%
      huxtable::set_all_borders(TRUE) %>%
      huxtable::set_caption("Wilcoxon Rank Sum Test") %>%
      huxtable::print_screen()
    if (!is.null(output))
      sink()

  }
}
