#' Print out a stats report 
#' 
#' @param data A Tibble in tidy data format 
#' @param class_var Column containing the class variable
#' @param class1 The first class to compare
#' @param class2 The second class to compare
#' @param value_var The column containing the values 
#' @param out The text file to output
#' @return A text file 
#' @export
#' @import tidyverse
#' @importFrom huxtable hux add_colnames set_bold set_all_borders 
#' @examples
bb_print_full_stats <- function(data, 
				class_var, 
				class1, 
				class2, 
				value_var, 
				out = NULL) {
	if (!is.null(out)) {
	  sink(out) 
	}
	data %>%
	      group_by(!!sym(class_var)) %>%
	      summarise(n = n(),
			mean = mean(!!sym(value_var)),
			sd = sd(!!sym(value_var)),
			se = se(!!sym(value_var)),
			median = median(!!sym(value_var)),
			IQR = IQR(!!sym(value_var)),
			shapiro_p = shapiro.test(!!sym(value_var))[[2]]) %>%
	      huxtable::hux() %>%
	      huxtable::add_colnames() %>%
	      huxtable::set_bold(row = 1, col = everywhere, value = TRUE) %>% 
	      huxtable::set_all_borders(TRUE)
	class1_data <- data %>%
		filter(!!sym(class_var) == class1) %>%
		pull(!!sym(value_var))
	class2_data <- data %>%
		filter(!!sym(class_var) == class2) %>%
		pull(!!sym(value_var))
	cat("_________________________________")
	cat("\n\n")
	cat(paste0("Class 1 is ",class1," and Class 2 is ",class2,"\n"))
	cat("_________________________________")
	cat("\n")
	print(shapiro.test(class1_data))
	cat("_________________________________")
	cat("\n")
	print(shapiro.test(class2_data))
	cat("_________________________________")
	cat("\n")
	print(t.test(class1_data, class2_data, alternative = "two.sided", var.equal = T))
	cat("_________________________________")
	cat("\n")
	print(t.test(class1_data, class2_data, alternative = "two.sided", var.equal = F))
	cat("_________________________________")
	cat("\n")
	print(wilcox.test(class1_data, class2_data))
	if (!is.null(out)){
	  sink()
	}
}
