#' Calculate standard error of the mean 
#' 
#' @param x A numeric vector 
#' @export
se <- function(x) sqrt(var(x, na.rm=TRUE)/length(x))

#' Generate median +/- se stat object for jitter ggplot 
#' 
#' @param x A numeric vector 
#' @export
data_median_se <- function(x) {
  m <- median(x)
  ymin <- m - se(x)
  ymax <- m + se(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

#' Generate median +/- iqr stat object for jitter ggplot 
#' 
#' @param x A numeric vector 
#' @export
data_summary_median_iqr <- function(x) {
  m <- median(x)
  ymin <- m - (IQR(x) / 2)
  ymax <- m + (IQR(x) / 2)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

#' Generate mean +/- se stat object for jitter ggplot 
#' 
#' @param x A numeric vector 
#' @export
data_summary_mean_se <- function(x) {
  m <- mean(x)
  ymin <- m - se(x)
  ymax <- m + se(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

#' Generate median +/- mad stat object for jitter ggplot 
#' 
#' @param x A numeric vector 
#' @export
data_summary_median_mad <- function(x) {
  m <- median(x)
  ymin <- m - mad(x)
  ymax <- m + mad(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}

#' Generate mean +/- sd stat object for jitter ggplot 
#' 
#' @param x A numeric vector 
#' @export
data_summary_mean_sd <- function(x) {
  m <- mean(x)
  ymin <- m - sd(x)
  ymax <- m + sd(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}



#' Add predefined sample-level cell metadata to cell data sets while importing 
#' 
#' @param cds A cell data set object 
#' @param columns_to_add A named vector where the name of each element becomes the name of the new colData column and the value is the value for that particular sample.  Best used when importing from a metadata table. 
#' @export
add_cds_factor_columns <- function(cds, columns_to_add) {
  for (i in 1:length(columns_to_add)) {
    colData(cds)$new <- unname(columns_to_add[i])
    names(colData(cds))[names(colData(cds)) == "new"] <-
      names(columns_to_add[i])
  }
  return(cds)
}
