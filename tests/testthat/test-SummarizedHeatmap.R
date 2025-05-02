library(SummarizedExperiment)
library(ggdendro)
library(patchwork)
set.seed(123)
mat <- matrix(rnorm(100), ncol=5)
colnames(mat) <- letters[1:5]
rownames(mat) <- letters[6:25]
test_sh <- SummarizedHeatmap(mat)
colData(test_sh)$sample_type <- c("vowel", "consonant", "consonant", "consonant", "vowel")
colData_test <- S4Vectors::DataFrame(
          sample_type = c("vowel", "consonant", "consonant", "consonant", "vowel"),
          row.names = letters[1:5])
isVowel <- function(char) char %in% c('a', 'e', 'i', 'o', 'u')
rowData(test_sh)$feature_type <- ifelse(isVowel(letters[6:25]), "vowel", "consonant")
rowData_test <- S4Vectors::DataFrame(
                          feature_type = ifelse(isVowel(letters[6:25]), "vowel", "consonant"),
                          row.names = letters[6:25])
test_sh1 <- SummarizedHeatmap(mat, colOrder = letters[1:5])
test_sh2 <- SummarizedHeatmap(mat, rowOrder = letters[6:25])

test_heatmap <- function() {
  p1 <- bb_plot_heatmap_main(test_sh, flip = FALSE)
  p2 <- bb_plot_heatmap_colDendro(test_sh, side = "bottom")
  p3 <- bb_plot_heatmap_colData(test_sh, side = "top")
  p4 <- bb_plot_heatmap_rowDendro(test_sh, side = "left")
  p5 <- bb_plot_heatmap_rowData(test_sh, side = "left")
  p6 <- guide_area()
  p7 <-
    bb_plot_heatmap_colHighlight(test_sh,
                                 highlights = c("a", "b", "c"),
                                 side = "bottom")
  p8 <-
    bb_plot_heatmap_rowHighlight(test_sh,
                                 highlights = c("w", "s", "v"),
                                 side = "right")

  design <- "
##2#6
##3#6
45186
##7##
"
p1 + p2 + free(p3, side = "r", type = "space") + p4 + free(p5, side = "t", type = "space") + p6 + p7 + p8 + plot_layout(design = design, guides = "collect")
}

test_heatmap()

test_that("Summarized Heatmap works", {
  expect_true(validObject(test_sh))
  expect_identical(rownames(colData(test_sh)), letters[1:5])
  expect_identical(rownames(rowData(test_sh)), letters[6:25])
  expect_identical(attributes(bb_plot_heatmap_main(test_sh))$class, c("gg", "ggplot"))
  expect_identical(attributes(bb_plot_heatmap_rowDendro(test_sh))$class, c("gg", "ggplot"))
  expect_identical(attributes(bb_plot_heatmap_colDendro(test_sh))$class, c("gg", "ggplot"))
  expect_identical(colData(test_sh), colData_test)
  expect_identical(rowData(test_sh), rowData_test)
  expect_true(validObject(test_sh1))
  expect_null(colDendro(test_sh1))
  expect_true(validObject(test_sh2))
  expect_null(rowDendro(test_sh2))
  expect_error(bb_plot_heatmap_rowDendro(test_sh2))
  expect_error(bb_plot_heatmap_colDendro(test_sh1))
})
