library(SummarizedExperiment)
library(ggdendro)
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
test_that("Summarized Heatmap works", {
  expect_true(validObject(test_sh))
  expect_identical(rownames(colData(test_sh)), letters[1:5])
  expect_identical(rownames(rowData(test_sh)), letters[6:25])
  expect_identical(attributes(bb_plot_heatmap_main(test_sh))$class, c("gg", "ggplot"))
  expect_identical(attributes(bb_plot_heatmap_rowDendro(test_sh))$class, c("gg", "ggplot"))
  expect_identical(attributes(bb_plot_heatmap_colDendro(test_sh))$class, c("gg", "ggplot"))
  expect_identical(colData(test_sh), colData_test)
  expect_identical(rowData(test_sh), rowData_test)
})
