devtools::load_all()
testthat::test_file("tests/testthat/test-SummarizedHeatmap.R")
library(patchwork)
sample_data <- tibble(
  the_letters = letters[1:5],
  sample_type = c("vowel", "consonant", "consonant", "consonant", "vowel"),
  sample_type2 = c("vowel2", "consonant2", "consonant2", "consonant2", "vowel2")
)
isVowel <- function(char) char %in% c('a', 'e', 'i', 'o', 'u')
feature_data <- tibble(
  the_letters = letters[6:25],
  feature_type = ifelse(isVowel(the_letters), "vowel", "consonant"),
  feature_type2  = paste0(feature_type, "2")
)

# generate a sample matrix
mat <- matrix(rnorm(100), ncol = 5)
colnames(mat) <- sample_data$the_letters
rownames(mat) <- feature_data$the_letters
mat

test_sh <- SummarizedHeatmap(mat)

colData(test_sh) <-
  DataFrame(sample_data, row.names = sample_data$the_letters)
rowData(test_sh) <- DataFrame(feature_data, row.names = feature_data$the_letters)
colData(test_sh)
rowData(test_sh)
ggdendro::dendro_data(test_sh@colDendro)

# plot showing everything
{
  p1 <- bb_plot_heatmap_main(test_sh, flip = FALSE)
  p2 <- bb_plot_heatmap_colDendro(test_sh, side = "top")
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

# flipped
{
  p1 <- bb_plot_heatmap_main(test_sh, flip = TRUE)
  p2 <- bb_plot_heatmap_colDendro(test_sh, side = "left")
  p3 <- bb_plot_heatmap_colData(test_sh, side = "left")
  p4 <- bb_plot_heatmap_rowDendro(test_sh, side = "top")
  p5 <- bb_plot_heatmap_rowData(test_sh, side = "top")
  p6 <- guide_area()
  p7 <-
    bb_plot_heatmap_colHighlight(test_sh,
                                 highlights = c("a", "b", "c"),
                                 side = "right")
  p8 <-
    bb_plot_heatmap_rowHighlight(test_sh,
                                 highlights = c("w", "s", "v"),
                                 side = "bottom")

  design <- "
##4#6
##5#6
23176
##8##
"
  p1 + p2 + free(p3, side = "t", type = "space") + p4 + free(p5, side = "r", type = "space") + p6 + p7 + p8 + plot_layout(design = design, guides = "collect")
}

# some customizations
{
  p1 <- bb_plot_heatmap_colDendro(test_sh)
  p2 <- bb_plot_heatmap_rowDendro(test_sh, side = "left")
  p3 <- bb_plot_heatmap_colData(
    test_sh,
    vars = c("Sample Type" = "sample_type", "Sample Type 2" = "sample_type2"),
    manual_pal = c(
      "consonant" = "ivory",
      "vowel" = "pink",
      "consonant2" = "cornsilk",
      "vowel2" = "pink3"
    )
  )
  p4 <- bb_plot_heatmap_rowData(
    test_sh,
    vars = c("Feature Type" = "feature_type",
             "Feature Type2" = "feature_type2"),
    manual_pal = c(
      "consonant" = "ivory3",
      "vowel" = "pink4",
      "consonant2" = "bisque1",
      "vowel2" = "lightpink3"
    )
  )


  p5 <- bb_plot_heatmap_main(test_sh, tile_color = "white")

  p6 <- guide_area()
  design <- "
##1#
##36
2456

"
  p1 + p2 + free(p3, side = "r",  type = "space") + free(p4, side = "t", type = "space") + p5 + p6 + plot_layout(design = design, widths = c(.1, .2, 1, 0.5), heights = c(.1, .2, 1),
                                                                                                                 guides = "collect")

}

test_sh1 <- SummarizedHeatmap(mat, colOrder = letters[1:5])
colData(test_sh1) <-
  DataFrame(sample_data, row.names = sample_data$the_letters)
rowData(test_sh1) <- DataFrame(feature_data, row.names = feature_data$the_letters)
# set col order
{
  p1 <- bb_plot_heatmap_main(test_sh1)
  p2 <- patchwork::plot_spacer()
  p3 <- bb_plot_heatmap_colData(test_sh1)
  p4 <- bb_plot_heatmap_rowDendro(test_sh1)
  p5 <- bb_plot_heatmap_rowData(test_sh1)
  p6 <- guide_area()
  p7 <-
    bb_plot_heatmap_colHighlight(test_sh1,
                                 highlights = c("a", "b", "c"),
                                 side = "bottom")
  p8 <-
    bb_plot_heatmap_rowHighlight(test_sh1,
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
colOrder(test_sh1)
