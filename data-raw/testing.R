devtools::load_all()
testthat::test_file("tests/testthat/test-SummarizedHeatmap.R")
library(patchwork)
mat <- matrix(rnorm(100), ncol=5)
colnames(mat) <- letters[1:5]
rownames(mat) <- letters[6:25]
test_sh <- SummarizedHeatmap(mat)
colData(test_sh)$sample_type <- c("vowel", "consonant", "consonant", "consonant", "vowel")
colData(test_sh)$sample_type2 <- c("vowel2", "consonant2", "consonant2", "consonant2", "vowel2")
isVowel <- function(char) char %in% c('a', 'e', 'i', 'o', 'u')
rowData(test_sh)$feature_type <- ifelse(isVowel(letters[6:25]), "vowel", "consonant")
rowData(test_sh)$feature_type2 <- paste0(rowData(test_sh)$feature_type, "2")


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
  ) * theme(plot.background = element_rect(color = "red"))
p4 <- bb_plot_heatmap_rowData(test_sh,
                          vars = c("Feature Type" = "feature_type",
                                   "Feature Type2" = "feature_type2"),
                          manual_pal = c("consonant" = "ivory3",
                                         "vowel" = "pink4",
                                         "consonant2" = "bisque1",
                                         "vowel2" = "lightpink3"))

# bb_plot_heatmap_main(test_sh, tile_color = "white") * scale_y_discrete(breaks = c("u", "m", "h", "v", "t"), position = "right", guide = guide_axis(n.dodge = 3))  +
# p5 <- bb_plot_heatmap_main(test_sh, tile_color = "white") * theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p5 <- bb_plot_heatmap_main(test_sh, tile_color = "white") + theme(plot.background = element_rect(color = "red"))

p6 <- guide_area()
design <- "
##1#
##36
2456

"
p1 + p2 + free(p3, side = "r",  type = "space") + free(p4, side = "t", type = "space") + p5 + p6 + plot_layout(
              design = design,
              guides = "collect")


library(patchwork)
devtools::load_all()
{
p1 <- bb_plot_heatmap_main(test_sh, flip = TRUE)
p2 <- bb_plot_heatmap_colDendro(test_sh, side = "left")
p3 <- bb_plot_heatmap_colData(test_sh, side = "left")
p4 <- bb_plot_heatmap_rowDendro(test_sh, side = "top")
p5 <- bb_plot_heatmap_rowData(test_sh, side = "top")
p6 <- guide_area()
p7 <- bb_plot_heatmap_colHighlight(test_sh, highlights = c("a", "b", "c"), side = "right")
p8 <- bb_plot_heatmap_rowHighlight(test_sh, highlights = c("w", "s", "v"), side = "bottom")

design <- "
##4#6
##5#6
23176
##8##
"
p1 + p2 + free(p3, side = "t", type = "space") + p4 + free(p5, side = "r", type = "space") + p6 + p7 + p8 + plot_layout(design = design, guides = "collect")
}


devtools::load_all()
{
p1 <- bb_plot_heatmap_main(test_sh, flip = FALSE) * theme(plot.background = element_rect(color = "red"))
p2 <- bb_plot_heatmap_colDendro(test_sh, side = "top")
p3 <- bb_plot_heatmap_colData(test_sh, side = "top") * theme(plot.background = element_rect(color = "red"))
p4 <- bb_plot_heatmap_rowDendro(test_sh, side = "left")
p5 <- bb_plot_heatmap_rowData(test_sh, side = "left") * theme(plot.background = element_rect(color = "red"))
p6 <- guide_area()
p7 <- bb_plot_heatmap_colHighlight(test_sh, highlights = c("a", "b", "c"), side = "bottom")
p8 <- bb_plot_heatmap_rowHighlight(test_sh, highlights = c("w", "s", "v"), side = "right")

design <- "
##2#6
##3#6
45186
##7##
"
p1 + p2 + free(p3, side = "r", type = "space") + p4 + free(p5, side = "t", type = "space") + p6 + p7 + p8 + plot_layout(design = design, guides = "collect")
}
