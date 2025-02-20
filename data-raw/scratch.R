devtools::load_all()
library(rtracklayer)
bb_import_bw(fs::path(
  "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/wantong/CUT&Tag/novogene_0904/igv_bw/K27Ac_R1.bigWig"
), group = "blah"
)

k27ac <-
  import.bw(
    fs::path(
      "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/wantong/CUT&Tag/novogene_0904/igv_bw/K27Ac_R1.bigWig"
    ),
  ) |>
  plyranges::mutate(group = "H3K27Ac") |>
  plyranges::mutate(coverage = score) |>
  plyranges::select(-score)

k27ac_peak <-
  bb_import_seacr_peaks(
    fs::path(
      "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/wantong/CUT&Tag/novogene_0904/igv_bw/K27Ac_R1.seacr.peaks.stringent.bed"
    ),
    group_variable = "group",
    group_value = "H3K27Ac"
  )


k27me3 <-
  import.bw(
    fs::path(
      "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/wantong/CUT&Tag/novogene_0904/igv_bw/K27me3_R1.bigWig"
    )
  ) |>
  plyranges::mutate(group = "H3K27me3") |>
  plyranges::mutate(coverage = score) |>
  plyranges::select(-score)

k27me3_peak <-
  bb_import_seacr_peaks(
    fs::path(
      "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/wantong/CUT&Tag/novogene_0904/igv_bw/K27me3_R1.seacr.peaks.stringent.bed"
    ),
    group_variable = "group",
    group_value = "H3K27me3"
  )

k4me3 <-
  import.bw(
    fs::path(
      "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/wantong/CUT&Tag/novogene_0904/igv_bw/K4me3_R2.bigWig"
    )
  ) |>
  plyranges::mutate(group = "H3K4me3") |>
  plyranges::mutate(coverage = score) |>
  plyranges::select(-score)


k4me3_peak <-
  bb_import_seacr_peaks(
    fs::path(
      "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/wantong/CUT&Tag/novogene_0904/igv_bw/K4me3_R2.seacr.peaks.stringent.bed"
    ),
    group_variable = "group",
    group_value = "H3K4me3"
  )


rel_bw <- import.bw(
    fs::path(
      "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/wantong/CUT&Tag/novogene_0904/igv_bw/Rel_R1.bigWig"
    )
  ) |>
  plyranges::mutate(group = "REL") |>
  plyranges::mutate(coverage = score) |>
  plyranges::select(-score)

rel_peak <-
  bb_import_seacr_peaks(
    fs::path(
      "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/wantong/CUT&Tag/novogene_0904/igv_bw/Rel_R1.seacr.peaks.stringent.bed"
    ),
    group_variable = "group",
    group_value = "REL"
  )

test_plotfun <- function(trace) {
  p1 <-
    bb_plot_trace_data(
      trace,
      pal = test_pal,
      group_filter = "H3K27Ac",
      group_variable = "group"
    ) +
    theme(strip.placement = "outside") +
    theme(axis.title.y.left = element_blank())
  p2 <-
    bb_plot_trace_peaks(
      trace,
      pal = test_pal,
      group_filter = "H3K27Ac",
      group_variable = "group"
    ) + theme_nothing()
  p3 <-
    bb_plot_trace_data(
      trace,
      pal = test_pal,
      group_filter = "H3K27me3",
      group_variable = "group"
    )+
    theme(strip.placement = "outside") +
    theme(axis.title.y.left = element_blank())
  p4 <-
    bb_plot_trace_peaks(
      trace,
      pal = test_pal,
      group_filter = "H3K27me3",
      group_variable = "group"
    ) + theme_nothing()
  p5 <-
    bb_plot_trace_data(
      trace,
      pal = test_pal,
      group_filter = "H3K4me3",
      group_variable = "group"
    )+
    theme(strip.placement = "outside") +
    theme(axis.title.y.left = element_blank())
  p6 <-
    bb_plot_trace_peaks(
      trace,
      pal = test_pal,
      group_filter = "H3K4me3",
      group_variable = "group"
    ) + theme_nothing()
  p7 <-
    bb_plot_trace_data(
      trace,
      pal = test_pal,
      group_filter = "REL",
      group_variable = "group"
    ) +
    theme(strip.placement = "outside") +
    theme(axis.title.y.left = element_blank())
  p8 <-
    bb_plot_trace_peaks(
      trace,
      pal = test_pal,
      group_filter = "REL",
      group_variable = "group"
    ) + theme_nothing()

  p9 <-
    bb_plot_trace_model(trace) + theme(axis.line.y = element_blank()) + labs(y = NULL)
  p10 <- bb_plot_trace_axis(trace)

  p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + plot_layout(heights = c(3,
                                                                  0.3,
                                                                  3,
                                                                  0.3,
                                                                  3,
                                                                  0.3,
                                                                  3,
                                                                  0.3,
                                                                  1,
                                                                  0.1))
}


test_pal <- c("H3K27Ac" = "green",
              "H3K27me3" = "red",
              "H3K4me3" = "blue",
              "REL" = "black")

test_dll4_trace <- bb_makeTrace(
  c(k27ac, k27me3, k4me3, rel_bw),
  gene_to_plot = "DLL4",
  genome = "hg38",
  extend_left = 10000,
  extend_right = 10000,
  fill_in = TRUE,
  peaks = c(k27ac_peak, k27me3_peak, k4me3_peak, rel_peak)
)

test_irf4_trace <- bb_makeTrace(
  c(k27ac, k27me3, k4me3, rel_bw),
  gene_to_plot = "IRF4",
  genome = "hg38",
  extend_left = 10000,
  extend_right = 10000,
  fill_in = TRUE,
  peaks = c(k27ac_peak, k27me3_peak, k4me3_peak, rel_peak)
)

test_EP300_trace <- bb_makeTrace(
  c(k27ac, k27me3, k4me3, rel_bw),
  gene_to_plot = "EP300",
  genome = "hg38",
  extend_left = 10000,
  extend_right = 10000,
  fill_in = TRUE,
  peaks = c(k27ac_peak, k27me3_peak, k4me3_peak, rel_peak)
)

test_dll4_trace@trace_data
test_irf4_trace@trace_data
test_EP300_trace@trace_data
test_dll4_trace@peaks
test_plotfun(trace = test_dll4_trace)
test_plotfun(trace = test_irf4_trace)
test_plotfun(trace = test_EP300_trace)


# load("~/brad_workspace/huvec.multiome.datapkg.shared_data/signac_main.rda")
# signac_main <- bb_fragment_replacement(obj = signac_main,
#                                              new_path = fs::path(c("~/network/X/Labs/Blaser/staff/single_cell/multiome_march2023/pipestance/B1_C1/atac_fragments.tsv.gz",
#                                                                    "~/network/X/Labs/Blaser/staff/single_cell/multiome_jan2024/pipestance/batch2_sample1/atac_fragments.tsv.gz")))
test_dll4_trace_signac <- bb_makeTrace(
  signac_main,
  gene_to_plot = "DLL4",
  genome = "hg38",
  extend_left = 10000,
  extend_right = 10000, fill_in = TRUE
)
test_dll4_trace_signac@trace_data




