library("blaseRtools")
library("tidyverse")
prkcda_ape <- bb_grcz11_ape("prkcda")
out <- NULL
FEATURES(prkcda_ape)
fimo_feature <- "prkcda_exon_10"
# make jaspar object
# convert jaspar sequences to meme format
cmd <-
  paste0(
    "jaspar2meme -bundle /workspace/wantong_workspace/pkc.cxcl8.datapkg/inst/extdata/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt > meme"
  )
message(cmd, "\n")
system(cmd)

meme <- read_lines("meme")
fastSave::save.pigz(meme, file = "brad/meme.rda")

