library(blaseRtools)
library(tidyverse)
prkcda_ape <- bb_grcz11_ape("prkcda")
Ape.DNA(prkcda_ape)
ape_gr <- Ape.granges(prkcda_ape)
Ape.setFeatures(prkcda_ape, ape_gr)
ape_fimo <- Ape.fimo(prkcda_ape, c("prkcda_exon_1", "prkcda_exon_2"))
ape_fimo
prkcda_ape2 <- Ape.setFeatures(prkcda_ape, c(ape_gr, ape_fimo))
Ape.save(prkcda_ape2, out = "~/network/X/Labs/Blaser/Brad/prkcda_ape2.ape")

Ape.fimo(prkcda_ape, fimo_feature = "prkcda_exon_1", out = "/workspace/workspace_pipelines/test/test_fimo")
