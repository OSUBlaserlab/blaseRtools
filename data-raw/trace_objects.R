library(tidyverse)
library(GenomicRanges)
hg38_gff <- read_delim(file = "~/workspace_pipelines/sc_refdata/ensembl_hg38_20211126/Homo_sapiens.GRCh38.104.chr.gff3.gz",
                       delim = "\t",
                       comment = "#",
                       col_names = c("seqid",
                                     "source",
                                     "type",
                                     "start",
                                     "end",
                                     "score",
                                     "strand",
                                     "phase",
                                     "attributes"))


hg38_megabiglookup <- read_delim(file = "~/workspace_pipelines/sc_refdata/ensembl_hg38_20211126/hg38_megabiglookup.txt", delim = "\t")

hg38_full_model <-
  hg38_gff %>%
  filter(type %in% c("five_prime_UTR", "CDS", "three_prime_UTR")) %>%
  mutate(parent_transcript = str_extract(attributes, "Parent=transcript:ENST[:digit:]*")) %>%
  select(-attributes) %>%
  mutate(parent_transcript = str_replace(parent_transcript, "Parent=transcript:", "")) %>%
  left_join(
    hg38_megabiglookup %>% select(
      parent_transcript = `Transcript stable ID`,
      APPRIS = `APPRIS annotation`,
      length = `Transcript length (including UTRs and CDS)`,
      gene_name = `Gene name`
    )
  ) %>% mutate(APPRIS = factor(
    APPRIS,
    levels = c(
      "principal1",
      "principal2",
      "principal3",
      "principal4",
      "principal5",
      "alternative1",
      "alternative2"
    )
  )) %>%
  mutate(seqid = paste0("chr", seqid)) %>%
  select(-c(score, phase))

hg38_full_model_gr <- makeGRangesFromDataFrame(hg38_full_model, keep.extra.columns = T)
genome(hg38_full_model_gr) <- "hg38"


dr11_gff <- read_delim(file = "~/workspace_pipelines/sc_refdata/zfin_20211111/zfin_genes.gff3",
                       delim = "\t",
                       comment = "#",
                       col_names = c("seqid",
                                     "source",
                                     "type",
                                     "start",
                                     "end",
                                     "score",
                                     "strand",
                                     "phase",
                                     "attributes"))
dr11_megabiglookup <- read_delim("~/workspace_pipelines/sc_refdata/zfin_20211111/dr11_megabiglookup.txt", delim = "\t")


dr11_full_model <-
  dr11_gff %>%
  filter(type %in% c("five_prime_UTR", "CDS", "three_prime_UTR")) %>%
  mutate(parent_transcript = str_extract(attributes, "Parent=ENSDART[:digit:]*")) %>%
  select(-attributes) %>%
  mutate(parent_transcript = str_replace(parent_transcript, "Parent=", "")) %>%
  left_join(
    dr11_megabiglookup %>% select(
      parent_transcript = `Transcript stable ID`,
      APPRIS = `APPRIS annotation`,
      length = `Transcript length (including UTRs and CDS)`,
      gene_name = `Gene name`
    )
  ) %>% mutate(APPRIS = factor(
    APPRIS,
    levels = c(
      "principal1",
      "principal5",
      "alternative1",
      "alternative2"
    )
  )) %>%
  mutate(seqid = paste0("chr", seqid)) %>%
  select(-c(phase, score))

dr11_full_model_gr <- makeGRangesFromDataFrame(dr11_full_model, keep.extra.columns = T)
genome(dr11_full_model_gr) <- "danRer11"

