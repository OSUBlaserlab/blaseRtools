suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("tidyverse"))


# need:  hg38_granges, Hsapiens


# load the zfin gene model-------------------------------------------
hg38_gff <- read_delim(file = "~/workspace_pipelines/sc_refdata/ensembl_hg38_20211126/Homo_sapiens.GRCh38.104.chr.gff3.gz", delim = "\t", comment = "#", col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")) %>%
  mutate(ENS = str_extract(attributes, pattern = "ENS[:alpha:][:digit:]*")) %>%
  select(-attributes) %>%
  mutate(lookup_var = str_extract(ENS, "ENS[:alpha:]")) %>%
  filter(!is.na(lookup_var)) %>%
  mutate(lookup_var = case_when(lookup_var == "ENSG" ~ "gene_stable_ID",
                                lookup_var == "ENSP" ~ "protein_stable_ID",
                                lookup_var == "ENST" ~ "transcript_stable_ID"))

hg38_megalookup <- read_tsv("~/workspace_pipelines/sc_refdata/ensembl_hg38_20211126/hg38_megalookup.txt") %>%
  select(gene_stable_ID = `Gene stable ID`, transcript_stable_ID = `Transcript stable ID`, protein_stable_ID = `Protein stable ID`, gene_name = `Gene name`)

hg38_granges_df <- map_dfr(.x = c("gene_stable_ID", "protein_stable_ID", "transcript_stable_ID"),
        .f = function (x, data = hg38_gff, lookup = hg38_megalookup) {
          lookup <- pivot_longer(lookup, cols = c(gene_stable_ID, transcript_stable_ID, protein_stable_ID), values_to = "ENS") %>%
            filter(name == x) %>%
            group_by(gene_name, ENS) %>%
            summarise()
          data <- data %>% filter(lookup_var == x) %>% select(-lookup_var)
          res <- left_join(data, lookup)
          return(res)

        }) %>%
  arrange(gene_name) %>%
  select(-c(score, phase)) %>%
  mutate(seqid = paste0("chr", seqid))


hg38_exons_numbered <- bind_rows(
  hg38_granges_df %>%
    filter(type == "exon") %>%
    filter(strand == "+") %>%
    group_by(gene_name) %>%
    arrange(start, .by_group = T) %>%
    mutate(exon_number = paste0("exon_", row_number())),

  hg38_granges_df %>%
    filter(type == "exon") %>%
    filter(strand == "-") %>%
    group_by(gene_name) %>%
    arrange(desc(start), .by_group = T) %>%
    mutate(exon_number = paste0("exon_", row_number()))
)

hg38_granges_df1 <- left_join(hg38_granges_df, hg38_exons_numbered) %>%
  mutate(label = ifelse(type == "exon", exon_number, type)) %>%
  select(-c(exon_number, source, ENS))
hg38_granges_df1

# the granges object; made from zfin gene model
hg38_granges <- makeGRangesFromDataFrame(hg38_granges_df1, keep.extra.columns = T)

Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
