suppressPackageStartupMessages(library(BSgenome.Drerio.UCSC.danRer11))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(tidyverse))


# load the zfin gene model-------------------------------------------
zfin_gff <- read_delim(file = "~/workspace_pipelines/sc_refdata/zfin_20211111/zfin_genes.gff3", delim = "\t", skip = 35, col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")) %>%
  mutate(ENSDAR = str_extract(attributes, pattern = "ENSDAR[:alpha:][:digit:]*")) %>%
  select(-attributes) %>%
  mutate(lookup_var = str_extract(ENSDAR, "ENSDAR[:alpha:]")) %>%
  filter(!is.na(lookup_var)) %>%
  mutate(lookup_var = case_when(lookup_var == "ENSDARG" ~ "gene_stable_ID",
                                lookup_var == "ENSDARP" ~ "protein_stable_ID",
                                lookup_var == "ENSDART" ~ "transcript_stable_ID"))

zfin_megalookup <- read_tsv("~/workspace_pipelines/sc_refdata/zfin_20211111/zfin_megalookup.txt") %>%
  select(gene_stable_ID = `Gene stable ID`, transcript_stable_ID = `Transcript stable ID`, protein_stable_ID = `Protein stable ID`, gene_name = `Gene name`)

zfin_granges_df <- map_dfr(.x = c("gene_stable_ID", "protein_stable_ID", "transcript_stable_ID"),
        .f = function (x, data = zfin_gff, lookup = zfin_megalookup) {
          lookup <- pivot_longer(lookup, cols = c(gene_stable_ID, transcript_stable_ID, protein_stable_ID), values_to = "ENSDAR") %>%
            filter(name == x) %>%
            group_by(gene_name, ENSDAR) %>%
            summarise()
          data <- data %>% filter(lookup_var == x) %>% select(-lookup_var)
          res <- left_join(data, lookup)
          return(res)

        }) %>%
  arrange(gene_name) %>%
  select(-c(score, phase)) %>%
  mutate(seqid = paste0("chr", seqid))


zfin_exons_numbered <- bind_rows(
  zfin_granges_df %>%
    filter(type == "exon") %>%
    filter(strand == "+") %>%
    group_by(gene_name) %>%
    arrange(start, .by_group = T) %>%
    mutate(exon_number = paste0("exon_", row_number())),

  minus_exons <- zfin_granges_df %>%
    filter(type == "exon") %>%
    filter(strand == "-") %>%
    group_by(gene_name) %>%
    arrange(desc(start), .by_group = T) %>%
    mutate(exon_number = paste0("exon_", row_number()))
)

zfin_granges_df1 <- left_join(zfin_granges_df, zfin_exons_numbered) %>%
  mutate(label = ifelse(type == "exon", exon_number, type)) %>%
  select(-c(exon_number, source, ENSDAR))

# the granges object; made from zfin gene model
zfin_granges <- makeGRangesFromDataFrame(zfin_granges_df1, keep.extra.columns = T)

# the BS genome object for GRCz11
Drerio

