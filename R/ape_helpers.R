#' @importFrom purrr imap_chr map_chr
#' @importFrom stringr str_extract str_pad
dnastring_to_origin <-
   function(dna) {
      nucleotides <- as.character(dna)
      nucleotides <-
         str_split(gsub("(.{10})", "\\1 ", nucleotides),
                   pattern = " ",
                   n = Inf)
      nucleotides <- unlist(nucleotides)
      nucleotides <- nucleotides[nucleotides != ""]
      nucleotides <-
         purrr::imap_chr(.x = nucleotides, ~ paste0(.y, "_", .x))

      nucleotides <-
         purrr::map_chr(
            .x = nucleotides,
            .f = function(x) {
               digit <- as.numeric(stringr::str_extract(x, "^[:digit:]*"))
               new_digit <- (digit - 1) * 10 + 1
               chars <- stringr::str_extract(x, "[:alpha:]+")
               res <- paste0(new_digit, " ", chars)
               return(res)
            }
         )

      nucleotides <-
         purrr::map_chr(
            .x = nucleotides,
            .f = function(x) {
               digit <- as.numeric(stringr::str_extract(x, "^[:digit:]*"))
               if (digit %% 60 == 1) {
                  new_digit <- digit
                  new_digit <-
                     stringr::str_pad(as.character(new_digit),
                             side = "left",
                             width = 9)
               } else {
                  new_digit <- ""
               }
               if ((digit + 10) %% 60 == 1) {
                  end <- "\n"
               } else {
                  end <- ""
               }
               chars <- stringr::str_extract(x, "[:alpha:]+")
               res <- paste0(new_digit, " ", chars, end)
               return(res)
            }
         )


      origin <- c("ORIGIN\n", nucleotides)
      origin <- paste(origin, collapse = "")

      return(origin)
   }


#' @importFrom cli cli_abort
#' @importFrom GenomicRanges elementMetadata GRanges start end
#' @importFrom IRanges IRanges
parse_query <- function(query, genome, extend_left, extend_right, grange_use) {
  if (is.character(query)) {
    if (length(query) != 1) cli::cli_abort("You must provide only one gene name for the query.")
    if (query %notin% GenomicRanges::elementMetadata(grange_use)[, "gene_name"]) cli::cli_abort("Your gene name was not found.")
    query_grange <-
      grange_use[(GenomicRanges::elementMetadata(grange_use)[, "gene_name"] %in% query)]
    query_grange <-
      query_grange[(GenomicRanges::elementMetadata(query_grange)[, "type"] %in% "gene")]
    if (length(query_grange) != 1) cli::cli_abort("Your query returned multiple hits. This is not supported by this function.")
  } else if (is.numeric(query)) {
    if (length(query) != 3) cli::cli_abort("Your query must be a three-element numeric vector corresponding to chromosome, start and end.")
    if (query[2] >= query[3]) cli::cli_abort("Your query must be a three-element numeric vector corresponding to chromosome, start and end.")
    query_grange = GenomicRanges::GRanges(
      seqnames = paste0("chr", query["chr"]),
      ranges = IRanges::IRanges(start = query["start"],
                       end = query["end"])
    )
  } else {
    return("You must provide a unique genomic range or gene name for this function")
  }
  GenomicRanges::start(query_grange) <- GenomicRanges::start(query_grange) - extend_left
  GenomicRanges::end(query_grange) <- GenomicRanges::end(query_grange) + extend_right
  return(query_grange)

}

#' @importFrom BSgenome getSeq
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 Hsapiens
#' @importFrom BSgenome.Drerio.UCSC.danRer11 Drerio
#' @importFrom cli cli_abort
get_sequence <- function(genome, query_grange) {
  if (genome == "hg38") {
    dna <-
      suppressWarnings(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, query_grange))

  } else if (genome == "GRCz11") {
    dna <-
      suppressWarnings(BSgenome::getSeq(BSgenome.Drerio.UCSC.danRer11::Drerio, query_grange))


  } else {
    cli::cli_abort("The selected genome is not available.")
  }
  names(dna) <- "ape_seq"
  return(dna)
}

#' @importFrom lubridate day now month year
#' @importFrom Biostrings width
#' @importFrom stringr str_pad
make_locus_string <- function(dna, query) {
  date_string <-
    paste0(
      lubridate::day(lubridate::now()),
      "-",
      str_to_upper(lubridate::month(lubridate::now(), label = T)),
      "-",
      lubridate::year(lubridate::now())
    )
  length_string <- paste0(Biostrings::width(dna),
                          " bp ds-DNA")


  if (is.character(query)) {
    name_string <- query
  } else {
    name_string <- "custom range"
  }

  locus_string <-
    paste0(
      stringr::str_pad("LOCUS", width = 12, side = "right"),
      stringr::str_pad(name_string, width = 12, side = "right"),
      stringr::str_pad(length_string, width = 26, side = "left"),
      stringr::str_pad("linear", width = 11, side = "left"),
      stringr::str_pad(date_string, width = 19, side = "left")
    )
  return(locus_string)

}


#' @importFrom GenomicRanges seqnames start end
make_comment_string <-
  function(query,
           extend_left,
           extend_right,
           query_grange,
           genome) {
    if (is.numeric(query)) {
      query_text <-
        paste0("Locus is chr", query["chr"], " ", query["start"], "-", query["end"])
    }

    if (is.character(query)) {
      query_text <- paste0("Gene is ", query)
    }
    if (genome == "hg38") {
      species <- "Homo sapiens"
    } else if (genome == "GRCz11") {
      species  <- "Danio rerio"
    }
    comment_string <-
      paste0(
        "COMMENT     Sequence is from ", species , " ", genome, ".\nCOMMENT     Gene models are from ensembl.\nCOMMENT     ",
        query_text,
        "\nCOMMENT     Extensions ",
        extend_left,
        " bp left and ",
        extend_right,
        " bp right.",
        "\nCOMMENT     Final genomic coordinates are:\nCOMMENT     ",
        as.vector(GenomicRanges::seqnames(query_grange)),
        " ",
        GenomicRanges::start(query_grange),
        "-",
        GenomicRanges::end(query_grange)

      )
    return(comment_string)


  }

make_transcript_comment <- function(species, transcriptome, gene, query) {
  comment_string <-
    paste0(
      "COMMENT     cDNA sequence is from ",
      species ,
      " ",
      transcriptome,
      ".\nCOMMENT     Gene is ",
      gene,
      ".\nCOMMENT     Transcript is ",
      query

    )
  return(comment_string)

}

#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges elementMetadata start shift mcols
#' @importFrom GenomeInfoDb seqlevels seqlevelsInUse
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate select
make_grange <-
  function(grange_use,
           query_grange,
           include_type,
           additional_granges) {
    dna_grange <-
      IRanges::subsetByOverlaps(grange_use, query_grange)
    dna_grange <-
      dna_grange[(GenomicRanges::elementMetadata(dna_grange)[, "type"] %in% include_type)]
    #transform the seqnames
    GenomeInfoDb::seqlevels(dna_grange) <- GenomeInfoDb::seqlevelsInUse(dna_grange)
    GenomeInfoDb::seqlevels(dna_grange) <- "ape_seq"
    if (!is.null(additional_granges)) {
      GenomeInfoDb::seqlevels(additional_granges) <- "ape_seq"
    }
    # find the overall start
    overall_start <- GenomicRanges::start(query_grange)
    # shift the coordinates
    dna_grange <-
      GenomicRanges::shift(dna_grange, shift = -1 * (overall_start - 1))
    if (!is.null(additional_granges)) {
      additional_granges <-
        GenomicRanges::shift(additional_granges, shift = -1 * (overall_start - 1))
      # rename the metadata
      GenomicRanges::mcols(additional_granges) <-
        GenomicRanges::mcols(additional_granges) |>
        tibble::as_tibble() |>
        dplyr::mutate(locus_tag = paste0(gene_name, "_", label)) |>
        dplyr::select(-c(gene_name, label))

    }

    GenomicRanges::mcols(dna_grange) <-
      GenomicRanges::mcols(dna_grange) |>
      tibble::as_tibble() |>
      dplyr::mutate(locus_tag = paste0(gene_name, "_", label)) |>
      dplyr::select(-c(gene_name, label)) |>
      dplyr::mutate(
        fwdcolor = recode(
          type,
          "ncRNA_gene" = "#deebf7",
          "rRNA" = "#deebf7",
          "exon" = "#3182bd",
          "pseudogene" = "#deebf7",
          "pseudogenic_transcript" = "#deebf7",
          "ncRNA" = "#deebf7",
          "gene" = "#deebf7",
          "CDS" = "#deebf7",
          "lnc_RNA" = "#deebf7",
          "mRNA" = "#deebf7",
          "three_prime_UTR" = "#9ecae1",
          "five_prime_UTR" = "#9ecae1",
          "unconfirmed_transcript" = "#deebf7",
          "scRNA" = "#deebf7",
          "C_gene_segment" = "#deebf7",
          "D_gene_segment" = "#deebf7",
          "J_gene_segment" = "#deebf7",
          "V_gene_segment" = "#deebf7",
          "miRNA" = "#deebf7",
          "tRNA" = "#deebf7",
          "snRNA" = "#deebf7",
          "snoRNA" = "#deebf7",
          "lincRNA_gene" = "#deebf7",
          "unconfirmed_transcript" = "#deebf7",
          "lncRNA_gene" = "#deebf7"

        )
      ) |>
      dplyr::mutate(
        revcolor = recode(
          type,
          "ncRNA_gene" = "#e5f5e0",
          "rRNA" = "#e5f5e0",
          "exon" = "#31a354",
          "pseudogene" = "#e5f5e0",
          "pseudogenic_transcript" = "#e5f5e0",
          "ncRNA" = "#e5f5e0",
          "gene" = "#e5f5e0",
          "CDS" = "#e5f5e0",
          "lnc_RNA" = "#e5f5e0",
          "mRNA" = "#e5f5e0",
          "three_prime_UTR" = "#a1d99b",
          "five_prime_UTR" = "#a1d99b",
          "unconfirmed_transcript" = "#e5f5e0",
          "scRNA" = "#e5f5e0",
          "C_gene_segment" = "#e5f5e0",
          "D_gene_segment" = "#e5f5e0",
          "J_gene_segment" = "#e5f5e0",
          "V_gene_segment" = "#e5f5e0",
          "miRNA" = "#e5f5e0",
          "tRNA" = "#e5f5e0",
          "snRNA" = "#e5f5e0",
          "snoRNA" = "#e5f5e0",
          "lincRNA_gene" = "#e5f5e0",
          "unconfirmed_transcript" = "#e5f5e0",
          "lncRNA_gene" = "#e5f5e0",
        )
      )


    dna_grange <- c(dna_grange, additional_granges)
    names(dna_grange) <- dna_grange$locus_tag
    return(dna_grange)


  }


#' @importFrom AnnotationDbi select
get_transcript_genename <- function(query, org, txdb) {
  gene_id <- AnnotationDbi::select(txdb,
                        keytype = "TXNAME",
                        keys = query,
                        columns = "GENEID")
  gene_name <- AnnotationDbi::select(org,
                                      keytype = "ENTREZID",
                                      columns = "SYMBOL",
                                      keys = gene_id$GENEID)
  return(gene_name$SYMBOL)

}



