#' Construct an APE Object from the hg38 Reference
#'
#' @description This function takes either a gene name or sequence coordinates and returns an ape object with DNA sequence from the GRCh38 reference. Features overlapping the query range (gene or sequence) are returned as a GRanges object and as features within the Ape object.  The features included can optionally be filtered using the "include_type" argument to the function.  User-defined features can be added at the point of extracting the sequence and creating the object.  This is useful if features are predefined with GRCh38 sequence coordinates since the coordinates are renamed relative to the query when the Ape object is made.  Using the "additional_granges" argument, you can provide additional features not present in the Zfin gene model which will be added to the Ape object GRanges slot and to the features slot.
#'
#' @param query Either a valid gene name or a named numeric vector of genome coordinates.  This vector should be of the form:  c(chr = 1, start = 1000, end = 2000).  The vector must be numeric an must have those names.  The chromosome number will be converted to "chr1" etc internally.
#' @param extend_left Number of bases to extend the query to the left or "upstream" relative to the + strand.
#' @param extend_right Number of bases to extend the query to the right or "downstream" relative to the + strand.
#' @param include_type The type of features to include from the zfin gene model.
#' @param additional_granges A GRanges object with features to add to the Ape Object.  Coordinates should all be relative to the GRCz11 reference, *NOT* the sequence extracted for the ape file.  The Granges object can be constructed with the following syntax:  GenomicRanges::makeGRangesFromDataFrame(data.frame(seqname = "chr6", start = 40523370, end = 40523380, strand = "+", type = "addl_feature", gene_name = "some_gene", label = "feature1"), keep.extra.columns = T).  The gene_name argument here is optional.  If you have defined features based on the extracted sequence, (i.e. relative to position 1 in the ORIGIN section of the Ape object), the best option is to use the feature setting function FEATURES(instance_of_Ape) <- GRanges_Object.
#' @return An Ape object
#' @export
#' @import Biostrings GenomicRanges tidyverse IRanges GenomeInfoDb
bb_hg38_ape <-
   function(query,
            extend_left = 0,
            extend_right = 0,
            include_type = c(
               "ncRNA_gene",
               "rRNA",
               "exon",
               "pseudogene",
               "pseudogenic_transcript",
               "ncRNA",
               "gene",
               "CDS",
               "lnc_RNA",
               "mRNA",
               "three_prime_UTR",
               "five_prime_UTR",
               "unconfirmed_transcript",
               "scRNA",
               "C_gene_segment",
               "D_gene_segment",
               "J_gene_segment",
               "V_gene_segment",
               "miRNA",
               "tRNA",
               "snRNA",
               "snoRNA"
            ),
            additional_granges = NULL) {
      # first extract the biostring from the genome
      if (is.character(query)) {
         stopifnot("You must provide only one gene name for the query." = length(query) == 1)
         stopifnot("Your gene name was not found." = any(elementMetadata(hg38_granges)[, "gene_name"] %in% query))
         query_grange <-
            hg38_granges[(elementMetadata(hg38_granges)[, "gene_name"] %in% query)]
         query_grange <-
            query_grange[(elementMetadata(query_grange)[, "type"] %in% "gene")]
         stopifnot(
            "Your query returned multiple hits. This is not supported by this function." = length(query_grange) == 1
         )
      } else if (is.numeric(query)) {
         stopifnot(
            "Your query must be a three-element numeric vector corresponding to chromosome, start and end." = length(query) == 3
         )
         stopifnot(
            "Your query must be a three-element numeric vector corresponding to chromosome, start and end." = query[2] < query[3]
         )
         query_grange = GenomicRanges::GRanges(
            seqnames = paste0("chr", query["chr"]),
            ranges = IRanges(start = query["start"],
                             end = query["end"])
         )
      } else {
         return("You must provide a unique genomic range or gene name for this function")
      }
      start(query_grange) <- start(query_grange) - extend_left
      end(query_grange) <- end(query_grange) + extend_right
      dna <-
         suppressWarnings(BSgenome::getSeq(Hsapiens, query_grange))
      names(dna) <- "ape_seq"
      dna_grange <-
         IRanges::subsetByOverlaps(hg38_granges, query_grange)
      dna_grange <-
         dna_grange[(elementMetadata(dna_grange)[, "type"] %in% include_type)]
      #transform the seqnames
      seqlevels(dna_grange) <- seqlevelsInUse(dna_grange)
      seqlevels(dna_grange) <- "ape_seq"
      if (!is.null(additional_granges)) {
         seqlevels(additional_granges) <- "ape_seq"
      }
      # find the overall start
      overall_start <- start(query_grange)
      # shift the coordinates
      dna_grange <-
         GenomicRanges::shift(dna_grange, shift = -1 * (overall_start - 1))
      if (!is.null(additional_granges)) {
         additional_granges <-
            GenomicRanges::shift(additional_granges, shift = -1 * (overall_start - 1))
      # rename the metadata
         additional_granges <- additional_granges %>%
            plyranges::mutate(locus_tag = paste0(gene_name, "_", label)) %>%
            plyranges::select(-c(gene_name, label))

      }


      dna_grange <- dna_grange %>%
         plyranges::mutate(locus_tag = paste0(gene_name, "_", label)) %>%
         plyranges::select(-c(gene_name, label)) %>%
         plyranges::mutate(
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
               "snoRNA" = "#deebf7"
            )
         ) %>%
         plyranges::mutate(
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
               "snoRNA" = "#e5f5e0"
            )
         )

      dna_grange <- c(dna_grange, additional_granges)
      names(dna_grange) <- dna_grange$locus_tag
      # construct the locus text
      date_string <-
         paste0(
            lubridate::day(lubridate::now()),
            "-",
            str_to_upper(lubridate::month(lubridate::now(), label = T)),
            "-",
            lubridate::year(lubridate::now())
         )

      length_string <- paste0(width(dna),
                              " bp ds-DNA")


      if (is.character(query)) {
         name_string <- query
      } else {
         name_string <- "custom range"
      }

      locus_string <-
         paste0(
            str_pad("LOCUS", width = 12, side = "right"),
            str_pad(name_string, width = 12, side = "right"),
            str_pad(length_string, width = 26, side = "left"),
            str_pad("linear", width = 11, side = "left"),
            str_pad(date_string, width = 19, side = "left")
         )

      if (is.numeric(query)) {
         query_text <-
            paste0("Locus is chr", query["chr"], " ", query["start"], "-", query["end"])
      }

      if (is.character(query)) {
         query_text <- paste0("Gene is ", query)
      }

      comment_string <-
         paste0(
            "COMMENT     Sequence is from Homo sapiens, GRCh38.\nCOMMENT     Gene models are from ensembl.\nCOMMENT     ",
            query_text,
            "\nCOMMENT     Extensions ",
            extend_left,
            " bp left and ",
            extend_right,
            " bp right.",
            "\nCOMMENT     Final genomic coordinates are:\nCOMMENT     ",
            as.vector(seqnames(query_grange)),
            " ",
            start(query_grange),
            "-",
            end(query_grange)

         )

      ape_instance <- Ape(
         LOCUS = locus_string,
         COMMENT = comment_string,
         FEATURES = granges_to_features(dna_grange),
         ORIGIN = dnastring_to_origin(dna),
         dna_biostring = dna,
         granges = dna_grange
      )

      return(ape_instance)
   }
