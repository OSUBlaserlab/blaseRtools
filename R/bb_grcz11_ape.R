#' Format a DNAString as a Character Vector For Inclusion in Genebank File
#'
#' @description This function takes a DNAString or DNAStringSet object and converts to a character string.  Spaces, newlines and line numbers are introduced to reproduce a 60-nucleotide/80-character-wide genebank ORIGIN section.
#'
#' @param dna A DNAString or a DNAStringSet of length 1.
#' @return A character vector formatted for an Ape object ORIGIN slot.
#' @import Biostrings tidyverse
#' @export
dnastring_to_origin <-
   function(dna) {
      nucleotides <- as.character(dna)
      nucleotides <- str_split(gsub("(.{10})", "\\1 ", nucleotides), pattern = " ", n = Inf)
      nucleotides <- unlist(nucleotides)
      nucleotides <- imap_chr(.x = nucleotides, ~paste0(.y,"_", .x))

      nucleotides <- map_chr(.x = nucleotides, .f = function(x) {
         digit <- as.numeric(str_extract(x, "^[:digit:]*"))
         new_digit <- (digit-1)*10+1
         chars <- str_extract(x, "[:alpha:]+")
         res <- paste0(new_digit, " ", chars)
         return(res)
      })

      nucleotides <- map_chr(.x = nucleotides, .f = function(x) {
         digit <- as.numeric(str_extract(x, "^[:digit:]*"))
         if(digit%%60 == 1) {
            new_digit <- digit
            new_digit <- str_pad(as.character(new_digit), side = "left", width = 9)
         } else {
            new_digit <- ""
         }
         if((digit+10)%%60 == 1) {
            end <- "\n"
         } else {
            end <- ""
         }
         chars <- str_extract(x, "[:alpha:]+")
         res <- paste0(new_digit, " ", chars, end)
         return(res)
      })


      origin <- c("ORIGIN\n",nucleotides)
      origin <- paste(origin, collapse = "")

      return(origin)
   }

#' Construct an APE Object from the GRCz11 Reference
#'
#' @description This function takes either a gene name or sequence coordinates and returns an ape object with DNA sequence from the GRCz11 database.  The current (as of April 2018) zfin gene model (https://zfin.org/downloads/zfin_genes.gff3) was used to define gene structural features.  Features overlapping the query range (gene or sequence) are returned as a GRanges object and as features within the Ape object.  The features included can optionally be filtered using the "include_type" argument to the function.  User-defined features can be added at the point of extracting the sequence and creating the object.  This is useful if features are predefined with GRCz11 sequence coordinates since the coordinates are renamed relative to the query when the Ape object is made.  Using the "additional_granges" argument, you can provide additional features not present in the Zfin gene model which will be added to the Ape object GRanges slot and to the features slot.
#'
#' @param query Either a valid gene name or a named numeric vector of genome coordinates.  This vector should be of the form:  c(chr = 1, start = 1000, end = 2000).  The vector must be numeric an must have those names.  The chromosome number will be converted to "chr1" etc internally.
#' @param extend_left Number of bases to extend the query to the left or "upstream" relative to the + strand.
#' @param extend_right Number of bases to extend the query to the right or "downstream" relative to the + strand.
#' @param include_type The type of features to include from the zfin gene model.
#' @param additional_granges A GRanges object with features to add to the Ape Object.  Coordinates should all be relative to the GRCz11 reference, *NOT* the sequence extracted for the ape file.  The Granges object can be constructed with the following syntax:  GenomicRanges::makeGRangesFromDataFrame(data.frame(seqname = "chr6", start = 40523370, end = 40523380, strand = "+", type = "addl_feature", gene_name = "prkcda", label = "feature1"), keep.extra.columns = T).  The gene_name argument here is optional.  If you have defined features based on the extracted sequence, (i.e. relative to position 1 in the ORIGIN section of the Ape object), the best option is to use the feature setting function FEATURES(instance_of_Ape) <- GRanges_Object.
#' @return An Ape object
#' @export
#' @import Biostrings GenomicRanges tidyverse plyranges
bb_grcz11_ape <-
   function(query,
            extend_left = 0,
            extend_right = 0,
            include_type = c(
               "gene",
               "exon",
               "five_prime_UTR",
               "gene",
               "J_gene_segment",
               "lincRNA_gene",
               "lnc_RNA",
               "lncRNA_gene",
               "mRNA",
               "pseudogene",
               "pseudogenic_transcrpt",
               "three_prime_UTR",
               "unconfirmed_transcript"
            ),
            additional_granges = NULL) {
      # first extract the biostring from the genome
      if (is.character(query)) {
         stopifnot("You must provide only one gene name for the query." = length(query) == 1)
         stopifnot("Your gene name was not found." = any(elementMetadata(zfin_granges)[, "gene_name"] %in% query))
         query_grange <-
            zfin_granges[(elementMetadata(zfin_granges)[, "gene_name"] %in% query)]
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
      dna <- BSgenome::getSeq(Drerio, query_grange)
      names(dna) <- "ape_seq"
      dna_grange <- subsetByOverlaps(zfin_granges, query_grange)
      dna_grange <-
         dna_grange[(elementMetadata(dna_grange)[, "type"] %in% include_type)]
      dna_grange <- c(dna_grange, additional_granges)
      #transform the seqnames
      seqlevels(dna_grange) <- seqlevelsInUse(dna_grange)
      seqlevels(dna_grange) <- "ape_seq"
      # find the overall start
      overall_start <- start(query_grange)
      # shift the coordinates
      dna_grange <-
         GenomicRanges::shift(dna_grange, shift = -1 * (overall_start - 1))
      # rename the metadata
      dna_grange <- dna_grange %>%
         mutate(locus_tag = paste0(gene_name, "_", label)) %>%
         select(-c(gene_name, label)) %>%
         mutate(fwdcolor = "red", revcolor = "green")
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

      if(is.numeric(query)) {
         query_text <- paste0("Locus is chr", query["chr"], " ", query["start"], "-", query["end"])
      }

      if(is.character(query)) {
         query_text <- paste0("Gene is ", query)
      }

      comment_string <-
         paste0(
            "COMMENT     Sequence is from Danio rerio, GRCz11.\nCOMMENT     Gene models are from Zfin.\nCOMMENT     ",
            query_text,
            "\nCOMMENT     Extensions ",
            extend_left,
            " bp left and ",
            extend_right,
            " bp right.",
            "\nCOMMENT     Final genomic coordinates are:\nCOMMENT     ",
            as.vector(seqnames(query_grange)), " ", start(query_grange), "-", end(query_grange)

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
