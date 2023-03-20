#' @title Make an Ape Genome Object
#' @description This function takes either a gene name or sequence coordinates and returns an ape object with DNA sequence from the the selected genome reference.  You can choose from hg38 and GRCz11.  Features are generated from ensembl GFF files.  Features overlapping the query range (gene or sequence) are returned as a GRanges object and as features within the Ape object.  The features included can optionally be filtered using the "include_type" argument to the function. Using the "additional_granges" argument, you can provide additional features not present in the standard gene model which will be added to the Ape object GRanges slot and to the features slot.
#'
#' @param query Either a valid gene name or a named numeric vector of genome coordinates.  This vector should be of the form:  c(chr = 1, start = 1000, end = 2000).  The vector must be numeric an must have those names.  The chromosome number will be converted to "chr1" etc internally.
#' @param genome The genome to pull from, Default: c("hg38", "GRCz11")
#' @param extend_left Number of bases to extend the query to the left or "upstream" relative to the + strand.
#' @param extend_right Number of bases to extend the query to the right or "downstream" relative to the + strand.
#' @param include_type The type of features to include from the standard gene model. Default: c("ncRNA_gene", "rRNA", "exon", "pseudogene", "pseudogenic_transcript",
#'    "ncRNA", "gene", "CDS", "lnc_RNA", "mRNA", "three_prime_UTR",
#'    "five_prime_UTR", "unconfirmed_transcript", "scRNA", "C_gene_segment",
#'    "D_gene_segment", "J_gene_segment", "V_gene_segment", "miRNA",
#'    "tRNA", "snRNA", "snoRNA", "lincRNA_gene", "lncRNA_gene",
#'    "unconfirmed_transcript")
#' @param additional_granges A GRanges object with features to add to the Ape Object.  Coordinates should all be relative to the reference, *NOT* the sequence extracted for the ape file.  The Granges object can be constructed with the following syntax:  GenomicRanges::makeGRangesFromDataFrame(data.frame(seqname = "chr6", start = 40523370, end = 40523380, strand = "+", type = "addl_feature", gene_name = "prkcda", label = "feature1"), keep.extra.columns = T).  The gene_name argument here is optional.  If you have defined features based on the extracted sequence, (i.e. relative to position 1 in the ORIGIN section of the Ape object), the best option is to use the feature setting function FEATURES(instance_of_Ape) <- GRanges_Object.
#' @return An APE object
#' @seealso
#'  \code{\link[blaseRdata]{hg38_granges}}, \code{\link[blaseRdata]{zfin_granges}}
#'  \code{\link[cli]{cli_abort}}
#' @rdname bb_make_ape_genomic
#' @export
#' @importFrom blaseRdata hg38_granges zfin_granges
#' @importFrom cli cli_abort
bb_make_ape_genomic <-
  function(query,
           genome = c("hg38", "GRCz11"),
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
             "snoRNA",
             "lincRNA_gene",
             "lncRNA_gene",
             "unconfirmed_transcript"
           ),
           additional_granges = NULL) {
    genome <- match.arg(genome)
    include_type <- match.arg(include_type, several.ok = TRUE)
    if (genome == "hg38") {
      granges_use <- blaseRdata::hg38_granges
    } else if (genome == "GRCz11") {
      granges_use <- blaseRdata::zfin_granges
    } else {
      cli::cli_abort("You must choose either hg38 or GRCz11 genomes.")
    }
    # return(granges_use)

    query_grange <- parse_query(query,
                                genome,
                                extend_left,
                                extend_right,
                                granges_use)
    dna <-
      get_sequence(genome, query_grange)

    locus_string <- make_locus_string(dna, query)

    dna_grange <- make_grange(grange_use,
                              query_grange,
                              include_type,
                              additional_granges)

    ape_features <- granges_to_features(dna_grange)

    ape_origin <- dnastring_to_origin(dna)

    comment_string <-
      make_comment_string(query, extend_left, extend_right, query_grange, genome)

    ape_instance <- Ape(
      LOCUS = locus_string,
      COMMENT = comment_string,
      FEATURES = ape_features,
      ORIGIN = ape_origin,
      dna_biostring = dna,
      granges = dna_grange
    )

    return(ape_instance)

  }

#' @title Make an Ape Transcriptome Object
#' @description This function takes a specific ensembl transcript identifier, such as ENST00000348343.11, and gets the cDNA sequence from the corresponding transcriptomic reference.  This is returned as an Ape object with the UTR's and the CDS annotated as features.
#' @param query A specific ensembl transcript identifier.
#' @param transcriptome Genome/transcriptome reference to use, Default: c("hg38", "GRCz11")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[TxDb.Hsapiens.UCSC.hg38.knownGene]{TxDb.Hsapiens.UCSC.hg38.knownGene}}
#'  \code{\link[BSgenome.Hsapiens.UCSC.hg38]{BSgenome.Hsapiens.UCSC.hg38}}
#'  \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}
#'  \code{\link[blaseRdata]{character(0)}}
#'  \code{\link[BSgenome.Drerio.UCSC.danRer11]{BSgenome.Drerio.UCSC.danRer11}}
#'  \code{\link[org.Dr.eg.db]{org.Dr.eg.db}}
#'  \code{\link[cli]{cli_abort}}
#'  \code{\link[Biostrings]{matchPattern}}
#' @rdname bb_make_ape_transcript
#' @export
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom blaseRdata TxDb.Drerio.UCSC.danRer11.ensGene
#' @importFrom BSgenome.Drerio.UCSC.danRer11 BSgenome.Drerio.UCSC.danRer11
#' @importFrom org.Dr.eg.db org.Dr.eg.db
#' @importFrom cli cli_abort
#' @importFrom Biostrings matchPattern
#' @importFrom GenomicFeatures exonsBy cdsBy extractTranscriptSeqs
bb_make_ape_transcript <- function(query,
                                   transcriptome = c("hg38", "GRCz11")
                                   ) {
  transcriptome <- match.arg(transcriptome)
  if (transcriptome == "hg38") {
    txdb <-
      TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    genome <-
      BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    org <- org.Hs.eg.db::org.Hs.eg.db
    species <- "Homo sapiens"
  } else if (transcriptome == "GRCz11") {
    txdb <- blaseRdata::TxDb.Drerio.UCSC.danRer11.ensGene
    genome <-
      BSgenome.Drerio.UCSC.danRer11::BSgenome.Drerio.UCSC.danRer11
    org <- org.Dr.eg.db::org.Dr.eg.db
    species <- "Danio rerio"
  } else {
    cli::cli_abort("You must choose either hg38 or GRCz11 genomes.")
  }
  gene_name <- get_transcript_genename(query, org, txdb)
  transcripts <-
    suppressWarnings(GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE))
  transcript <- transcripts[query]
  tx_seq <- GenomicFeatures::extractTranscriptSeqs(genome, transcript)
  cdss <- suppressWarnings(GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE))
  cds <- cdss[query]
  cds_seq <- GenomicFeatures::extractTranscriptSeqs(genome, cds)
  overlap <-
    Biostrings::matchPattern(pattern = cds_seq[[1]], subject = tx_seq[[1]])
  cds_grange <- GRanges(
    seqnames = "ape_transcript",
    ranges = IRanges(start = start(overlap),
                     width = width(overlap)),
    locus_tag = paste0(gene_name,"_cds"),
    type = "cds",
    fwdcolor = "#deebf7",
    revcolor = "#deebf7"
  )
  utr_5_grange <- GRanges(
    seqnames = "ape_transcript",
    ranges = IRanges(start = 1,
                     end = start(overlap) - 1),
    locus_tag = paste0(gene_name, "_five_prime_UTR"),
    type = "five_prime_UTR",
    fwdcolor = "#9ecae1",
    revcolor = "#9ecae1"
  )
  utr_3_grange <- GRanges(
    seqnames = "ape_transcript",
    ranges = IRanges(
      start = start(overlap) +
        width(overlap),
      end = width(tx_seq)
    ),
    locus_tag = paste0(gene_name, "_three_prime_UTR"),
    type = "three_prime_UTR",
    fwdcolor = "#9ecae1",
    revcolor = "#9ecae1"
  )
  tx_grange <- c(utr_5_grange, cds_grange, utr_3_grange)
  seqs <- c(tx_seq, cds_seq)
  names(seqs) <- c("cDNA", "cds")

  ape_features <- granges_to_features(tx_grange)

  ape_origin <- dnastring_to_origin(tx_seq)


  comment_string <- make_transcript_comment(species,
                                            transcriptome,
                                            gene_name,
                                            query)

  locus_string <- make_locus_string(dna = seqs["cDNA"],
                                    query = query)

  ape_instance <- Ape(
    LOCUS = locus_string,
    COMMENT = comment_string,
    FEATURES = ape_features,
    ORIGIN = ape_origin,
    dna_biostring = seqs,
    granges = tx_grange#
  )

  return(ape_instance)


}
