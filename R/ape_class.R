#' An S4 class to hold genebank/APE file data.
#'
#' @description An instance of this class is best created by calling "bb_parseape()" on a genebank or APE-formatted file.  That function will parse the file, correctly format the sections and place them in the slots of the Ape Object.  Technically only "LOCUS" is a required slot for the Ape object, however there is no point without having "ORIGIN" (sequence data), and so bb_parseape() will fail without an "ORIGIN" section.  Other slots are optional.  Additional slots will be ignored by the constructor function.  DNA sequence will be stored in a DNAStringSet object and features in a GRanges object.  See https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html for genebank file specification.
#'
#' @slot LOCUS The LOCUS line of the genebank formatted as a character string.
#' @slot DEFINITION The DEFINITION line of the genebank file formatted as a character string.
#' @slot ACCESSION The ACCESSION section of the genebank file formatted as a character string.
#' @slot VERSION The VERSION section of the genebank file formatted as a character string.
#' @slot SOURCE The SOURCE section of the genebank file formatted as a character string.
#' @slot COMMENT The COMMENT section of the genebank file formatted as a character string.
#' @slot FEATURES The FEATURES section of the genebank file formatted as a character string.  Created internally from the GRanges object.  Caution:  some FEATURE attributes may be lost in conversion.
#' @slot ORIGIN The DNA sequence
#' @slot end_of_file The end of the file signal.
#' @slot dna_biostring The entire ORIGIN sequence formatted as a DNAStringSet of length 1.
#' @slot granges Genebank features formatted as a GRanges object.
#' @import methods
#' @export Ape
#' @exportClass Ape
Ape <- setClass(
  "Ape",
  slots = list(
    LOCUS = "character",
    DEFINITION = "character",
    ACCESSION = "character",
    VERSION = "character",
    SOURCE = "character",
    COMMENT = "character",
    FEATURES = "character",
    ORIGIN = "character",
    end_of_file = "character",
    dna_biostring = "DNAStringSet",
    granges = "GRanges"
  ),
  prototype = list(
    LOCUS = NA_character_,
    DEFINITION = NA_character_,
    ACCESSION = NA_character_,
    VERSION = NA_character_,
    SOURCE = NA_character_,
    COMMENT = NA_character_,
    FEATURES = NA_character_,
    ORIGIN = NA_character_,
    end_of_file = "//"
  )
)


#' Show an Ape Object
#'
#' @export
setMethod("show",
          "Ape",
          function(object) {
            if(!all(is.na(object@LOCUS))) cat(object@LOCUS, "\n")
            if(!all(is.na(object@DEFINITION))) cat(object@DEFINITION, sep = "\n")
            if(!all(is.na(object@ACCESSION))) cat(object@ACCESSION, sep = "\n")
            if(!all(is.na(object@VERSION))) cat(object@VERSION, sep = "\n")
            if(!all(is.na(object@SOURCE))) cat(object@SOURCE, sep = "\n")
            if(!all(is.na(object@COMMENT))) cat(object@COMMENT, sep = "\n")
            if(!all(is.na(object@FEATURES))) cat(object@FEATURES, sep = "\n")
            if(!all(is.na(object@ORIGIN))) cat(object@ORIGIN, sep = "\n")
            cat(object@end_of_file)
          }
)

#' Get the GRanges Slot from an Ape Object
#'
#' @export
setGeneric("Ape.granges", function(ape) standardGeneric("Ape.granges"))
#' @export
setMethod("Ape.granges", "Ape", function(ape) ape@granges)

#' Get the DNASringSet Slot from an Ape Object
#'
#' @export
setGeneric("Ape.DNA", function(ape) standardGeneric("Ape.DNA"))
#' @export
setMethod("Ape.DNA", "Ape", function(ape) ape@dna_biostring)

#' Get the Features Slot from an Ape Object
#'
#' @description This function cats the features slot from an Ape object.
#'
#' @export
setGeneric("FEATURES", function(ape) standardGeneric("FEATURES"))
#' @export
setMethod("FEATURES", "Ape", function(ape) cat(ape@FEATURES, sep = "\n"))

#' Get the Locus Slot from an Ape Object
#'
#' @description This function cats the Locus slot from an Ape object.
#'
#' @export
setGeneric("LOCUS", function(ape) standardGeneric("LOCUS"))
#' @export
setMethod("LOCUS", "Ape", function(ape) cat(ape@LOCUS, sep = "\n"))

#' Get the Comments Slot from an Ape Object
#'
#' @description This function gets the comments slot from an Ape object.
#'
#' @export
setGeneric("COMMENTS", function(ape) standardGeneric("COMMENTS"))
#' @export
setMethod("COMMENTS", "Ape", function(ape) cat(ape@COMMENTS, sep = "\n"))

#' Set the FEATURES Slot of a GRanges Object
#'
#' @param ape An ape object
#' @param gr A GRanges object.  This object will become the new FEATURES and granges slots for the Ape object.  So if you want to keep the old features, the new features need to be appended using c(old_gr, new_gr) as the value for the gr argument.
#' @export
setGeneric("Ape.setFeatures", function(ape, gr)
  standardGeneric("Ape.setFeatures"))
#' @export
setMethod("Ape.setFeatures", "Ape", function(ape, gr) {
  ape@granges <- gr
  ape@FEATURES <- blaseRtools::granges_to_features(gr)
  validObject(ape)
  ape
})


# other methods
#' Save an Ape Instance as a Genebank Format File
#'
#' @param out Name of genebank/APE file to write
#' @export
setGeneric("Ape.save", function(ape, out) standardGeneric("Ape.save"))
#' @export
setMethod("Ape.save", "Ape", function(ape, out) capture.output(ape, file = out))

#' Save an Ape Instance as a Fasta File
#'
#' @param out Name of FASTA file to write
#' @param feature Name of feature to select when writing FASTA file.  If null (default), the whole biostring will be saved as a fasta.
#' @export
setGeneric("Ape.fasta", function(ape, feature = NULL, out) standardGeneric("Ape.fasta"))
#' @export
setMethod("Ape.fasta", "Ape", function(ape, feature = NULL, out) {
  if (is.null(feature)) {
   Biostrings::writeXStringSet(Ape.DNA(ape), file = out, format = "fasta")
  } else {
   stopifnot("Requested feature is not in GRanges" = feature %in% names(ape@granges))
   seq <- getSeq(Ape.DNA(ape), ape@granges[feature])
   Biostrings::writeXStringSet(seq, file = out, format = "fasta")
  }

  })

#' Run FIMO on Selected Ape Object Features
#'
#' @description For the supplied Ape object, run FIMO to identify putative transcription factor binding sites in a DNA subsequence.
#'
#' @param ape An Ape instance
#' @param fimo_feature A character vector of features from the Ape object that will be used to run fimo.
#' @param out Directory that will be created to hold the fimo results.  A date/time stamp will be appended.  If null, the objects will not be saved and the function will only return a GRanges object
#' @import tidyverse GenomicRanges circlize
#' @export
setGeneric("Ape.fimo", function(ape, fimo_feature, out = NULL)
  standardGeneric("Ape.fimo"))
#' @import tidyverse GenomicRanges circlize
#' @export
setMethod("Ape.fimo", "Ape", function(ape, fimo_feature, out = NULL) {
  # check if fimo exists
  fimocheck <- Sys.which("fimo")

  stopifnot(
    "Go to:  https://meme-suite.org/meme/meme_5.3.2/doc/install.html?man_type=web.\nInstall meme suite tools using the Quick Install Instructions.\nThen try again." = str_detect(fimocheck, "meme/bin/fimo")
  )

  # get the fasta from the ape
  if (is.null(out)) {
    outdir <- "fimo_temp"
  } else {
    ts <-
      str_replace_all(Sys.time(), "[:punct:]|[:space:]|[:alpha:]", "")
    outdir <- paste0(out, "_", ts)
  }

  dir.create(outdir, recursive = T)

  # run fimo for all given features
  readr::write_lines(meme, file = paste0(outdir, "/meme.motif"))

  all_fimo_res <- map_dfr(
    .x = fimo_feature,
    .f = function(x,
                  the_ape = ape,
                  od = outdir,
                  fc = fimocheck) {

      dir.create(file.path(od, x), recursive = T, showWarnings = F)

      Ape.fasta(
        the_ape,
        feature = x,
        out = paste0(od, "/", x, "/", x, ".fa")
      )

      cmd <- paste0(
        fc,
        " -o ",
        od,
        "/",
        x,
        "/fimo ",
        od,
        "/meme.motif ",
        od,
        "/",
        x,
        "/",
        x,
        ".fa"
      )
      message(cmd, "\n")
      system(cmd)

      # ingest the tsv
      fimo_res <-
        read_tsv(
          str_glue("{od}/{x}/fimo/fimo.tsv"),
          skip_empty_rows = T,
          comment = "#",
          show_col_types = F
        ) %>%
        filter(!is.na(motif_alt_id))

      # calculate the startpoint relative to the original Ape object
      ape_gr <- Ape.granges(the_ape)
      fimo_feature_start <- start(ape_gr[x])
      # calculate the number of bases to add
      bases_to_add <- fimo_feature_start - 1
      # modify the start and stop in the dataframe
      fimo_res <- fimo_res %>%
        mutate(start = start + bases_to_add) %>%
        mutate(stop = stop + bases_to_add) # strand should be unchanged

      return(fimo_res)
    }
  ) %>%
    mutate(feature = sequence_name) %>%
    mutate(sequence_name = "ape_seq") %>%
    mutate(type = "fimo_feature") %>%
    mutate(strand_1 = recode(strand, "+" = "plus", "-" = "minus")) %>%
    mutate(locus_tag = paste(motif_alt_id, start, stop, strand_1, sep = "_"))

  # make the color function
  fwd_col_fun <- circlize::colorRamp2(breaks = c(min(all_fimo_res$`q-value`), max(all_fimo_res$`q-value`)),
                                      colors = c("red", "white"))

  rev_col_fun <- circlize::colorRamp2(breaks = c(min(all_fimo_res$`q-value`), max(all_fimo_res$`q-value`)),
                                      colors = c("purple", "white"))

  all_fimo_res <- all_fimo_res %>%
    mutate(fwdcolor = str_sub(fwd_col_fun(all_fimo_res$`q-value`), start = 1, end = 7)) %>%
    mutate(revcolor = str_sub(rev_col_fun(all_fimo_res$`q-value`), start = 1, end = 7)) %>%
    select(-c(motif_alt_id, score, `p-value`, matched_sequence, strand_1))

  # remove meme.motif
  unlink(file.path(outdir, "meme.motif"))

  # remove the temp dir
  if (is.null(out)) {
    unlink(outdir, recursive = T)
  }

  all_fimo_res <- GenomicRanges::makeGRangesFromDataFrame(df = all_fimo_res, keep.extra.columns = T, seqnames.field = "sequence_name")
  names(all_fimo_res) <- all_fimo_res$locus_tag
  return(all_fimo_res)
})


# Validity Check
#' @export
setValidity("Ape", function(object) {
  if (is.na(object@LOCUS))
    return("LOCUS slot must be defined")
  # object@FEATURES <- granges_to_features(object@granges)
  if (any(object@FEATURES != granges_to_features(object@granges)))
    return("The granges slot must resolve to FEATURES slot using granges_to_features")
  TRUE

})
