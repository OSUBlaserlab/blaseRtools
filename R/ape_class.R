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

# Getter methods
#' Get the GRanges Slot from an Ape Object
#'
#' @export
setGeneric("GRanges", function(x) standardGeneric("GRanges"))
#' @export
setMethod("GRanges", "Ape", function(x) x@granges)

#' Get the DNASringSet Slot from an Ape Object
#'
#' @export
setGeneric("DNAStringSet", function(x) standardGeneric("DNAStringSet"))
#' @export
setMethod("DNAStringSet", "Ape", function(x) x@dna_biostring)

#' Get the Features Slot from an Ape Object
#'
#' @description This function gets the features slot from an Ape object.  Depending on console settings this may not be easily readable.  Better to simply show the whole Ape Object.  This may be useful for debugging.
#'
#' @export
setGeneric("FEATURES", function(x) standardGeneric("FEATURES"))
#' @export
setMethod("FEATURES", "Ape", function(x) x@FEATURES)

# Setter methods
#' Set the FEATURES Slot of a GRanges Object
#'
#' @param x An ape object
#' @param gr A GRanges object.  This object will become the new FEATURES and granges slots for the Ape object.  So if you want to keep the old features, the new features need to be appended using c(old_gr, new_gr) as the value for the gr argument.
#' @export
setGeneric("FEATURES<-", function(x, gr) standardGeneric("FEATURES<-"))
#' @export
setReplaceMethod("FEATURES", "Ape", function(x, gr) {
  # value should be a GRanges object formatted exactly how you want the features
  x@granges <- gr
  # write out the string in the necessary format
  x@FEATURES <- granges_to_features(gr)
  validObject(x)
  x
})


# other methods
#' Save an Ape Instance as a Genebank Format File
#'
#' @param out Name of genebank/APE file to write
#' @export
setGeneric("Ape.save", function(x, out) standardGeneric("Ape.save"))
#' @export
setMethod("Ape.save", "Ape", function(x, out) capture.output(x, file = out))

#' Save an Ape Instance as a Fasta File
#'
#' @param out Name of FASTA file to write
#' @param feature Name of feature to select when writing FASTA file.  If null (default), the whole biostring will be saved as a fasta.
#' @export
setGeneric("Ape.fasta", function(x, feature = NULL, out) standardGeneric("Ape.fasta"))
#' @export
setMethod("Ape.fasta", "Ape", function(x, feature = NULL, out) {
  if (is.null(feature)) {
   Biostrings::writeXStringSet(DNAStringSet(x), file = out, format = "fasta")
  } else {
   stopifnot("Requested feature is not in GRanges" = feature %in% names(x@granges))
   seq <- getSeq(DNAStringSet(x), x@granges[feature])
   Biostrings::writeXStringSet(seq, file = out, format = "fasta")
  }

  })

# Validity Check
#' @export
setValidity("Ape", function(object) {
  if (is.na(object@LOCUS))
    return("LOCUS slot must be defined")
  if (any(object@FEATURES != granges_to_features(object@granges)))
    return("The granges slot must resolve to FEATURES slot using granges_to_features")
  TRUE

})
