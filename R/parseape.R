#' Format a GRanges Object as a Character Vector For Inclusion in Genebank File
#'
#' @description This function takes a GRanges object and returns a character vector.  Only 4 metadata fields from the GRanges object will be included:  locus_tag, type, fwdcolor and revcolor.  Locus_tag must be unique.  This will be checked by the Ape constructor.  This function should mostly be used internally in the construction and FEATURE-setting of instances of the Ape class.
#'
#' @param gr A GRanges object.
#' @return A character vector.
#' @import Biostrings GenomicRanges tidyverse
#' @export
granges_to_features <- function(gr) {
  grlist <-
    list(
      locus_tag = gr@elementMetadata$locus_tag,
      start = gr@ranges@start,
      end = gr@ranges@width + gr@ranges@start - 1,
      type = gr@elementMetadata$type,
      fwdcolor = gr@elementMetadata$fwdcolor,
      revcolor = gr@elementMetadata$revcolor
    )
  res <-
    pmap_chr(
      .l = grlist,
      .f = function(locus_tag,
                    start,
                    end,
                    type,
                    fwdcolor,
                    revcolor) {
        paste0(
          str_pad(
            paste0("     ", type),
            width = 21,
            side = "right",
            pad = " "
          ),
          start,
          "..",
          end,
          "\n                     /locus_tag=",
          "\"",
          locus_tag,
          "\"\n                     /ApEinfo_fwdcolor=",
          "\"",
          fwdcolor,
          "\"\n                     /ApEinfo_revcolor=",
          "\"",
          revcolor,
          "\n"
        )
      }
    )
  res <- c("FEATURES             Location/Qualifiers\n", res)
  return(res)
}

#' Parse a Genebank File and Construct an APE Object
#'
#' @description This is the main function for reading file in genebank/ape/equivalent format and generating an instance of the Ape class.  String manipulations are used to parse the input ape file.  Biostrings and GRanges functions are called to generate DNAStringSet and GRanges objects to store sequence and feature data, respectively.  The Ape constructor function is called internally at the end.
#'
#' @param input_file The genebank/ape file to parse and construct into an instance of the Ape class.
#' @return An Ape object
#' @export
#' @import Biostrings GenomicRanges tidyverse
bb_parseape <- function(input_file) {
  ape <- read_lines(input_file)

  # get first level items
  toplevel <-
    str_extract(ape, "^[:alpha:]*")[which(str_extract(ape, "^[:alpha:]*") %notin% c(""))]
  toplevel_nocomments <-
    str_extract(ape, "^[:alpha:]*")[which(str_extract(ape, "^[:alpha:]*") %notin% c("", "COMMENT"))]

  # get the indices of the top level items
  level0 <- which(str_detect(ape, "^[:alpha:]"))

  # first level chunking
  # make a list of rows that contain metadata
  chunked_ape <- map2(
    .x = level0,
    .y = c(level0[-1], length(ape)),
    .f = function(x, y, obj = ape) {
      start <- x
      end <- y - 1
      chunk <- obj[start:end]
      return(chunk)
    }
  )

  # some of the comments are coming up without the "COMMENTS" prefix.  Unclear why.  Probably we will need to fix this so that we can squish them all back in the same way later on.
  # start by flattening out the comments list items
  flattened_comments <- map(
    .x = chunked_ape,
    .f = function(x) {
      # if(str_detect(x[[1]], "COMMENT", negate = T)) {
      #   x <- NULL
      # }
      string <- x[1]
      test <- str_detect(string = string, pattern = "COMMENT")
      if (test) {
        return(x)
      } else {
        return(NULL)
      }
    }
  ) %>%
    compact() %>%
    flatten()

  # append the "COMMENT" tags if missing
  repaired_comments <-
    map_chr(
      .x = flattened_comments,
      .f = function(x) {
        if (str_detect(x, "^COMMENT     ", negate = T)) {
          x <- str_replace(x, "^[:blank:]*", "")
          x <- paste0("COMMENT     ", x)
        }
        return(x)
      }
    )

  # name all of the elements present.
  chunked_ape_nocomments <- map(
    .x = chunked_ape,
    .f = function(x, levels = toplevel_nocomments) {
      newlevels <- paste0("^", levels)
      object <- x
      map(
        .x = newlevels,
        .f = function(x, obj = object) {
          first <- obj[1]
          test <- str_detect(string =  first, pattern = x)
          if (test) {
            return(object)
          } else {
            return(NULL)
          }
        }
      ) %>%
        compact()
    }

  ) %>%
    compact() %>%
    flatten() %>%
    set_names(toplevel_nocomments)

  # get the indices of teh subelements of the feature section
  feature_index <-
    c(1, which(str_detect(
      chunked_ape_nocomments$FEATURES, "^     [:alpha:]"
    )))

  # chunk them into their own objects and map into a list
  chunked_features <- map2(
    .x = feature_index,
    .y = c(feature_index[-1], length(chunked_ape_nocomments$FEATURES)),
    .f = function(x, y, obj = chunked_ape_nocomments$FEATURES) {
      start <- x
      end <- y - 1
      chunk <- obj[start:end]
      return(chunk)
    }
  )

  # map into character vectors
  chunked_features <-
    map_chr(chunked_features, ~ paste(.x, collapse = "\n"))

  #  removed the unchunked features from the main object
  chunked_ape_nocomments$FEATURES <- NULL

  # add the repaired comments and chunked features back in
  chunked_ape_repaired <-
    c(
      chunked_ape_nocomments,
      "COMMENT" = list(repaired_comments),
      "FEATURES" = list(chunked_features)
    )

  # reorder
  chunked_ape_repaired <- chunked_ape_repaired[unique(toplevel)]

  # get nucleotide sequence
  origin <- paste(chunked_ape_repaired$ORIGIN[-1], collapse = "")
  nucleotides <-
    str_replace_all(string = origin,
                    pattern = "[:blank:]|[:digit:]",
                    replacement = "")

  # make a biostrings object
  dna_biostring <- Biostrings::DNAStringSet(x = nucleotides)
  names(dna_biostring) <- "ape_seq"
  ape_plus <-
    c(chunked_ape_repaired, "dna_biostring" = list(dna_biostring))

  # make granges object from features
  granges_features <- ape_plus$FEATURES


  # reformat granges_features

  granges_df <-
    map_dfr(
      .x = granges_features[-1],
      .f = function(x) {
        firstline <- x[1]
        labelline <- x[which(str_detect(x, "/label="))]
        type <-
          str_extract(string = x, pattern = "^[:blank:]*[:graph:]*")
        type <- str_replace_all(type, "[:blank:]", "")
        start <- str_extract(string = x, pattern = "[:digit:]+")
        end <- str_extract(string = x, pattern = "\\.\\.[:digit:]+")
        end <- str_replace(end, "..", "")
        locus_tag <-
          str_extract(string = x, pattern = "locus_tag=\\\"[:graph:]*")
        locus_tag <-
          str_replace(string = locus_tag,
                      pattern = "locus_tag=",
                      replacement = "")
        locus_tag <-
          str_extract(string = locus_tag, pattern = "[[:alnum:]_()]+")
        fwdcolor <-
          str_extract(string = x, pattern = "fwdcolor=\\\"[:graph:]*")
        fwdcolor <-
          str_replace(string = fwdcolor,
                      pattern = "fwdcolor=",
                      replacement = "")
        fwdcolor <-
          str_extract(string = fwdcolor, pattern = "[[:alnum:]#]+")
        revcolor <-
          str_extract(string = x, pattern = "revcolor=\\\"[:graph:]*")
        revcolor <-
          str_replace(string = revcolor,
                      pattern = "revcolor=",
                      replacement = "")
        revcolor <-
          str_extract(string = revcolor, pattern = "[[:alnum:]#]+")
        ret <-
          tribble(
            ~seqname,
            ~ type,
            ~ start,
            ~ end,
            ~ locus_tag,
            ~ fwdcolor,
            ~ revcolor,
            "ape_seq",
            type,
            as.numeric(start),
            as.numeric(end),
            locus_tag,
            fwdcolor,
            revcolor
          )
        return(ret)
      }
    )

  granges <-
    GenomicRanges::makeGRangesFromDataFrame(df = granges_df, keep.extra.columns = TRUE)
  names(granges) <- granges_df$locus_tag
  ape_plus$granges <- granges

  ape_instance <-
    Ape(
      LOCUS = ape_plus$LOCUS,
      DEFINITION = ape_plus$DEFINITION,
      ACCESSION = ape_plus$ACCESSION,
      VERSION = ape_plus$VERSION,
      SOURCE = ape_plus$SOURCE,
      COMMENT = ape_plus$COMMENT,
      FEATURES = granges_to_features(ape_plus$granges),
      ORIGIN = ape_plus$ORIGIN,
      dna_biostring = ape_plus$dna_biostring,
      granges = ape_plus$granges
    )

  return(ape_instance)
}
