#' A Function to Read In Narrow Peaks Files
#'
#' @description This function reads the narrow peaks file (BED6 + 4 format) and turns it into a GRanges object.
#' @param file The file path to the narrow peaks file.
#' @return A GRanges object
#' @export
#' @import tidyverse GenomicRanges
bb_read_narrowpeak <- function(file) {
  data <-
    read_delim(
      file,
      delim = "\t",
      col_names = c(
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "signal",
        "pval",
        "qval",
        "peak_peak"
      )
    )
  res <-
    GenomicRanges::makeGRangesFromDataFrame(df = data, keep.extra.columns = T)
  return(res)
}


#' A Function To Merge Replicate Narrow Peaks GRanges
#'
#'  @description This function merges replicate peaks GRanges objects into 1.
#'  @param peaks_gr A regular list of GRanges objects (Not a GRangesList).
#'  @return A GRanges object
#'  @export
#'  @import IRanges GenomicRanges
bb_merge_narrowpeaks <- function(peaks_gr) {
  res <- Reduce(IRanges::subsetByOverlaps, x = peaks_gr)
  return(res)
}


#' Calculate The Overlap Between Peaks and Promoters
#'
#' @description For a given GRanges object containing peaks, determine how many peaks overlap promoters, how many promoters are overlapped by peaks and the significance of enrichment of query peaks relative to promoters.
#' @param query A GRanges object containing peaks
#' @param tss The tss data base to use.  Must be one of "hg38_tss" or "dr11_tss"
#' @param width The width around the tss to evaluate.  Defaults to 200 bp.
#' @return A list including overlap information and binomal test results.
#' @export
#' @import blaseRdata IRanges GenomicRanges tidyverse
bb_promoter_overlap <-
  function(query,
           tss = c("hg38_tss", "dr11_tss"),
           width = 200) {
    if (tss == "hg38_tss") {
      promoters <- blaseRdata::hg38_tss + width
      genome_length <-
        sum(seqlengths(Hsapiens)[c(paste0("chr", seq(1:22)), "chrX", "chrY")])
    } else if (tss == "dr11_tss") {
      promoters <- blaseRdata::dr11_tss + width
      genome_length <-
        sum(seqlengths(Drerio)[paste0("chr", seq(1:25))])
    } else {
      stop("You must select either 'hg38_tss' or 'dr11_tss'.")
    }
    # calculate the number of promoters that overlap query peaks
    overlap_qp <- IRanges::findOverlaps(query, promoters)
    res_text1 <-
      sprintf(
        "%d of %d promoters are overlapped by a query peak.",
        length(unique(subjectHits(overlap_qp))),
        length(promoters)
      )
    cat(paste0(res_text1, "\n"))

    # calcualte the number of query peaks that overlap promoters
    overlap_pq <- IRanges::findOverlaps(promoters, query)
    res_text2 <-
      sprintf("%d of %d query peaks overlap a promoter.\n",
              length(unique(subjectHits(overlap_pq))),
              length(query))
    cat(paste0(res_text2, "\n"))

    # determine enrichment of query peaks relative to promoters
    ## calculate total promoter length and fraction of total genome length
    promoter_total_length <-
      sum(IRanges::width(IRanges::reduce(promoters)))
    promoter_fraction  <- promoter_total_length / genome_length

    ## use binomial test to find significance.  Maybe this test is overpowered.  Not sure.
    binom_res <-
      binom.test(length(unique(subjectHits(overlap_pq))), length(query), promoter_fraction)


    return(list(res_text1, res_text2, binom_res))

  }

#' Read Bam Files
#'
#' @description This function reads a single sorted, dedupliated paired end bam file and returns either a GRanges object or a GenomicAlignmentPairs object.  The former requires much less memory but a the cost of retaining the outer boundaries of each read.  If read 1 has start S1 and end E1 and read 2 has start S2 and end S2, the Granges object spans S1-E2.
#' @param sortedBam File path to the bam file to load.
#' @param genome One of "hg38" or "danRer11".  This is used to clean up the granges object if necessary.
#' @param return_type Type of object to return.  GRanges is smaller.  GenomicAlignmentPairs retains read pair data.
#' @return An object according to return_type.
#' @export
#' @import tidyverse blaseRdata GenomicRanges GenomicAlignments Rsamtools
bb_read_bam <- function(sortedBam,
                        genome = c("hg38", "danRer11"),
                        return_type = c("GenomicAlignmentPairs", "GRanges")) {
  stopifnot(
    "You have chosen and invalide type to return." = return_type %in% c("GenomicAlignmentPairs", "GRanges")
  )
  if (genome == "hg38") {
    genome_use <- Hsapiens
  } else if (genome == "danRer11") {
    genome_use <- Drerio
  } else {
    stop("You must choose either 'hg38' or 'danRer11' for the genome")
  }
  reads <- GenomicAlignments::readGAlignmentPairs(sortedBam,
                                                  param = ScanBamParam(
                                                    mapqFilter = 1,
                                                    flag = scanBamFlag(isPaired = TRUE,
                                                                       isProperPair = TRUE),
                                                    what = c("qname",
                                                             "mapq",
                                                             "isize"),
                                                    which = GRanges(
                                                      seqnames = names(seqlengths(genome_use)[str_detect(names(seqlengths(genome_use)), pattern = "_|chrM", negate = T)]),
                                                      ranges = IRanges(start = 1,
                                                                       end = seqlengths(genome_use)[str_detect(names(seqlengths(genome_use)), pattern = "_|chrM", negate = T)])
                                                    )
                                                  ))
  if (return_type == "GenomicAlignmentPairs") {
    return(reads)
  } else {
    reads_gr <- granges(reads)
    reads_gr <- buff_granges(reads_gr, gen = genome)
    mcols(reads_gr)$insertSize <- width(reads_gr)
    return(reads_gr)
  }


}

#' A Function to Generate Data For Making MetaPlots
#'
#' @description Use this function to generate data for making TSS enrichment plots or other metafeature plots that are centered on a single genomic locus.  This function returns the data you need for the plot.  Use the tibble element that is returned to plot the enrichment plot and the matrix for the heatmap.  The problem currently is that the binwidths for the enrichment plot need to be smaller than the binwidths for the heatmap to look good.  If you use good binwidths for the enrichment plot, the heatmap will crash.  So either reduce the size of the heatmap matrix before plotting that or rerun the function with a different bin size.  This function allows sample names to be added, so several samples can be column-bound together for comparison.  Each gene is normalized to its own outer flanks so this should account for differences in sequencing depth to some degree.  You also have the option to include all possible TSS in the plot (i.e. including zeros) which you may want to do if comparing several samples.  To do this, set select_hits to FALSE.
#' @param query A GRanges object.  This should be from a bam file so you can plot read coverage across the metagene.
#' @param targets A GRanges object.  The targets you want to plot around.
#' @param select_hits Do you want to plot only the targets that have overlappign query reads?  Defaults to true.
#' @param width The width of the analysis in bp.
#' @param binwidth The binwidth in bp.  Width must be evenly divided by binwidth.
#' @param sample_id An optional sample id if you want to join this matrix up with another one.
#' @return A list including a matrix and a tibble.
#' @export
#' @import tidyverse blaseRdata GenomicRanges IRanges
bb_metafeature <- function(query,
           targets,
           select_hits = TRUE,
           width = 2000,
           binwidth = 10,
           sample_id = NULL) {
    stopifnot("Your width must be evenly divisible by your binwidth." = width %% binwidth == 0)
    stopifnot("Your query must be a GRanges class object." = "GRanges" %in% class(query))
    stopifnot("Your targets must be a GRanges class object." = "GRanges" %in% class(targets))

    # figure the number of bins
    nbins <- 2 * width / binwidth

    # first find which targets (e.g. tss) are overlapped by query
    if (select_hits) {
      hits <- subsetByOverlaps(targets, query)
    } else {
      hits <- targets
    }

    # convert to df
    hits <- data.frame(hits)

    # generate the tiling.
    # For genes on + strand, tiling goes from -width to +width minus 1 bin relative to each target
    # For genes on - strand, tiling goes from +width-1 to -width relative to each target
    tiles <- sapply(1:nrow(hits), function(i)
      if (hits$strand[i] == "+") {
        hits$start[i] + seq(-width, width - binwidth, length.out = nbins)
      } else {
        hits$start[i] + seq(width - binwidth,-width, length.out = nbins)
      })

    # make a granges object from the tiles
    tiles <- GRanges(
      seqnames = Rle(rep(hits$seqnames, each = nbins)),
      ranges = IRanges(start = as.vector(tiles),
                       width = binwidth),
      strand = Rle(rep("*", length(as.vector(
        tiles
      ))))
    )

    # split the query into reads mapping to positive and negative strands
    query_pos <- query[strand(query) == "+"]
    query_neg <- query[strand(query) == "-"]

    # make GPos to mark start and ends and shift to account for Tn5 insertion
    query_pos_start <- GPos(seqnames = seqnames(query_pos),
                            pos = start(query_pos) + 4,
                            strand = strand(query_pos),
                            seqinfo = seqinfo(query_pos))
    query_pos_end <- GPos(seqnames = seqnames(query_pos),
                            pos = end(query_pos) - 5,
                            strand = strand(query_pos),
                            seqinfo = seqinfo(query_pos))
    query_neg_start <- GPos(seqnames = seqnames(query_neg),
                            pos = start(query_neg) +4,
                            strand = strand(query_neg),
                            seqinfo = seqinfo(query_neg))
    query_neg_end <- GPos(seqnames = seqnames(query_neg),
                            pos = end(query_neg) - 5,
                            strand = strand(query_neg),
                            seqinfo = seqinfo(query_neg))


    # count the number of query ranges that overlap each tile
    tile_overlaps <-
      IRanges::countOverlaps(tiles, query_pos_start) +
      IRanges::countOverlaps(tiles, query_pos_end) +
      IRanges::countOverlaps(tiles, query_neg_start) +
      IRanges::countOverlaps(tiles, query_neg_end)

    # convert to matrix with 1 row for each target hit and 1 column for each bin
    tile_overlap_matrix <-
      matrix(
        tile_overlaps,
        nrow = nrow(hits),
        ncol = nbins,
        byrow = TRUE
      )
    #
    rownames(tile_overlap_matrix) <- hits$tx_id
    if (!is.null(sample_id)) {
      colnames(tile_overlap_matrix) <-
        paste0(sample_id, "_bin_", 1:ncol(tile_overlap_matrix))
    } else {
      colnames(tile_overlap_matrix) <-
        paste0("bin_", 1:ncol(tile_overlap_matrix))
    }

    # normalize the matrix to background (Greenleaf method)
    # calculate the number of bins for background
    nbb <- ceiling(0.05 * nbins)

    # calculate the background value
    # background value is the mean of a subset of the full matrix.
    # this subset is 2.5% of the total number of bins (rounded up) on the outer ends of the matrix

    background <-
      mean(tile_overlap_matrix[, c(1:nbb / 2, (ncol(tile_overlap_matrix) - (nbb /
                                                                              2)):ncol(tile_overlap_matrix))])
    # normalize the matrix to background
    normalized_matrix <- tile_overlap_matrix / background

    mean_normalized_counts <- colMeans(normalized_matrix)
    names(mean_normalized_counts) <-
      paste0("bin", 1:length(mean_normalized_counts))

    # convert to tidy df for ggplotting
    tile_overlap_tbl <- enframe(mean_normalized_counts) %>%
      mutate(bin = as.numeric(str_extract(name, "[:digit:]+"))) %>%
      mutate(target_distance = -width + (bin - 1) * binwidth) %>%
      dplyr::select(mean_normalized_count = value, target_distance)
    return_list <- list(normalized_matrix, tile_overlap_tbl)
    names(return_list) <- c("matrix", "tibble")

    return(return_list)


  }
