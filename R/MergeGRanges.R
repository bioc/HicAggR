#' Merge GRanges.
#'
#' MergeGRanges
#' @description Merge GRanges or a list of GRanges
#' @param ... <GRanges or GRangesList or list[GRanges]>:
#' Some GRanges or a list of GRanges or a GRangesList.
#' @param sortRanges <logical>: Whether the result should be
#' sorted. (Default FALSE)
#' @param reduceRanges <logical>: Whether the result should be
#' reduced. See GenomicRanges::reduce for more details. (Default FALSE)
#' @return a GRange object.
#' @export
#' @examples
#' GRange_1.grn <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr1", "chr2", "chr1"), c(1, 3, 1)),
#'     ranges = IRanges::IRanges(101:105, end = 111:115,
#'         names = letters[seq_len(5)]),
#'     strand = S4Vectors::Rle(BiocGenerics::strand(c("-", "+", "*", "+")),
#'         c(1, 1, 2, 1)),
#'     score = seq_len(5)
#' )
#' GRange_2.grn <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr1", "chr3"), c(1, 4)),
#'     ranges = IRanges::IRanges(106:110, end = 116:120, names = letters[6:10]),
#'     strand = S4Vectors::Rle(BiocGenerics::strand(c("*", "+", "-")),
#'         c(2, 1, 2)),
#'     score = 6:10
#' )
#' GRange_1.grn
#' GRange_2.grn
#' MergeGRanges(GRange_1.grn, GRange_2.grn)
#' GRange.lst <- list(GRange_1.grn, GRange_2.grn)
#' MergeGRanges(GRange.lst)
#' MergeGRanges(GRange.lst, reduceRanges = TRUE)
#'
MergeGRanges <- function(
    ...,
    sortRanges = FALSE,
    reduceRanges = FALSE
) {
    mergedGrange.grn <- unlist(GenomicRanges::GRangesList(...))
    if (sortRanges) {
        mergedGrange.grn <- sort(mergedGrange.grn)
    }
    if (reduceRanges) {
        mergedGrange.grn <- GenomicRanges::reduce(mergedGrange.grn)
    }
    return(mergedGrange.grn)
}
