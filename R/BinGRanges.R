#' Bin a GRanges.
#'
#' BinGRanges
#' @description Bin a GRanges and apply a summary method (e.g: 'mean',
#'  'median', 'sum', 'max, 'min' ...) to a chosen numerical variable of ranges
#'  in the same bin.
#' @param gRange <GRanges>: A GRanges to bin.
#' @param chromSizes <data.frame>: A data.frame where first colum corresponds to
#'  the chromosomes names, and the second column corresponds to the chromosomes
#'  lengths in base pairs.
#' @param binSize <numerical>: Width of the bins.
#' @param method <character>: Name of a summary method as 'mean', 'median',
#'  'sum', 'max, 'min'. (Default 'mean')
#' @param metadataColName <character> : A character vector that specify the
#'  metadata columns of GRanges on which apply the summary method.
#' @param na.rm <logical> : A logical value indicating whether 'NA' values
#'  should be stripped before the computation proceeds. (Default TRUE)
#' @param cores <numerical> : The number of cores. (Default 1)
#' @param reduceRanges <logical> : Whether duplicated Bin must been reduced
#'  with de summary method. (Default TRUE)
#' @param verbose <logical>: If TRUE show the progression in console.
#'  (Default FALSE)
#' @return A binned GRanges.
#' @importFrom stats na.omit
#' @export
#' @examples
#' GRange.gnr <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr1", "chr2"), c(3, 1)),
#'     ranges = IRanges::IRanges(
#'         start = c(1, 201, 251, 1),
#'         end = c(200, 250, 330, 100),
#'         names = letters[seq_len(4)]
#'     ),
#'     strand = S4Vectors::Rle(BiocGenerics::strand(c("*")), 4),
#'     score = c(50, NA, 100, 30)
#' )
#' GRange.gnr
#' BinGRanges(
#'     gRange = GRange.gnr,
#'     chromSizes = data.frame(c("chr1", "chr2"), c(350, 100)),
#'     binSize = 100,
#'     method = "mean",
#'     metadataColName = "score",
#'     na.rm = TRUE
#' )
#'
BinGRanges <- function(
    gRange = NULL, chromSizes = NULL, binSize = NULL,
    method = "mean", metadataColName = NULL, na.rm = TRUE,
    cores = 1, reduceRanges = TRUE, verbose = FALSE
) {
    if (is.null(chromSizes)) {
        seqlengths.lst <- GenomeInfoDb::seqlengths(gRange)
    } else {
        seqlengths.lst <- dplyr::pull(chromSizes, 2) |>
            stats::setNames(dplyr::pull(chromSizes, 1))
        seqlengths.lst <- seqlengths.lst[stats::na.omit(match(
            levels(GenomeInfoDb::seqnames(gRange)@values),
            names(seqlengths.lst)
        ))]
        gRange <- GenomeInfoDb::keepSeqlevels(
            gRange, value = names(seqlengths.lst),
            "coarse"
        )
        GenomeInfoDb::seqlengths(gRange) <- seqlengths.lst
    }
    binnedGenome.gnr <- GenomicRanges::tileGenome(
        seqlengths.lst, tilewidth = binSize, cut.last.tile.in.chrom = TRUE
    )
    ovlp.dtf <- GenomicRanges::findOverlaps(binnedGenome.gnr, gRange)
    binnedGRanges.gnr <- binnedGenome.gnr[ovlp.dtf@from]
    S4Vectors::mcols(binnedGRanges.gnr) <- S4Vectors::mcols(
        gRange[ovlp.dtf@to])
    binnedGRanges.gnr$bin <- paste0(
        GenomeInfoDb::seqnames(binnedGRanges.gnr), ":",
        ceiling(BiocGenerics::start(binnedGRanges.gnr)/binSize)
    )
    dupplicated.lgk <- duplicated(binnedGRanges.gnr$bin)
    dupplicated.id <- binnedGRanges.gnr$bin[dupplicated.lgk]
    if (reduceRanges && length(dupplicated.id)) {
        binnedGRange.tbl <- tibble::tibble(data.frame(binnedGRanges.gnr))
        nodup_binnedGRange.tbl <- dplyr::slice(
            binnedGRange.tbl,
            which(NotIn(binnedGRange.tbl$bin, dupplicated.id))
        )
        dup_binnedGRange.tbl <- dplyr::slice(
            binnedGRange.tbl,
            which(binnedGRange.tbl$bin %in% dupplicated.id)
        )
        dup_binnedGRange.tbl <- dplyr::group_by(
            dup_binnedGRange.tbl,
            bin = dup_binnedGRange.tbl$bin) |>
            tidyr::nest()
        multicoreParam <- MakeParallelParam(
            cores = cores,
            verbose = verbose)
        dup_binnedGRange.lst <- BiocParallel::bplapply(
            BPPARAM = multicoreParam, seq_len(nrow(dup_binnedGRange.tbl)),
            function(row.ndx) {
                rowName.chr <- dup_binnedGRange.tbl$bin[[row.ndx]]
                row <- dup_binnedGRange.tbl$data[[row.ndx]]
                row.lst <- lapply(seq_along(row), function(col.ndx) {
                    col <- dplyr::pull(row, col.ndx)
                    colName.chr <- names(row)[col.ndx]
                    if (is.numeric(col) &
                    colName.chr %in% metadataColName) {
                        return(
                            as.numeric(
                                eval(parse(text = method))(
                                    as.numeric(col),
                                    na.rm = na.rm)
                            )
                        )
                    } else if (length(unique(col)) == 1) {
                        return(unique(col))
                    } else {
                        return(list(col))
                    }
                }) |>
                stats::setNames(names(row))
                row.tbl <- tibble::as_tibble(row.lst) |>
                    tibble::add_column(bin = rowName.chr)
                return(row.tbl)
            }
        )
        dup_binnedGRange.dtf <- BindFillRows(dup_binnedGRange.lst)
        dup_binnedGRange.tbl <- tibble::tibble(dup_binnedGRange.dtf)
        dup_binnedGRange.tbl <- dplyr::mutate(
            dup_binnedGRange.tbl,
            strand = as.factor(dup_binnedGRange.tbl$strand)
        )
        for (colName.chr in names(dup_binnedGRange.tbl)) {
            method <- dplyr::pull(
                dup_binnedGRange.tbl,
                dplyr::all_of(colName.chr)) |>
                class()
            method.fun <- eval(parse(text = paste0("as.", method)))
            nodup_binnedGRange.tbl <- nodup_binnedGRange.tbl |>
                dplyr::mutate(
                    dplyr::across(
                        dplyr::all_of(colName.chr),
                        method.fun
                    )
                )
        }
        binnedGRange.tbl <- dplyr::bind_rows(
            dup_binnedGRange.tbl,
            nodup_binnedGRange.tbl)
        binnedGRanges.gnr <- methods::as(binnedGRange.tbl, "GRanges")
    }
    binnedGRanges.gnr <- sort(binnedGRanges.gnr)
    GenomeInfoDb::seqinfo(binnedGRanges.gnr) <- GenomeInfoDb::seqinfo(
        sort(GenomeInfoDb::sortSeqlevels(gRange))
    )
    return(binnedGRanges.gnr)
}
