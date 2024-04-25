#' Cut HiC map in chunks.
#'
#' CutHiC
#' @description Cut a mega contactMatrix (joint from multiple chromosomic maps)
#'  into a list of contactMatrix.
#' @param megaHic <contactMatrix>: The HiC megamap.
#' @param verbose <logical>: If TRUE,
#' show the progression in console. (Default FALSE)
#' @return A matrices list.
#' @importFrom checkmate assertClass
#' @export
#' @examples
#' data(HiC_Ctrl.cmx_lst)
#' Mega_Ctrl.cmx <- JoinHiC(HiC_Ctrl.cmx_lst)
#' CutHiC(Mega_Ctrl.cmx)
#'
CutHiC <- function(
megaHic, verbose = FALSE
) {
    checkmate::assertClass(
        x = megaHic,  classes = "ContactMatrix",
        null.ok = FALSE)
    hicResolution <- megaHic@metadata$resolution
    mtx.chr <- megaHic@metadata$mtx
    chromSizes <- megaHic@metadata$chromSize
    binnedGenome.grn <- chromSizes |>
        dplyr::pull("length") |>
        stats::setNames(chromSizes$name) |>
        GenomicRanges::tileGenome(
            tilewidth = hicResolution,
            cut.last.tile.in.chrom = TRUE
        )
    GenomeInfoDb::seqlengths(binnedGenome.grn) <- chromSizes$length |>
        stats::setNames(chromSizes$name)
    attributes.tbl <- megaHic@metadata$matricesKind
    chromComb.lst <- attributes.tbl$name
    hic.lst_cmx <- BiocParallel::bplapply(
        BPPARAM = BiocParallel::SerialParam(progressbar = verbose),
        seq_along(chromComb.lst),
        function(ele.ndx) {
            # Chromosomes
            ele.lst <- unlist(strsplit(chromComb.lst[[ele.ndx]], "_"))
            row.regions <- binnedGenome.grn[which(
                as.vector(binnedGenome.grn@seqnames) == ele.lst[[1]]
            )]
            col.regions <- binnedGenome.grn[which(
                as.vector(binnedGenome.grn@seqnames) == ele.lst[[2]]
            )]
            # Extract matrix
            hic.spm <- megaHic@matrix[
                GenomicRanges::findOverlaps(
                    InteractionSet::anchors(megaHic)$row,
                    row.regions
                )@from,
                GenomicRanges::findOverlaps(
                    InteractionSet::anchors(megaHic)$column,
                    col.regions
                )@from
            ]
            hic <- InteractionSet::ContactMatrix(
                hic.spm,
                row.regions,
                col.regions
            )
            # removedCounts
            removedCounts <- NULL
            if (!is.null(megaHic@metadata$removedCounts)) {
                removedCounts <- list(
                    removedCounts = megaHic@metadata$removedCounts[
                        GenomicRanges::findOverlaps(
                            InteractionSet::anchors(megaHic)$row,
                            row.regions
                        )@from,
                        GenomicRanges::findOverlaps(
                            InteractionSet::anchors(megaHic)$column,
                            col.regions
                        )@from
                    ]
                )
            }
            # observed
            observed <- NULL
            if (!is.null(megaHic@metadata$observed)) {
                observed.spm <- megaHic@matrix
                observed.spm@x <- megaHic@metadata$observed
                observed <- list(
                    observed = observed.spm[GenomicRanges::findOverlaps(
                        InteractionSet::anchors(megaHic)$row,
                        row.regions
                    )@from,
                    GenomicRanges::findOverlaps(
                        InteractionSet::anchors(megaHic)$column,
                        col.regions
                    )@from]@x
                )
            }
            # normalizer
            normalizer <- NULL
            if (!is.null(megaHic@metadata$normalizer)) {
                normalizer.spm <- megaHic@matrix
                normalizer.spm@x <- megaHic@metadata$normalizer
                normalizer <- list(
                    normalizer = normalizer.spm[
                        GenomicRanges::findOverlaps(
                            InteractionSet::anchors(megaHic)$row,
                            row.regions
                        )@from,
                        GenomicRanges::findOverlaps(
                            InteractionSet::anchors(megaHic)$column,
                            col.regions
                        )@from
                    ]@x
                )
            }
            # expected
            expected <- NULL
            if (!is.null(megaHic@metadata$expected)) {
                expected.spm <- megaHic@matrix
                expected.spm@x <- megaHic@metadata$expected
                expected <- list(
                    expected = expected.spm[
                        GenomicRanges::findOverlaps(
                            InteractionSet::anchors(megaHic)$row,
                            row.regions
                        )@from,
                        GenomicRanges::findOverlaps(
                            InteractionSet::anchors(megaHic)$column,
                            col.regions
                        )@from
                    ]@x
                )
            }
            # Metadata
            hic@metadata <- dplyr::filter(
                attributes.tbl,
                attributes.tbl$name == chromComb.lst[[ele.ndx]]) |>
                tibble::add_column(resolution = hicResolution) |>
                as.list() |>
                c(removedCounts, observed, normalizer, expected)
            return(hic)
        }
    )
    hic.lst_cmx <- hic.lst_cmx |>
        stats::setNames(chromComb.lst) |>
        AddAttr(
            attrs = list(
                resolution = hicResolution, mtx = mtx.chr,
                chromSize = tibble::as_tibble(chromSizes),
                matricesKind = attributes.tbl
            )
        )
    return(hic.lst_cmx)
}
