#' Import called loops in .bedpe format to use in HicAggR
#'
#' import_loops
#' @description Imports bedpe file formats as GInteractions object usable
#' to perform submatrix extraction with `ExtractSubmatrix()`
#' @param file_bedpe bedpe file path (Default NULL)
#' @param genomicConstraint <GRanges>: GRanges object of
#' constraint regions. If NULL chromosomes in chromSizes are used as
#' constraints (Default NULL).
#' @param discard_trans <logical>: If TRUE discard loops where
#' anchor and bait are in different genomic constraint elements,
#' either different TADs or chromosomes, if `genomicConstraint = NULL`.
#' (Default FALSE)
#' @param chromSizes <data.frame>: A data.frame containing
#' chromosomes names and lengths in base pairs (see example). (Default NULL)
#' @param binSize <integer>: Bin size in bp - corresponds to
#' matrix resolution.
#' @param minDist <numeric>: Minimal distance between anchors
#' and baits. (Default NULL)
#' @param maxDist <numeric>: Maximal distance between anchors
#' and baits. (Default NULL)
#' @param verbose <logical>: Show the progression in console?
#'  (Default FALSE)
#' @param cores <integer>: Number of cores to use. (Default 1)
#'
#' @return A GInteractions object.
#' @importFrom S4Vectors first second mcols
#' @import InteractionSet
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom data.table as.data.table setkey foverlaps
#' @export
#'
#' @examples
#' 
#' bedpe_path <- system.file("extdata",
#'     "postprocessed_pixels_5000.bedpe",
#'     package = "HicAggR", mustWork = TRUE
#' )
#' loops <- import_loops(
#'     file_bedpe = bedpe_path,
#'     chromSizes = data.frame(c("1","2"),
#'         c(249250621,243199373)),
#'     binSize = 5000, minDist = 105000
#' )
#'
import_loops <- function(
    file_bedpe = NULL,
    genomicConstraint = NULL,
    discard_trans = FALSE,
    chromSizes = NULL,
    binSize = NULL,
    minDist = NULL,
    maxDist = NULL,
    verbose = FALSE,
    cores = 1) {
    loops <- rtracklayer::import(file_bedpe, format = "bedpe")
    if ("ALL" %in% toupper(chromSizes[[1]])){
        chromSizes <- chromSizes[-which(toupper(chromSizes[[1]]) == "ALL"), ]
        if(verbose){
            message("ALL removed from chromSizes")
        }
    }
    anchor_bins <- BinGRanges(
        S4Vectors::first(loops),
        chromSizes = chromSizes,
        cores = cores,
        binSize = binSize,
        na.rm = FALSE
    )
    bait_bins <- BinGRanges(
        S4Vectors::second(loops),
        chromSizes = chromSizes,
        cores = cores,
        binSize = binSize,
        na.rm = FALSE
    )
    # Constraint Informations
    if (is.null(genomicConstraint)) {
        genomicConstraint <- GenomicRanges::GRanges(
            seqnames = chromSizes[[1]],
            ranges = IRanges::IRanges(
                start = rep(1, nrow(chromSizes)),
                end = chromSizes[[2]]
            ),
            strand = "*", name = chromSizes[[1]]
        )
    } else {
        if (is.null(genomicConstraint$name) ||
            length(which(!is.na(genomicConstraint$name))) == 0) {
            genomicConstraint$name <- paste0(
                "Constraint_",
                seq_along(genomicConstraint)
            )
        }
    }
    seqLevelsStyle.chr <- GenomeInfoDb::seqlevelsStyle(genomicConstraint)
    if (length(seqLevelsStyle.chr) > 1) {
        seqLevelsStyle.chr <- seqLevelsStyle.chr[[1]]
        GenomeInfoDb::seqlevelsStyle(genomicConstraint) <- seqLevelsStyle.chr
    }
    genomicConstraint <- GenomeInfoDb::sortSeqlevels(genomicConstraint)
    binned_constraint <- BinGRanges(
        gRange = genomicConstraint,
        chromSizes = chromSizes,
        binSize = binSize,
        verbose = verbose,
        reduceRanges = FALSE,
        cores = cores
    )

    loops_gni <- InteractionSet::GInteractions(
        S4Vectors::first(loops),
        S4Vectors::second(loops)
    )
    loops_gni$anchor.bin <- data.table::foverlaps(
        x = data.table::as.data.table(S4Vectors::first(loops)),
        y = data.table::as.data.table(anchor_bins) |>
            data.table::setkey("seqnames", "start", "end"),
        by.x = c("seqnames", "start", "end"),
        by.y = c("seqnames", "start", "end"), mult="first"
    )$bin
    loops_gni$bait.bin <- data.table::foverlaps(
        x = data.table::as.data.table(S4Vectors::second(loops)),
        y = data.table::as.data.table(bait_bins) |>
            data.table::setkey("seqnames", "start", "end"),
        by.x = c("seqnames", "start", "end"),
        by.y = c("seqnames", "start", "end"), mult="first"
    )$bin
    loops_gni$anchor_constraint <- data.table::foverlaps(
        x = data.table::as.data.table(S4Vectors::first(loops)),
        y = data.table::as.data.table(binned_constraint) |>
            data.table::setkey("seqnames", "start", "end"),
        by.x = c("seqnames", "start", "end"),
        by.y = c("seqnames", "start", "end"),mult="first"
    )$name
    loops_gni$bait_constraint <- data.table::foverlaps(
        x = data.table::as.data.table(S4Vectors::second(loops)),
        y = data.table::as.data.table(binned_constraint) |>
            data.table::setkey("seqnames", "start", "end"),
        by.x = c("seqnames", "start", "end"),
        by.y = c("seqnames", "start", "end"),mult="first"
    )$name

    if (discard_trans) {
        loops_gni <- loops_gni[which(
            S4Vectors::mcols(loops_gni)$anchor_constraint ==
                S4Vectors::mcols(loops_gni)$bait_constraint
        )]
        S4Vectors::mcols(loops_gni)$constraint <-
            S4Vectors::mcols(loops_gni)$anchor_constraint
    } else {
        S4Vectors::mcols(loops_gni)$constraint <-
            S4Vectors::mcols(loops_gni)$anchor_constraint
    }

    loops_gni$distance <- InteractionSet::pairdist(loops_gni)
    if (!is.null(minDist)) {
        loops_gni <- loops_gni[which(
            loops_gni@elementMetadata$distance >= minDist
        )]
    }
    if (!is.null(maxDist)) {
        loops_gni <- loops_gni[which(
            loops_gni@elementMetadata$distance <= maxDist
        )]
    }
    S4Vectors::mcols(loops_gni)$name <- paste0(
        S4Vectors::mcols(loops_gni)$anchor.bin, "_",
        S4Vectors::mcols(loops_gni)$bait.bin
    )
    S4Vectors::mcols(loops_gni)$orientation <- (
        loops_gni == InteractionSet::swapAnchors(loops_gni)
    )
    S4Vectors::mcols(loops_gni)$submatrix.name <-
        S4Vectors::mcols(loops_gni)$name
    S4Vectors::mcols(loops_gni)$submatrix.name[
        !S4Vectors::mcols(loops_gni)$orientation
    ] <- paste0(
        S4Vectors::mcols(loops_gni)$bait.bin[
            !S4Vectors::mcols(loops_gni)$orientation
        ], "_",
        S4Vectors::mcols(loops_gni)$anchor.bin[
            !S4Vectors::mcols(loops_gni)$orientation
        ]
    )
    S4Vectors::mcols(loops_gni)$anchor.name <- paste0(
        "anchor", "_",
        S4Vectors::mcols(loops_gni)$anchor.bin
    )
    S4Vectors::mcols(loops_gni)$bait.name <- paste0(
        "bait", "_",
        S4Vectors::mcols(loops_gni)$bait.bin
    )
    colum_order <- unique(c(
        "name", "constraint", "distance", "orientation", "submatrix.name",
        "anchor.bin", "anchor.name", "bait.bin", "bait.name",
        names(S4Vectors::mcols(loops_gni))[grep(
            x = names(S4Vectors::mcols(loops_gni)),
            pattern = "^anchor.(.*)$"
        )], names(S4Vectors::mcols(loops_gni))[grep(
            x = names(S4Vectors::mcols(loops_gni)),
            pattern = "^bait.(.*)$"
        )]
    ))
    S4Vectors::mcols(loops_gni) <- dplyr::select(
        data.frame(S4Vectors::mcols(loops_gni)),
        dplyr::all_of(colum_order)
    )
    names(loops_gni) <- S4Vectors::mcols(loops_gni)$name
    GenomeInfoDb::seqlevels(loops_gni)<-
    GenomeInfoDb::seqlevels(GenomeInfoDb::Seqinfo(
            seqnames = as.character(chromSizes[[1]]),
            seqlengths = as.numeric(chromSizes[[2]])
        ))
    GenomeInfoDb::seqinfo(loops_gni) <-
        GenomeInfoDb::Seqinfo(
            seqnames = as.character(chromSizes[[1]]),
            seqlengths = as.numeric(chromSizes[[2]])
        )

    return(loops_gni)
}
