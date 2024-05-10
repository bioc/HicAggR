#' Creates pairs from genomic index.
#'
#' SearchPairs
#' @description Creates pairs of coordinates from indexed anchor and bait
#'  genomic coordinates according to distance constraints.
#' @param indexAnchor <GRanges>: A first indexed GRanges object
#' used as pairs anchor (must be indexed using IndexFeatures()).
#' @param indexBait <GRanges>: A second indexed GRanges object
#' used as pairs bait (must be indexed using IndexFeatures()).
#' If NULL, indexAnchor is used instead (Default NULL)
#' @param minDist <numeric>: Minimal distance between anchors
#' and baits. (Default NULL)
#' @param maxDist <numeric>: Maximal distance between anchors
#' and baits. (Default NULL)
#' @param exclude_duplicates <logical> Should duplicated
#' pairs ("2L:100_2L:150" & "2L:150_2L:100") be removed?
#' (Default: TRUE)
#' @param exclude_self_interactions <logical> Should pairs
#' between the same bin ("2L:100_2L:100") be excluded? (Default: TRUE)
#' @param cores <integer> : Number of cores to use. (Default 1)
#' @param verbose <logical>: Show the progression in console?
#'  (Default FALSE)
#' @return A GInteractions object.
#' @export
#' @importFrom S4Vectors mcols
#' @importFrom dplyr select
#' @examples
#' # Data
#' data(Beaf32_Peaks.gnr)
#'
#' # Index Beaf32
#' Beaf32_Index.gnr <- IndexFeatures(
#'     gRangeList = list(Beaf = Beaf32_Peaks.gnr),
#'     chromSizes = data.frame(seqnames = c("2L", "2R"),
#'         seqlengths = c(23513712, 25286936)),
#'     binSize = 100000
#' )
#'
#' # Beaf32 <-> Beaf32 Pairing
#' Beaf_Beaf.gni <- SearchPairs(indexAnchor = Beaf32_Index.gnr)
#'
SearchPairs <- function(
    indexAnchor = NULL, indexBait = NULL, minDist = NULL,
    maxDist = NULL, exclude_duplicates = TRUE,
    exclude_self_interactions = TRUE,
    verbose = FALSE, cores = 1
) {
    .validGranges(indexAnchor)
    if (is.character(minDist)) {
        minDist <- GenomicSystem(minDist)
    }
    if (is.character(maxDist)) {
        maxDist <- GenomicSystem(maxDist)
    }
    if (is.null(indexBait)) {
        indexBait <- indexAnchor
    }else {
        .validGranges(indexBait)
    }
    checkmate::assertLogical(
        x = c(exclude_duplicates,exclude_self_interactions,verbose),
        any.missing = FALSE, null.ok = FALSE
    )
    checkmate::assertNumeric(
        x = c(minDist, maxDist),
        null.ok = TRUE
    )
    checkmate::assertNumeric(
        x = cores,
        lower = 1, null.ok = FALSE
    )
    commonConstraint.lst <- intersect(
        indexAnchor$constraint,
        indexBait$constraint
    )
    multicoreParam <- MakeParallelParam(
        cores = cores,
        verbose = verbose
    )
    pairs.gni_lst <- BiocParallel::bplapply(
        BPPARAM = multicoreParam,
        seq_along(commonConstraint.lst),
        function(constraint.ndx) {
            commonConstraint.chr <- commonConstraint.lst[[constraint.ndx]]
            subIndexAnchor.gnr <- indexAnchor[which(
                indexAnchor@elementMetadata$constraint ==
                commonConstraint.chr
            )]
            subIndexBait.gnr <- indexBait[which(
                indexBait@elementMetadata$constraint ==
                commonConstraint.chr
            )]
            pairsCombination.dtf <- expand.grid(
                seq_along(subIndexAnchor.gnr),
                seq_along(subIndexBait.gnr)
            )
            subIndexAnchor.gnr[pairsCombination.dtf[, "Var1"]]
            subIndexBait.gnr[pairsCombination.dtf[, "Var2"]]
            subPairs.gni <- InteractionSet::GInteractions(
                subIndexAnchor.gnr[pairsCombination.dtf[,"Var1"]],
                subIndexBait.gnr[pairsCombination.dtf[, "Var2"]]
            )
            subPairs.gni$distance <- InteractionSet::pairdist(subPairs.gni)
            if (!is.null(minDist)) {
                subPairs.gni <- subPairs.gni[which(
                    subPairs.gni@elementMetadata$distance >= minDist
                )]
            }
            if(exclude_self_interactions){
                subPairs.gni <- subPairs.gni[which(
                    subPairs.gni@elementMetadata$distance > 0
                )]
            }
            if (!is.null(maxDist)) {
                subPairs.gni <- subPairs.gni[which(
                    subPairs.gni@elementMetadata$distance <= maxDist
                )]
            }
            return(subPairs.gni)
        }
    )
    pairs.gni <- do.call(c, pairs.gni_lst)
    pairs.gni$anchor2.constraint <- NULL
    S4Vectors::mcols(pairs.gni) <- dplyr::rename(
        data.frame(S4Vectors::mcols(pairs.gni)),
        constraint = "anchor1.constraint"
    )
    names(S4Vectors::mcols(pairs.gni)) <- gsub(
        "anchor1.",
        "anchor.",
        gsub(
            "anchor2.",
            "bait.",
            names(S4Vectors::mcols(pairs.gni))
            )
        )
    S4Vectors::mcols(pairs.gni)$name <- paste0(
        S4Vectors::mcols(pairs.gni)$anchor.bin, "_",
        S4Vectors::mcols(pairs.gni)$bait.bin
    )
    S4Vectors::mcols(pairs.gni)$orientation <- (
        pairs.gni == InteractionSet::swapAnchors(pairs.gni)
    )
    S4Vectors::mcols(pairs.gni)$submatrix.name <- 
        S4Vectors::mcols(pairs.gni)$name
    S4Vectors::mcols(pairs.gni)$submatrix.name[
        !S4Vectors::mcols(pairs.gni)$orientation
    ] <- paste0(
        S4Vectors::mcols(pairs.gni)$bait.bin[
            !S4Vectors::mcols(pairs.gni)$orientation
        ],"_",
        S4Vectors::mcols(pairs.gni)$anchor.bin[
            !S4Vectors::mcols(pairs.gni)$orientation
        ]
    )
    columOrder.chr <- unique(c(
        "name", "constraint", "distance", "orientation", "submatrix.name",
        "anchor.bin", "anchor.name", "bait.bin", "bait.name",
        names(S4Vectors::mcols(pairs.gni))[grep(
            x = names(S4Vectors::mcols(pairs.gni)),
            pattern = "^anchor.(.*)$"
        )], names(S4Vectors::mcols(pairs.gni))[grep(
            x = names(S4Vectors::mcols(pairs.gni)),
            pattern = "^bait.(.*)$"
        )]
    ))
    S4Vectors::mcols(pairs.gni) <- dplyr::select(
        data.frame(S4Vectors::mcols(pairs.gni)),
        columOrder.chr
    )
    if(exclude_duplicates){
        dupSubMatrices <- which(
            duplicated(apply(
                X = (S4Vectors::mcols(pairs.gni) |> 
                as.data.frame() |>
                dplyr::select(c("anchor.bin","bait.bin"))),
                MARGIN = 1,
                FUN = function(x) paste(sort(x),collapse = "-"))))
        if(length(dupSubMatrices) > 0) pairs.gni <- pairs.gni[-dupSubMatrices]
    }
    names(pairs.gni) <- S4Vectors::mcols(pairs.gni)$name
    ## To avoid error when using interactions built 
    ## with grangeslist for both the anchor and bait. eg:
    ## Error in new_Rle(values, lengths) : Rle of type 'NULL' is not supported
    if(is.null(GenomeInfoDb::seqinfo(pairs.gni))){
        GenomeInfoDb::seqinfo(pairs.gni) <- GenomeInfoDb::seqinfo(indexAnchor)
    }
    
    return(pairs.gni)
}
