#' Submatrix or Interactions filtering.
#'
#' FilterInteractions
#' @description Search in a GInteraction object which interactions correspond
#'  to a target list and return a list of index or filter a matrices list
#'  according to target and a selectionFunction.
#' @param matrices <List[matrix]>: The matrices list
#' to filter. If is not `NULL`, the function will return the filtred
#' matrices list, else it will return a list of index.
#' @param genomicInteractions <GInteractions>: The GInteraction
#' object on which to compute the filter.
#' @param targets <List>: A named list that describe the target.
#' @param selectionFun <function>: A function
#' that defines how the target variables must be crossed.
#' (Default intersection of all targets)
#' @return A list of elements index or a filtred matrices list with attributes
#'  updates.
#' @export
#' @importFrom S4Vectors mcols
#' @examples
#' # Data
#' data(Beaf32_Peaks.gnr)
#' data(HiC_Ctrl.cmx_lst)
#'
#' # Index Beaf32
#' Beaf32_Index.gnr <- IndexFeatures(
#'     gRangeList = list(Beaf = Beaf32_Peaks.gnr),
#'     chromSizes  = data.frame(seqnames = c("2L", "2R"),
#'         seqlengths = c(23513712, 25286936)),
#'     binSize    = 100000
#' )
#'
#' # Beaf32 <-> Beaf32 Pairing
#' Beaf_Beaf.gni <- SearchPairs(indexAnchor = Beaf32_Index.gnr)
#' Beaf_Beaf.gni <- Beaf_Beaf.gni[seq_len(2000)] # subset 2000 first for exemple
#'
#' # Matrices extractions center on Beaf32 <-> Beaf32 point interaction
#' interactions_PF.mtx_lst  <- ExtractSubmatrix(
#'     genomicFeature         = Beaf_Beaf.gni,
#'     hicLst        = HiC_Ctrl.cmx_lst,
#'     referencePoint = "pf"
#' )
#'
#' # Create a target
#' targets <- list(
#'     anchor.Beaf.name = c("Beaf32_108", "Beaf32_814"),
#'     distance         = function(dist) {
#'         dist < 300000
#'     }
#' )
#' # We target the Beaf32<->Beaf32 interactions that are less than 300Kb away
#' # and have peak Beaf32_2 and Beaf32_191 as anchors (i.e. left).
#'
#' # Create a selection
#' selectionFun <- function() {
#'     intersect(anchor.Beaf.name, distance)
#' }
#' # We select the Beaf32<->Beaf32 interactions that satisfy both targeting
#' # criteria (intersection).
#'
#' # Filtration on InteractionSet (Beaf32 <-> Beaf32 Pairs)
#' FilterInteractions(
#'     genomicInteractions = Beaf_Beaf.gni,
#'     targets        = targets,
#'     selectionFun     = NULL
#' ) |> str(max.level = 1)
#' # Returns a named list (the names match the targeting criteria).
#' # Each element is an index vector of Beaf32<->Beaf32 interactions
#' # that satisfy the given criteria.
#'
#' # Filtration on Matrices List (Beaf32 <-> Beaf32 Extracted matrices)
#' FilterInteractions(
#'     matrices      = interactions_PF.mtx_lst,
#'     targets        = targets,
#'     selectionFun     = NULL
#' ) |> str(max.level = 1)
#' # Return the same kind of result.
#'
#' # Add the selection on InteractionSet Filtration
#' FilterInteractions(
#'     genomicInteractions = Beaf_Beaf.gni,
#'     targets = targets,
#'     selectionFun = selectionFun
#' ) |> str(max.level = 1)
#' # This return the intersection of the index vector that satisfy both
#' # targeting criteria.
#'
#' # Add the selection on Matrices List Filtration
#' FilterInteractions(
#'     matrices = interactions_PF.mtx_lst,
#'     targets = targets,
#'     selectionFun = selectionFun
#' ) |> str(max.level = 1)
#' # This return the filtred matrices, i.e the matrices for which
#' # the Beaf32<->Beaf32 interactions satisfy both targeting criteria.
#'
#' # Filtration with InteractionsSet as filtration criteria
#' targets <- list(interactions = Beaf_Beaf.gni[seq_len(2)])
#' FilterInteractions(
#'     genomicInteractions = Beaf_Beaf.gni,
#'     targets = targets,
#'     selectionFun = NULL
#' ) |> str(max.level = 1)
#'
#'
#' # Filtration with GRanges as filtration criteria
#' targets <- list(first =
#'     InteractionSet::anchors(Beaf_Beaf.gni)[["first"]][seq_len(2)])
#' FilterInteractions(
#'     genomicInteractions = Beaf_Beaf.gni,
#'     targets = targets,
#'     selectionFun = NULL
#' ) |> str(max.level = 1)
#'
FilterInteractions <- function(
    matrices = NULL, genomicInteractions = NULL, targets = NULL,
    selectionFun = function() {Reduce(intersect, interarctions.ndx_lst)}
) {
    if (!is.null(matrices) &&
        !is.null(attributes(matrices)$interactions)) {
        genomicInteractions <- attributes(matrices)$interactions
    }
    interarctions.ndx_lst <- lapply(
        seq_along(targets),
        function(target.ndx) {
            columnName.chr <- names(targets)[target.ndx]
            if (is.function(targets[[target.ndx]])) {
                lapply(
                    S4Vectors::mcols(genomicInteractions)[,columnName.chr],
                    targets[[target.ndx]]
                ) |>
                    unlist() |>
                    which()
            } else if (is.character(targets[[target.ndx]])) {
                lapply(
                    S4Vectors::mcols(genomicInteractions)[,columnName.chr],
                    function(columnElement) {
                        res.lgk <- intersect(
                        as.character(columnElement),
                        targets[[columnName.chr]]
                        ) |>
                        length() |>
                        as.logical()
                        return(res.lgk)
                    }
                ) |>
                    unlist() |>
                    which()
            } else if (inherits(targets[[target.ndx]], "GRanges")) {
                GenomicRanges::findOverlaps(
                    InteractionSet::anchors(
                        genomicInteractions)[[columnName.chr]],
                    targets[[target.ndx]]
                )@from
            } else if (inherits(targets[[target.ndx]], "GInteractions")) {
                InteractionSet::findOverlaps(
                    genomicInteractions,
                    targets[[target.ndx]]
                )@from
            }
        }
    ) |>
        stats::setNames(names(targets))
    if (length(targets) == 1) {
        interarctions.ndx <- unlist(interarctions.ndx_lst)
    } else if (!is.null(selectionFun)) {
        for (target.ndx in seq_along(interarctions.ndx_lst)) {
            assign(
                names(interarctions.ndx_lst)[target.ndx],
                interarctions.ndx_lst[[target.ndx]],
                envir = parent.frame()
            )
        }
        interarctions.ndx <- selectionFun()
    } else {
        return(interarctions.ndx_lst)
    }
    if (!is.null(matrices)) {
        matrices.filt.lst <- matrices[interarctions.ndx]
        attributes(matrices.filt.lst)$interactions <-
            attributes(matrices)$interactions[interarctions.ndx]
        attributes(matrices.filt.lst)$resolution <-
            attributes(matrices)$resolution
        attributes(matrices.filt.lst)$referencePoint <-
            attributes(matrices)$referencePoint
        attributes(matrices.filt.lst)$matriceDim <-
            attributes(matrices)$matriceDim
        attributes(matrices.filt.lst)$target <- targets
        attributes(matrices.filt.lst)$selection <- selectionFun
        return(matrices.filt.lst)
    } else {
        return(interarctions.ndx)
    }
}
