#' Matrix orientation
#'
#' OrientateMatrix
#' @description Oriente extracted Matrix according to the anchors and bait order. Apply a 180° rotation follow with a transposation on a matrix or on matricies in a list according to the interactions attributes of the list.
#' @param mtx <matrix or List[matrix]>: Matrix or matricies list to oriente
#' @return Oriented matrix or matricies list
#' @examples
#' # Data
#' data(Beaf32_Peaks.gnr)
#' data(HiC_Ctrl.cmx_lst)
#'
#' # Index Beaf32 in TADs domains
#' Beaf32_Index.gnr <- IndexFeatures(
#'     gRangeList = list(Beaf = Beaf32_Peaks.gnr),
#'     chromSizes = data.frame(
#'         seqnames = c("2L", "2R"),
#'         seqlengths = c(23513712, 25286936)
#'     ),
#'     binSize = 100000
#' )
#'
#' # Beaf32 <-> Beaf32 Pairing
#' Beaf_Beaf.gni <- SearchPairs(indexAnchor = Beaf32_Index.gnr)
#' Beaf_Beaf.gni <- Beaf_Beaf.gni[seq_len(2000)] # subset 2000 first for exemple
#'
#' # Matrices extractions center on Beaf32 <-> Beaf32 point interaction
#' interactions_PF.mtx_lst <- ExtractSubmatrix(
#'     genomicFeature = Beaf_Beaf.gni,
#'     hicLst = HiC_Ctrl.cmx_lst,
#'     referencePoint = "pf"
#' )
#'
#' # Matrices Orientation
#' oriented_Interactions_PF.mtx_lst <- OrientateMatrix(interactions_PF.mtx_lst)
#'
OrientateMatrix <- function(
    mtx
) {
    if (is.list(mtx) &&
        !is.null(attributes(mtx)$interactions)
    ) {
        orientedMatrice.mtx <- mtx
        message(
            length(which(!attributes(mtx)$interactions$orientation)),
            " matrices are oriented"
        )
        orientedMatrice.mtx[which(
            !attributes(mtx)$interactions$orientation
        )] <- lapply(
            orientedMatrice.mtx[which(
                !attributes(mtx)$interactions$orientation
            )],
            OrientateMatrix
        )
        orientedMatrice.mtx <- AddAttr(
            x = orientedMatrice.mtx,
            attrs = attributes(mtx)
        )
        attributes(orientedMatrice.mtx)$interactions$orientation <- TRUE
        attributes(orientedMatrice.mtx)$interactions$submatrix.name <-
            attributes(orientedMatrice.mtx)$interactions$name
        names(orientedMatrice.mtx) <-
            attributes(orientedMatrice.mtx)$interactions$name
        return(orientedMatrice.mtx)
    } else {
        return(t(apply(
            as.data.frame(apply(mtx, 1, rev)),
            1,
            rev
        )))
    }
}
