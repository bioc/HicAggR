#' Prepare matrices list for further analysis.
#'
#' PrepareMtxList
#' @description Prepares matrices list for further analysis (eg. Aggregation or GetQuantif). Orientation can be corrected, and per matrix transformation can be performed.
#' @param matrices <list[matrix]>: The matrices list to prepare.
#' @param minDist <numeric>: The minimal distance between anchor and bait.
#' @param maxDist <numeric>: The maximal distance between anchor and bait.
#' @param rm0 <logical>: Whether 0 should be replaced with NA. (Default FALSE)
#' @param transFun <function or chracter>: The function used to transform or scale values in each submatrix before aggregation. The following characters can be submitted:
#' \itemize{
#' \item "quantile" or "qtl" apply function dplyr::ntile(x,500)
#' \item "percentile" or "prct" apply percentile.
#' \item "rank" apply a ranking.
#' \item "zscore" apply a scaling.
#' \item "minmax" apply a HicAggR::MinMaxScale.
#' \item "mu" apply a HicAggR::MeanScale.
#' \item other or NULL don't apply transformation (Default).
#' }
#' @param orientate <logical>: Whether matrices must be orientate before the aggregation.
#' @return A matrix list ready for aggregation of values extraction.
#' @examples
#' # Data
#' data(Beaf32_Peaks.gnr)
#' data(HiC_Ctrl.cmx_lst)
#' data(HiC_HS.cmx_lst)
#'
#' # Index Beaf32
#' Beaf32_Index.gnr <- IndexFeatures(
#'     gRangeList = list(Beaf = Beaf32_Peaks.gnr),
#'     chromSizes = data.frame(seqnames = c("2L", "2R"), seqlengths = c(23513712, 25286936)),
#'     binSize = 100000
#' )
#'
#' # Beaf32 <-> Beaf32 Pairing
#' Beaf_Beaf.gni <- SearchPairs(indexAnchor = Beaf32_Index.gnr)
#' Beaf_Beaf.gni <- Beaf_Beaf.gni[seq_len(2000)] # subset 2000 first for exemple
#'
#' # Matrices extractions center on Beaf32 <-> Beaf32 point interaction
#' interactions_Ctrl.mtx_lst <- ExtractSubmatrix(
#'     genomicFeature = Beaf_Beaf.gni,
#'     hicLst = HiC_Ctrl.cmx_lst,
#'     referencePoint = "pf"
#' )
#' interactions_HS.mtx_lst <- ExtractSubmatrix(
#'     genomicFeature = Beaf_Beaf.gni,
#'     hicLst = HiC_HS.cmx_lst,
#'     referencePoint = "pf"
#' )
#' interactions_Ctrl.mtx_lst <- PrepareMtxList(
#'     matrices = interactions_Ctrl.mtx_lst
#' )
#'
#' # Aggregate matrices in one matrix
#' aggreg.mtx <- Aggregation(interactions_Ctrl.mtx_lst)
#'
#'
#' interactions_HS.mtx_lst <- PrepareMtxList(
#'     matrices = interactions_HS.mtx_lst
#' )
#'
#' # Differential Aggregation
#' aggregDiff.mtx <- Aggregation(
#'     ctrlMatrices = interactions_Ctrl.mtx_lst,
#'     matrices = interactions_HS.mtx_lst
#' )
#'
PrepareMtxList <- function(
    matrices, minDist = NULL, maxDist = NULL, rm0 = FALSE,
    transFun = NULL, orientate=FALSE
) {
    # Get attributes
    attributes.lst <- attributes(matrices)
    # Transformation Function
    if (!is.function(transFun) &
        !is.null(transFun)) {
        transFun <- dplyr::case_when(
            tolower(transFun) %in% c("quantile", "qtl") ~
                "function(mat.mtx){
                    matrix(
                        # vectorize matrix
                        dplyr::ntile(c(mat.mtx),500),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            tolower(transFun) %in% c("percentile", "prct") ~
                "function(mat.mtx){
                    matrix(
                        # vectorize matrix
                        rank(c(mat.mtx),na.last='keep')/
                        length(c(mat.mtx[!is.na(mat.mtx)])),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            tolower(transFun) %in% c("rank") ~
                "function(mat.mtx){
                    matrix(
                        # vectorize matrix
                        rank(c(mat.mtx),na.last='keep'),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            tolower(transFun) %in% c("zscore") ~
                "function(mat.mtx){
                    matrix(
                        scale(c(mat.mtx)),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            tolower(transFun) %in% c("minmax") ~
                "function(mat.mtx){
                    matrix(
                        MinMaxScale(c(mat.mtx)),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            tolower(transFun) %in% c("mu") ~
                "function(mat.mtx){
                    matrix(
                        MeanScale(c(mat.mtx)),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            TRUE ~
                "NULL"
        )
        transFun <- WrapFunction(transFun)
    }
    # Prepare Matrix List
    if (is.null(minDist)) {
        minDist <- NA
    }
    if (is.null(maxDist)) {
        maxDist <- NA
    }
    interactions.gni <- attributes(matrices)$interactions
    # Filter on distances
    totMtx <- length(matrices)
    if (!is.na(minDist)) {
        if (is.character(minDist)) {
            minDist <- GenomicSystem(minDist)
        }
        matrices <- matrices[
            S4Vectors::mcols(interactions.gni)$submatrix.name[which(
                S4Vectors::mcols(interactions.gni)$distance >= minDist
            )]
        ]
    }
    if (!is.na(maxDist)) {
        if (is.character(maxDist)) {
            maxDist <- GenomicSystem(maxDist)
        }
        matrices <- matrices[
            S4Vectors::mcols(interactions.gni)$submatrix.name[which(
                S4Vectors::mcols(interactions.gni)$distance <= maxDist
            )]
        ]
    }
    # Here changed order between the filtering in line 169 with orientation because somehow the filtering removes the attributes. This is a cause for error in OrientateMatrix
    # Orientation
    if(orientate){
        matrices  <- OrientateMatrix(matrices)
    }
    
    matrices <- matrices[!is.na(names(matrices))]
    
    # Convert sparse matrix in dense matrix and convert 0 in NA
    # if rm0 is TRUE
    matrices <- lapply(
        matrices, function(spMtx) {
            mat.mtx <- as.matrix(spMtx)
            if (rm0) {
                mat.mtx[mat.mtx == 0] <- NA
            }
            if (!is.null(transFun)) {
                mat.mtx <- transFun(mat.mtx)
            }
            return(mat.mtx)
        }
    )
    # if transFun is null, it can not be added as attribute to the matrices list object
    if (is.null(transFun)) {
        transFun <- NA
    }
    matrices <- AddAttr(
        x = matrices,
        attrs = c(attributes.lst,c(totalMatrixNumber = totMtx,
            filteredMatrixNumber = length(matrices),
            minimalDistance = minDist, 
            maximalDistance = maxDist,
            transformationMethod = transFun,
            zeroRemoved = rm0
        )
        )
    )
    return(matrices)
}