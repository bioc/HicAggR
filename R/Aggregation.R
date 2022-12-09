#' Aggregation of matrices list.
#'
#' Aggregation
#' @description Aggregates all the matrices of a list (or two lists in case of differential aggregation) into a single matrix. This function allows to apply different aggregation (average, sum, ...), transformation (rank, percentage, ...) and differential (subtraction, ratio, ...) functions.
#' @param ctrlMatrices <list[matrix]>: The matrices list to aggregate as control.
#' @param matrices <list[matrix]>: The matrices list to aggregate.
#' @param minDist <numeric>: The minimal distance between anchor and bait.
#' @param maxDist <numeric>: The maximal distance between anchor and bait.
#' @param aggFun <function or chracter>: The function use to aggregate each pixel in matrix. If the parameter is a character so:
#' \itemize{
#' \item "50%" or "median" apply the median
#' \item "+" or "sum" apply the sum
#' \item other (Default) apply the mean
#' }
#' @param rm0 <logical>: Whether 0 are replace with NA. (Default FALSE)
#' @param diffFun <function or chracter>: The function use to compute differential. If the parameter is character so:
#' \itemize{
#' \item "-", "substract" or "substraction" apply a substraction
#' \item "/" or "ratio" apply a ratio
#' \item "log2","log2-","log2/" or "log2ratio" apply a log2 on ratio
#' \item other (Default) apply a log2 on 1+ratio
#' }
#' @param transFun <function or chracter>: The function use to transforme or scale values in each submatrix before aggregation. If the parameter is character so:
#' \itemize{
#' \item "quantile" or "qtl" apply function dplyr::ntile(x,500)
#' \item "percentile" or "prct" apply percentile.
#' \item "rank" apply a ranking.
#' \item "zscore" apply a scaling.
#' \item "minmax" apply a HicAggR::MinMaxScale.
#' \item "mu" apply a HicAggR::MeanScale.
#' \item other or NULL don't apply transformation (Default).
#' }
#' @param scaleCorrection <logical>: Whether a correction should be done on the median value take in ane noising area. (Default TRUE)
#' @param correctionArea <list>: Nested list of indice that define a noising area fore correction. List muste contain in first an element "i" (row indices) then an element called "j" (columns indices). If NULL automatically take in upper left part of aggregated matrices. (Default NULL)
#' @param statCompare <logical>: Whether a t.test must be apply to each pxl of the differential aggregated matrix.
#' @param orientate <logical>: Whether matrices must be orientate before the aggregation.
#' @return A matrix
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
#'
#' # Aggregate matrices in one matrix
#' aggreg.mtx <- Aggregation(interactions_Ctrl.mtx_lst)
#'
#' # Differential Aggregation
#' aggregDiff.mtx <- Aggregation(
#'     ctrlMatrices = interactions_Ctrl.mtx_lst,
#'     matrices = interactions_HS.mtx_lst
#' )
#'
Aggregation <- function(
    ctrlMatrices = NULL, matrices = NULL, minDist = NULL,
    maxDist = NULL, transFun = NULL, aggFun = "mean", rm0 = FALSE,
    diffFun = "substraction", scaleCorrection = FALSE,
    correctionArea = NULL, statCompare = FALSE, orientate=TRUE
) {
    # subFunctions
    .PrepareMtxList <- function(
        matrices, minDist = NULL, maxDist = NULL, rm0 = FALSE,
        transFun = NULL, orientate=FALSE
    ) {
        interactions.gni <- attributes(matrices)$interactions
        # Filter on distances
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
        matrices <- matrices[!is.na(names(matrices))]
        # Orientation
        if(orientate){
            matrices  <- OrientateMatrix(matrices)
        }
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
        return(matrices)
    }
    # Put list on correct variable
    if (!is.null(ctrlMatrices) && is.null(matrices)) {
        matrices <- ctrlMatrices
        ctrlMatrices <- NULL
    }
    # Get attributes
    matDim <- attributes(matrices)$matriceDim
    totMtx <- length(matrices)
    attributes.lst <- attributes(matrices)
    if ("names" %in% names(attributes.lst)) {
        attributes.lst <- attributes.lst[-which(
            names(attributes.lst) == "names"
        )]
    }
    # Differential Function
    if (!is.function(diffFun) &&
        !is.null(ctrlMatrices)) {
        diffFun <- dplyr::case_when(
            tolower(diffFun) %in% c("-", "substract", "substraction") ~
                "function(mat.mtx,ctrl.mtx){mat.mtx - ctrl.mtx}",
            tolower(diffFun) %in% c("/", "ratio") ~
                "function(mat.mtx,ctrl.mtx){mat.mtx / ctrl.mtx}",
            tolower(diffFun) %in% c("log2", "log2-", "log2/", "log2ratio") ~
                "function(mat.mtx,ctrl.mtx){log2(mat.mtx)-log2(ctrl.mtx)}",
            TRUE ~
                "function(mat.mtx,ctrl.mtx){log2(mat.mtx+1)-log2(ctrl.mtx+1)}"
        )
        diffFun <- WrapFunction(diffFun)
    }
    # Aggregation Function
    if (!is.function(aggFun)) {
        aggFun <- dplyr::case_when(
            tolower(aggFun) %in% c("50%", "median") ~
                "function(pxl){stats::median(pxl,na.rm=TRUE)}",
            tolower(aggFun) %in% c("+", "sum") ~
                "function(pxl){sum(pxl,na.rm=TRUE)}",
            TRUE ~
                "function(pxl){mean(pxl,na.rm=TRUE,trim=0.01)}"
        )
        aggFun <- WrapFunction(aggFun)
    }
    # Transformation Function
    if (!is.function(transFun) &
        !is.null(transFun)) {
        transFun <- dplyr::case_when(
            tolower(transFun) %in% c("quantile", "qtl") ~
                "function(mat.mtx){
                    matrix(
                        dplyr::ntile(mat.mtx,500),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            tolower(transFun) %in% c("percentile", "prct") ~
                "function(mat.mtx){
                    matrix(
                        rank(mat.mtx,na.last='keep')/
                        length(mat.mtx[!is.na(mat.mtx)]),
                        dim(mat.mtx)[[1]],
                        dim(mat.mtx)[[2]]
                    )
                }",
            tolower(transFun) %in% c("rank") ~
                "function(mat.mtx){
                    matrix(
                        rank(mat.mtx,na.last='keep'),
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
    matrices <- .PrepareMtxList(
        matrices = matrices, minDist = minDist,
        maxDist = maxDist, rm0 = rm0, transFun = transFun
    )
    # Aggregate
    agg.mtx <- apply(
        simplify2array(matrices),
        seq_len(2),
        aggFun
    )
    gc()
    # Differential Case else Return
    if (!is.null(ctrlMatrices)) {
        # Prepare Matrix List
        ctrlMatrices <- .PrepareMtxList(
            matrices = ctrlMatrices, minDist = minDist,
            maxDist = maxDist, rm0 = rm0, transFun = transFun
        )
        # Aggregate
        aggCtrl.mtx <- apply(
            simplify2array(ctrlMatrices),
            seq_len(2),
            aggFun
        )
        gc()
        # Scale mat on Ctrl median
        if (scaleCorrection) {
            if (is.null(correctionArea) ||
            sum(unlist(correctionArea) > matDim)) {
                correctionArea <- list(
                    i = seq_len(round(matDim * 0.3)),
                    j = (matDim - round(matDim * 0.3) + 1):matDim
                )
            }
            correctionValue.num <- stats::median(
                    aggCtrl.mtx[correctionArea[[1]],
                    correctionArea[[2]]]) -
                stats::median(
                    agg.mtx[correctionArea[[1]],
                    correctionArea[[2]]])
            aggCorrected.mtx <- agg.mtx + correctionValue.num
        } else {
            correctionValue.num <- NULL
            aggCorrected.mtx <- NULL
        }
        # Stat compare
        if (statCompare) {
            mtx.nlst <- simplify2array(lapply(matrices, c))
            ctrlMtx.nlst <- simplify2array(lapply(ctrlMatrices, c))
            pval.mtx <- lapply(
                seq_len(dim(mtx.nlst)[[1]]),
                function(ndx) {
                    WT.vec <- mtx.nlst[ndx, ]
                    KD.vec <- ctrlMtx.nlst[ndx, ]
                    WT.vec <- unlist(WT.vec[!is.na(WT.vec)])
                    KD.vec <- unlist(KD.vec[!is.na(KD.vec)])
                    if (length(WT.vec) > 10 & length(KD.vec) > 10) {
                        return(stats::t.test(
                            WT.vec,
                            KD.vec,
                            var = FALSE
                        )$p.value)
                    } else {
                        return(NA)
                    }
                }
            )
            gc()
            pval.mtx <- matrix(
                stats::p.adjust(
                    c(pval.mtx),
                    method = "fdr"
                ),
                nrow = matDim,
                ncol = matDim
            )
            pval.mtx[pval.mtx > 0.05] <- NA
            pval.mtx[pval.mtx < 1e-16] <- 1e-16
            pval.mtx <- -log10(pval.mtx)
        }
        # Differential at submatrix and aggregated scale
        diffmatrices <- lapply(
            seq_along(ctrlMatrices),
            function(mtx.ndx) {
                diff.mtx <- diffFun(
                    matrices[[mtx.ndx]],
                    ctrlMatrices[[mtx.ndx]])
                diff.mtx[is.infinite(diff.mtx)] <- NA
                return(as.matrix(diff.mtx))
            }
        )
        aggDelta.mtx <- diffFun(agg.mtx, aggCtrl.mtx)
        if (!is.null(correctionArea)) {
            aggCorrectedDelta.mtx <- diffFun(aggCorrected.mtx, aggCtrl.mtx)
        } else {
            aggCorrectedDelta.mtx <- NULL
        }
        # Aggregation of differential list
        aggDiff.mtx <- apply(
            simplify2array(diffmatrices),
            seq_len(2),
            aggFun
        )
        gc()
        # Filter aggregated differential by pval.mtx
        if (statCompare) {
            aggDiff.vec <- c(aggDiff.mtx)
            aggDiff.vec[is.na(c(pval.mtx))] <- diffFun(1, 1)
            aggDiffPvalFilt.mtx <- matrix(
                aggDiff.vec,
                nrow = matDim,
                ncol = matDim)
        } else {
            pval.mtx <- NULL
            aggDiffPvalFilt.mtx <- NULL
        }
        # Return
        aggDiff.mtx <- AddAttr(
            x = aggDiff.mtx,
            overwrite = TRUE,
            attrs = c(
                totalMatrixNumber = totMtx,
                filteredMatrixNumber = length(matrices),
                minimalDistance = minDist, maximalDistance = maxDist,
                transformationMethod = transFun, aggregationMethod = aggFun,
                differentialMethod = diffFun, zeroRemoved = rm0,
                correctedFact = correctionValue.num,
                correctionArea = correctionArea,
                matrices = list(list(
                    agg = agg.mtx, aggCtrl = aggCtrl.mtx,
                    aggDelta = aggDelta.mtx, aggCorrected = aggCorrected.mtx,
                    aggCorrectedDelta = aggCorrectedDelta.mtx, pVal = pval.mtx,
                    aggDiffPvalFilt = aggDiffPvalFilt.mtx
                )),
                attributes.lst
            )
        )
        return(aggDiff.mtx)
    } else {
        agg.mtx <- AddAttr(
            x = agg.mtx,
            attrs = c(
                totalMatrixNumber = totMtx,
                filteredMatrixNumber = length(matrices),
                minimalDistance = minDist, maximalDistance = maxDist,
                transformationMethod = transFun, aggregationMethod = aggFun,
                zeroRemoved = rm0, attributes.lst
            )
        )
        return(agg.mtx)
    }
}
