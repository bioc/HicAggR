#' Compute Iterative Correction.
#'
#' ICEnorm
#' @description Compute Iterative Correction (Vanilla Count) on hic maps.
#' @param hic <contactMatrix>: The HiC maps chunk to normalize.
#' @param qtlTh <numerical>: The threshold quantile below which the bins
#'  will be ignored. (Default 0.15)
#' @param maxIter <numerical>: The maximum iteration number.
#' @return A normalized contactMatrix
#' @export
#' @examples
#' data(HiC_Ctrl.cmx_lst)
#' HiC_Ctrl_ICE.cmx <- ICEnorm(HiC_Ctrl.cmx_lst[['2L_2L']])
#'
ICEnorm <- function(hic, qtlTh = 0.15, maxIter = 50) {
    # Removed Low counts bins
    if (qtlTh) {
        if (hic@metadata$symmetric) {
            rowBias.num <- Matrix::rowSums(hic@matrix, na.rm = TRUE) +
                Matrix::colSums(hic@matrix, na.rm = TRUE) -
                Matrix::diag(hic@matrix)
            colBias.num <- rowBias.num
        } else {
            rowBias.num <- Matrix::rowSums(hic@matrix, na.rm = TRUE)
            colBias.num <- Matrix::colSums(hic@matrix, na.rm = TRUE)
        }
        row.ndx <- which(
            rowBias.num < stats::quantile(rowBias.num, qtlTh) &
            rowBias.num > 0
        )
        col.ndx <- which(
            colBias.num < stats::quantile(colBias.num, qtlTh) &
            colBias.num > 0
        )
        meltedHic.dtf <- MeltSpm(hic@matrix)
        removedHic.dtf <- dplyr::filter(
            meltedHic.dtf,
            meltedHic.dtf$i %in% row.ndx |
            meltedHic.dtf$j %in% col.ndx
        )
        hic@metadata$removedCounts <- Matrix::sparseMatrix(
            i = removedHic.dtf$i,
            j = removedHic.dtf$j,
            x = removedHic.dtf$x,
            dims = dim(hic@matrix)
        )
        hic.dtf <- dplyr::filter(
            meltedHic.dtf,
            NotIn(meltedHic.dtf$i, row.ndx) &
            NotIn(meltedHic.dtf$j, col.ndx)
        )
        hic@matrix <- Matrix::sparseMatrix(
            i = hic.dtf$i,
            j = hic.dtf$j,
            x = hic.dtf$x,
            dims = dim(hic@matrix)
        )
    }
    observed.num <- hic@matrix@x
    bias <- Matrix::rowSums(hic@matrix, na.rm = TRUE) +
        Matrix::colSums(hic@matrix, na.rm = TRUE) -
        Matrix::diag(hic@matrix)
    iter.num <- 1
    fit.lm <- stats::lm(c(
        stats::var(bias[which(bias != 0)]),0) ~ c(iter.num,maxIter)
    )
    slope.num <- stats::coef(fit.lm)[2]
    angle.num <- 90 - tan(slope.num) *
        pi/180
    intercept.num <- stats::coef(fit.lm)[1]
    max_gain.num <- -Inf
    while (iter.num < maxIter) {
        hic <- VCnorm(hic, qtlTh = 0, vcsqrt = FALSE)
        bias <- Matrix::rowSums(hic@matrix, na.rm = TRUE) +
            Matrix::colSums(hic@matrix, na.rm = TRUE) -
            Matrix::diag(hic@matrix)
        vertical_distances.num <- abs(
            slope.num * iter.num +
            intercept.num -
            stats::var(bias[which(bias != 0)])
        )
        perpendicular_distance.num <- sin(angle.num) *
            vertical_distances.num
        gain.num <- vertical_distances.num - perpendicular_distance.num
        if (gain.num > max_gain.num) {
            max_gain.num <- gain.num
            i_max <- iter.num
            hicNorm.cmx <- hic
        } else if (i_max + floor(maxIter * 0.1) <= iter.num) {
            break
        }
        iter.num <- iter.num + 1
    }
    hicNorm.cmx@metadata$normalizer <- hicNorm.cmx@matrix@x/observed.num
    hicNorm.cmx@metadata$observed <- observed.num
    return(hicNorm.cmx)
}
