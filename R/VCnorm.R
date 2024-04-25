#' Compute Vanilla Count Correction.
#'
#' VCnorm
#' @description Compute Vanilla Count or Vanilla Count square root correction
#'  normalization on hic maps.
#' @param hic <contactMatrix>: The HiC maps chunk to normalize.
#' @param qtlTh <numerical>: The threshold quantile below which
#' the bins will be ignored. (Default 0.15)
#' @param vcsqrt <logical>: Whether the square root should be
#' applied. (Default TRUE)
#' @return A matrices list.
#' @importFrom checkmate assertClass
#' @export
#' @examples
#' # Data
#' data(HiC_Ctrl.cmx_lst)
#'
#' HiC_Ctrl_VC.cmx <- VCnorm(HiC_Ctrl.cmx_lst[["2L_2L"]])
#' HiC_Ctrl_VC_SQRT.cmx <- VCnorm(HiC_Ctrl.cmx_lst[["2L_2L"]], vcsqrt = TRUE)
#'
VCnorm <- function(
    hic = NULL,
    qtlTh = 0.15,
    vcsqrt = TRUE
) {
    checkmate::assertClass(
        x = hic,
        classes = "ContactMatrix",
        null.ok = FALSE
    )
    pow.num <- ifelse(vcsqrt, 0.5, 1)
    hic.spm <- hic@matrix
    # Removed Low counts bins
    if (qtlTh) {
        if (hic@metadata$symmetric) {
            rowBias.num <- Matrix::rowSums(hic.spm, na.rm = TRUE) +
                Matrix::colSums(hic.spm, na.rm = TRUE) -
                Matrix::diag(hic.spm)
            colBias.num <- rowBias.num
        } else {
            rowBias.num <- Matrix::rowSums(hic.spm, na.rm = TRUE)
            colBias.num <- Matrix::colSums(hic.spm, na.rm = TRUE)
        }
        row.ndx <- which(
            rowBias.num < stats::quantile(rowBias.num, qtlTh) &
            rowBias.num > 0
        )
        col.ndx <- which(
            colBias.num < stats::quantile(colBias.num, qtlTh) &
            colBias.num > 0
        )
        meltedHic.dtf <- MeltSpm(hic.spm)
        removedHic.dtf <- dplyr::filter(
            meltedHic.dtf,
            meltedHic.dtf$i %in% row.ndx |
            meltedHic.dtf$j %in% col.ndx
        )
        hic@metadata$removedCounts <- Matrix::sparseMatrix(
            i = removedHic.dtf$i,
            j = removedHic.dtf$j,
            x = removedHic.dtf$x,
            dims = dim(hic.spm)
        )
        hic.dtf <- dplyr::filter(
            meltedHic.dtf,
            NotIn(meltedHic.dtf$i,row.ndx) &
            NotIn(meltedHic.dtf$j, col.ndx)
        )
        hic.spm <- Matrix::sparseMatrix(
            i = hic.dtf$i,
            j = hic.dtf$j,
            x = hic.dtf$x,
            dims = dim(hic.spm)
        )
    }
    hic@metadata$observed <- hic.spm@x
    # Bias computation
    coords.tbl <- MeltSpm(hic.spm)
    rowNormalizer.num <- Matrix::rowSums(hic.spm, na.rm = TRUE)
    if (hic@metadata$symmetric) {
        if (!is.na(hic@metadata$kind)) {
            rowNormalizer.num <- rowNormalizer.num +
            Matrix::colSums(hic.spm,na.rm = TRUE) -
            Matrix::diag(hic.spm)
        }
        colNormalizer.num <- rowNormalizer.num
    } else {
        colNormalizer.num <- Matrix::colSums(hic.spm, na.rm = TRUE)
    }
    hic@metadata$normalizer <-
        (
            (
                rowNormalizer.num[coords.tbl$i] *
                colNormalizer.num[coords.tbl$j]
            )^(-pow.num)
        ) *
        sum(hic.spm@x) /
        sum(
            hic.spm@x *
            (
                (
                    rowNormalizer.num[coords.tbl$i] *
                    colNormalizer.num[coords.tbl$j]
                )^(-pow.num)
            )
        )
    X.num <- hic.spm@x *
        (1 / (rowNormalizer.num[coords.tbl$i]^pow.num)) *
        (1 / (colNormalizer.num[coords.tbl$j]^pow.num))
    hic.spm@x <- X.num * sum(hic.spm@x) / sum(X.num)
    hic@matrix <- hic.spm
    hic@metadata$mtx <- "norm"
    return(hic)
}
