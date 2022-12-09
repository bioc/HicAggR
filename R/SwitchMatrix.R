#' Change values of HiC map.
#'
#' SwitchMatrix
#' @description Change values in matrix with observed, balanced, observed/expected or expected values according to what are be done in hic.
#' @param hicLst <List[contactMatrix]>: The HiC maps list.
#' @param matrixKind.chr <character>: The kind of matrix you want.
#' @return A contactMatrix list.
#' @examples
#' # Data
#' data(HiC_Ctrl.cmx_lst)
#'
#' # Preprocess HiC
#' HiC.cmx_lst <- HiC_Ctrl.cmx_lst |>
#'     BalanceHiC(
#'         interactionType = "cis",
#'         method = "ICE"
#'     ) |>
#'     OverExpectedHiC()
#'
#' # Switch values in matrix
#' HiC_Ctrl_Obs.cmx_lst <- SwitchMatrix(HiC.cmx_lst, matrixKind.chr = "obs")
#' HiC_Ctrl_Norm.cmx_lst <- SwitchMatrix(HiC.cmx_lst, matrixKind.chr = "norm")
#' HiC_Ctrl_oe.cmx_lst <- SwitchMatrix(HiC.cmx_lst, matrixKind.chr = "o/e")
#' HiC_Ctrl_exp.cmx_lst <- SwitchMatrix(HiC.cmx_lst, matrixKind.chr = "exp")
#'
SwitchMatrix <- function(
    hicLst,
    matrixKind.chr = c("obs", "norm", "o/e", "exp")
) {
    if (!(matrixKind.chr %in% c("obs", "norm", "o/e", "exp"))) {
        err.chr <- paste0(
            "matrixKind.chr must be one of",
            "\"obs\", \"norm\", \"o/e\", \"exp\"."
        )
        stop(err.chr)
    }
    if (attributes(hicLst)$mtx != matrixKind.chr) {
        lapply(names(hicLst), function(hicName.chr) {
            hicLst[[hicName.chr]]@matrix@x <<- dplyr::case_when(
                matrixKind.chr == "obs" ~
                    (hicLst[[hicName.chr]]@metadata$observed),
                matrixKind.chr == "norm" ~
                    (hicLst[[hicName.chr]]@metadata$observed *
                    hicLst[[hicName.chr]]@metadata$normalizer),
                matrixKind.chr == "o/e" ~
                    (hicLst[[hicName.chr]]@metadata$observed *
                    hicLst[[hicName.chr]]@metadata$normalizer /
                    hicLst[[hicName.chr]]@metadata$expected),
                matrixKind.chr == "exp" ~
                    (hicLst[[hicName.chr]]@metadata$expected)
            )
        })
        attributes(hicLst)$mtx <- matrixKind.chr
        return(hicLst)
    } else {
        message("\nhicLst is already ", matrixKind.chr, ".\n")
    }
}
