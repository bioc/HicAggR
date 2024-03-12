#' Change values of HiC map.
#'
#' SwitchMatrix
#' @description Change values in matrix with observed, balanced,
#'  observed/expected or expected values according to what are be done in hic.
#' @param hicLst <List[ContactMatrix][InteractionSet::ContactMatrix()]>:
#'  The HiC maps list.
#' @param matrixKind <character>: The kind of matrix you want.
#' @return A ContactMatrix list.
#' @export
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
#' HiC_Ctrl_Obs.cmx_lst <- SwitchMatrix(HiC.cmx_lst, matrixKind = "obs")
#' HiC_Ctrl_Norm.cmx_lst <- SwitchMatrix(HiC.cmx_lst, matrixKind = "norm")
#' HiC_Ctrl_oe.cmx_lst <- SwitchMatrix(HiC.cmx_lst, matrixKind = "o/e")
#' HiC_Ctrl_exp.cmx_lst <- SwitchMatrix(HiC.cmx_lst, matrixKind = "exp")
#'
SwitchMatrix <- function(
    hicLst,
    matrixKind = c("obs", "norm", "o/e", "exp")
) {
    if (!(matrixKind %in% c("obs", "norm", "o/e", "exp"))) {
        err.chr <- paste0(
            "matrixKind must be one of",
            "\"obs\", \"norm\", \"o/e\", \"exp\"."
        )
        stop(err.chr)
    }
    if (attributes(hicLst)$mtx != matrixKind) {
        # https://stackoverflow.com/questions/26228625/
        # updating-columns-using-lapply
        # this is to avoid the use <<-
        list2env(lapply(hicLst, function(hicName.chr) {
        hicName.chr@matrix@x <- dplyr::case_when(
            matrixKind == "obs" ~
            (hicName.chr@metadata$observed),
            matrixKind == "norm" ~
            (hicName.chr@metadata$observed *
                hicName.chr@metadata$normalizer),
            matrixKind == "o/e" ~
            (hicName.chr@metadata$observed *
                hicName.chr@metadata$normalizer /
                hicName.chr@metadata$expected),
            matrixKind == "exp" ~
            (hicName.chr@metadata$expected)
        )
        }),envir = .GlobalEnv)
        attributes(hicLst)$mtx <- matrixKind
        return(hicLst)
    } else {
        message("\nhicLst is already ", matrixKind, ".\n")
        # Should return the hicLst as is!!
        return(hicLst)
    }
}
