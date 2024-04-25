#' Change values of HiC map.
#'
#' SwitchMatrix
#' @description Change values in matrix with observed, balanced,
#'  observed/expected or expected values according to what are be done in hic.
#' @param hicLst <List[ContactMatrix][InteractionSet::ContactMatrix()]>:
#'  The HiC maps list.
#' @param matrixKind <character>: The kind of matrix you want.
#' @return A ContactMatrix list.
#' @importFrom purrr map
#' @importFrom dplyr case_when
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
    .validHicMatrices(matrices = hicLst)
    if (!(matrixKind %in% c("obs", "norm", "o/e", "exp"))) {
        err.chr <- paste0(
            "matrixKind must be one of",
            "\"obs\", \"norm\", \"o/e\", \"exp\"."
        )
        stop(err.chr)
    }
    if (attributes(hicLst)$mtx != matrixKind) {
        ## this is to avoid the use of <<-
        ## replaced lapply by map
        hicLst <- purrr::map(
            .x = hicLst,
            .f =~{
                .x@matrix@x <- dplyr::case_when(
                    matrixKind == "obs" ~
                        (.x@metadata$observed),
                    matrixKind == "norm" ~
                        (.x@metadata$observed *
                        .x@metadata$normalizer),
                    matrixKind == "o/e" ~
                        ((.x@metadata$observed *
                        .x@metadata$normalizer) /
                        .x@metadata$expected),
                    matrixKind == "exp" ~
                        (.x@metadata$expected)
                ); .x}) |>
            `attributes<-`(attributes(hicLst))
        attributes(hicLst)$mtx <- matrixKind
        return(hicLst)
    } else {
        message("\nhicLst is already ", matrixKind, ".\n")
        # Should return the hicLst as is!!
        return(hicLst)
    }
}
