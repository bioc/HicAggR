#' Genomic distance bias correction.
#'
#' OverExpectedHiC
#' @description Function that normalises HiC matrices by expected values
#'  computed per genomic distance.
#' @param hicLst <List[ContactMatrix][InteractionSet::ContactMatrix()]>:
#'  The HiC maps list.
#' @param cores <numerical> : Number of cores to be used. (Default 1)
#' @param verbose <logical>: If TRUE show the progression in console.
#'  (Default FALSE)
#' @return A matrices list.
#' @examples
#' # Note: run HicAggR::BalanceHiC before OverExpectedHiC calculation.
#' data(HiC_Ctrl.cmx_lst)
#' OverExpectedHiC(HiC_Ctrl.cmx_lst)
#'
OverExpectedHiC <- function(
    hicLst, verbose = FALSE, cores = 1
) {
    resolution.num <- attributes(hicLst)$resolution
    matricesNames.chr <- names(hicLst)
    multicoreParam <- MakeParallelParam(
        cores = cores,
        verbose = verbose
    )
    expected.lst <- BiocParallel::bplapply(
        BPPARAM = multicoreParam, seq_along(matricesNames.chr),
        function(matrixName.ndx) {
            matrixName.chr <- matricesNames.chr[matrixName.ndx]
            if (hicLst[[matrixName.chr]]@metadata$type == "cis") {
                hic.dtf <- MeltSpm(hicLst[[matrixName.chr]]@matrix)
                hic.dtf <- dplyr::mutate(
                    hic.dtf,
                    distance = 1 + (hic.dtf$j - hic.dtf$i) * resolution.num) |>
                    dplyr::select(c("x", "distance"))
                expected.dtf <- dplyr::group_by(
                    hic.dtf, distance = hic.dtf$distance) |>
                    dplyr::summarise_at(
                        .vars = "x", .funs = list(expected = mean)
                    )
                expected.num <- dplyr::left_join(
                    hic.dtf,
                    expected.dtf,
                    by = "distance"
                )
            } else {
                expected.num <- list(
                    expected = sum(hicLst[[matrixName.chr]]@matrix@x)/
                        Reduce(`*`, hicLst[[matrixName.chr]]@matrix@Dim)
                )
            }
            expected.lst <- list(expected.num) |>
                stats::setNames(matrixName.chr)
            return(expected.lst)
        }
    ) |> do.call(what = "c")
    cisNames.chr <- NULL
    for (matrixName.chr in matricesNames.chr) {
        hicLst[[matrixName.chr]]@metadata$expected <-
            expected.lst[[matrixName.chr]]$expected
        hicLst[[matrixName.chr]]@matrix@x <-
            hicLst[[matrixName.chr]]@matrix@x/
            expected.lst[[matrixName.chr]]$expected
        if (hicLst[[matrixName.chr]]@metadata$type == "cis") {
            cisNames.chr <- c(cisNames.chr, matrixName.chr)
        }
    }
    attributes(hicLst)$mtx <- "o/e"
    expected.dtf <- expected.lst[cisNames.chr] |>
        do.call(what = rbind)
    attributes(hicLst)$expected <- dplyr::group_by(
        expected.dtf,
        distance = expected.dtf$distance) |>
        dplyr::summarise_at(.vars = "expected", .funs = list(expected = mean))
    return(hicLst)
}