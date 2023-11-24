#' Configure parallel parameters.
#'
#' MakeParallelParam
#' @keywords internal
#' @description Create BiocParallel parameter according to OS.
#' @param cores <numerical> : An integer to specify the number of cores. (Default 1)
#' @param verbose <logical>: A logical value. If TRUE show the progression in console. (Default TRUE)
#' @return Parrallel parameter according number of cores and OS to use with BiocParallel package.
#' @examples
#' multicoreParam <- MakeParallelParam(2)
#' BiocParallel::bplapply(BPPARAM = multicoreParam, seq_len(3), sqrt)
#'
MakeParallelParam <- function(
    cores = 1, verbose = FALSE
) {
    if (!is.numeric(cores) |
        cores < 2 |
        .Platform$OS.type == "windows"
    ) {
        return(BiocParallel::SerialParam(progressbar = verbose))
    } else {
        return(
            BiocParallel::MulticoreParam(
                workers = cores,
                progressbar = verbose
            )
        )
    }
}
