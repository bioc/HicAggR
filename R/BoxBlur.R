#' Blur a matrix.
#'
#' BoxBlur
#' @keywords internal
#' @description Blur a matrix with a one dimensional kernel.
#' @param mtx <matrix>: Numerical matrix.
#' @param boxKernel <numeric>: The numerical vector for kernel. If NULL apply a GaussBox (see 'GaussBox' function) (Default NULL)
#' @param kernSize <numeric>: If boxKernel is NULL, size of kernel for 'GaussBox' function. (Default NULL)
#' @param stdev <numeric>: If boxKernel is NULL, standard deviation parameter for 'GaussBox' function. (Default NULL)
#' @return Blurred matrix.
#' @examples
#' set.seed(981643)
#' mtx <- rnorm(10000, 50, 10)**3 |> matrix(100, 100)
#' heatmap(mtx, Rowv = NA, Colv = NA)
#' heatmap(BoxBlur(mtx), Rowv = NA, Colv = NA)
#'
BoxBlur <- function(
    mtx, boxKernel = NULL, kernSize = NULL,
    stdev = 1
) {
    if (is.null(boxKernel)) {
        boxKernel <- GaussBox(
            stdev = stdev, kernScale = "1",
            kernSize = kernSize
        )
    }
    pad.num <- (length(boxKernel) - 1)/2
    mtx <- PadMtx(
        mtx = mtx, padSize = pad.num,
        val = NULL, side = c("top", "bot", "right", "left")
    )
    matVsmth.mtx2 <- sapply(
        ((1 + pad.num):(dim(mtx)[2] - pad.num)),
        function(j) {
            (t(mtx[, (j - pad.num):(j + pad.num)]) *
                boxKernel) |>
                apply(2, Plus)
        }
    )
    matHsmth.mtx2 <- t(sapply(
        ((1 + pad.num):(dim(matVsmth.mtx2)[1] - pad.num)),
        function(i) {
            (matVsmth.mtx2[(i - pad.num):(i + pad.num), ] *
                boxKernel) |>
                apply(2, Plus) |>
                t()
        }
))
    indices <- which(matHsmth.mtx2 == 0)
    if (length(indices)) {
        matHsmth.mtx2 <- Rise0(matHsmth.mtx2, indices = indices)
    }
    return(matHsmth.mtx2)
}
