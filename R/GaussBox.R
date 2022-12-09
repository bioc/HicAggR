#' One dimension Gaussian kernel.
#'
#' GaussBox
#' @keywords internal
#' @description One dimension Gaussian kernel.
#' @param stdev <numeric>: Standard deviation parameter of the gaussian. (Default 1)
#' @param kernSize <numeric>: Kernel size. If NULL size is 1+4*stdev (Default NULL)
#' @param kernScale <character>: Scaling kind of box. If "1" sum of kernel equal 1. If "int" Minimal value of kernel is 1 and all entry are integer. If "none", kernel is not scale. (Default "1")
#' @return numerical vector.
#' @examples
#' GaussBox(stdev = 5, kernScale = "none")
#' GaussBox(kernScale = "1")
#' GaussBox(kernScale = "int")
#'
GaussBox <- function(
    stdev = 1, kernSize = NULL, kernScale = "1"
) {
    if (is.null(kernSize)) {
        kernSize <- 1 + 4 * stdev
    }
    x <- as.vector(
        scale(seq_len(kernSize),
        scale = FALSE,
        center = TRUE)
    )
    box <- lapply(x, function(x) {
        xInterval.num <- seq((x - 0.5), (x + 0.5), by = 0.01) |>
            lapply(function(xi) {
                Gauss(x = xi, stdev = stdev)
            }) |>
            unlist() |>
            mean()
        return(xInterval.num)
    }) |>
        unlist()
    if (kernScale == "1") {
        box <- {
            box / sum(abs(box))
        }
    } else if (kernScale == "int") {
        box <- ceiling(box / box[1])
    }
    return(box)
}
