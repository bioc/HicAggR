#' Scale values by mean.
#'
#' MeanScale
#' @keywords internal
#' @description Scale values with mean.
#' @param x <numeric>: Numerical vector.
#' @return Scaled numeric vector.
#' @examples
#' set.seed(655213)
#' x <- rnorm(500, 500)
#' set.seed(522613)
#' y.num <- rnorm(500, 100)
#' plot(density(x), col = "red", xlim = c(min(y.num), max(x)))
#' lines(density(y.num), col = "green")
#' plot(density(MeanScale(x)), col = "red", xlim = c(min(MeanScale(y.num)), max(MeanScale(x))))
#' lines(density(MeanScale(y.num)), col = "green")
#'
MeanScale <- function(
    x
) {
    (x - mean(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x,na.rm = TRUE))
}
