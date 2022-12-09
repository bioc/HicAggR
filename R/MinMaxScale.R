#' Scales values on min-max range.
#'
#' MinMaxScale
#' @keywords internal
#' @description Scale values on min-max range.
#' @param x <numeric>: Numerical vector.
#' @param x_min <numeric>: Minimal value after scaling.
#' @param x_max <numeric>: Maximal value after scaling.
#' @return Scaled numeric vector.
#' @examples
#' set.seed(655213)
#' x <- rnorm(500, 500)
#' set.seed(522613)
#' y.num <- rnorm(500, 100)
#' plot(density(x), col = "red", xlim = c(min(y.num), max(x)))
#' lines(density(y.num), col = "green")
#' plot(
#'     density(MinMaxScale(x)),
#'     col = "red",
#'     xlim = c(min(MinMaxScale(y.num)), max(MinMaxScale(x)))
#' )
#' lines(density(MinMaxScale(y.num)), col = "green")
#'
MinMaxScale <- function(
    x, x_min = (0), x_max = 1
) {
    x_min +
    (
        ((x - min(x, na.rm = TRUE)) * (x_max - x_min)) /
        (max(x, na.rm = TRUE ) - min(x, na.rm = TRUE))
    )
}
