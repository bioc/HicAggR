#' Scales values on min-max range.
#'
#' MinMaxScale
#' @keywords internal
#' @description Scale values on min-max range.
#' @param x.num <numeric>: Numerical vector.
#' @param x_min <numeric>: Minimal value after scaling.
#' @param x_max <numeric>: Maximal value after scaling.
#' @return Scaled numeric vector.
#' @examples
#' set.seed(655213)
#' x.num <- rnorm(500, 500)
#' set.seed(522613)
#' y.num <- rnorm(500, 100)
#' plot(density(x.num), col = "red", xlim = c(min(y.num), max(x.num)))
#' lines(density(y.num), col = "green")
#' plot(
#'     density(MinMaxScale(x.num)),
#'     col = "red",
#'     xlim = c(min(MinMaxScale(y.num)), max(MinMaxScale(x.num)))
#' )
#' lines(density(MinMaxScale(y.num)), col = "green")
#'
MinMaxScale <- function(
    x.num, x_min = (0), x_max = 1
) {
    x_min +
    (
        ((x.num - min(x.num, na.rm = TRUE)) * (x_max - x_min)) /
        (max(x.num, na.rm = TRUE ) - min(x.num, na.rm = TRUE))
    )
}
