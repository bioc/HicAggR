#' Cut a vector.
#'
#' BreakVector
#' @keywords internal
#' @description Compute the n+1 breaks of a vector in a linear or density based way with the possibility to fix minimal, center and maximal values.
#' @param x <numeric>: Numerical vector.
#' @param x_min <numeric>: Minimal fixed value.
#' @param center <numeric>: Center fixed value.
#' @param x_max <numeric>: Maximal fixed value.
#' @param n <numeric>: Number of tile (return n+1 breaks).
#' @param method <character>: Kind of breaking. "linear" or "density". (Default "linear")
#' @return Numerical vector of breaks.
#' @examples
#' set.seed(31415)
#' BreakVector(x = rnorm(100, 50, 200), n = 9)
#'
BreakVector <- function(
    x = NULL, x_min = NULL,
    center = NULL, x_max = NULL,
    n = 10, method = "linear"
) {
    n <- n + 1
    if (method == "linear") {
        x <- x[which(!is.na(x))]
        if (is.null(x_min)) {
            x_min <- min(x, na.rm = TRUE)
        }
        if (is.null(x_max)) {
            x_max <- max(x, na.rm = TRUE)
        }
        if (is.null(center)) {
            breaks.num <- seq(
                x_min, x_max,
                length.out = n
            )
        } else if (x_min < center & center < x_max) {
            breaks.num <- c(
                seq(
                    x_min, center,
                    length.out = n%/%2 + 1
                ),
                seq(
                    center, x_max,
                    length.out = n%/%2 + 1
                )
            )
        } else {
            center <- stats::median(x, na.rm = TRUE)
            breaks.num <- c(
                seq(
                    x_min, center,
                    length.out = n%/%2 + 1
                ),
                seq(
                    center, x_max,
                    length.out = n%/%2 + 1
                )
            )
        }
    } else if (method == "density") {
        breaks.num <- stats::quantile(
            x, na.rm = TRUE,
            probs = seq(0, 1, length.out = n)
        )
    } else {
        stop("Method.chr muste be one of 'linear' or 'density'.\n")
    }
    return(breaks.num[!duplicated(breaks.num)])
}
