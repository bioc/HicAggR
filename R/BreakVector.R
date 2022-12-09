#' Cut a vector.
#'
#' BreakVector
#' @keywords internal
#' @description Compute the n+1 breaks of a vector in a linear or density based way with the possibility to fix minimal, center and maximal values.
#' @param x.num <numeric>: Numerical vector.
#' @param x_min <numeric>: Minimal fixed value.
#' @param center.num <numeric>: Center fixed value.
#' @param x_max <numeric>: Maximal fixed value.
#' @param n.num <numeric>: Number of tile (return n.num+1 breaks).
#' @param method <character>: Kind of breaking. "linear" or "density". (Default "linear")
#' @return Numerical vector of breaks.
#' @examples
#' set.seed(31415)
#' BreakVector(x.num = rnorm(100, 50, 200), n.num = 9)
#'
BreakVector <- function(
    x.num = NULL, x_min = NULL,
    center.num = NULL, x_max = NULL,
    n.num = 10, method = "linear"
) {
    n.num <- n.num + 1
    if (method == "linear") {
        x.num <- x.num[which(!is.na(x.num))]
        if (is.null(x_min)) {
            x_min <- min(x.num, na.rm = TRUE)
        }
        if (is.null(x_max)) {
            x_max <- max(x.num, na.rm = TRUE)
        }
        if (is.null(center.num)) {
            breaks.num <- seq(
                x_min, x_max,
                length.out = n.num
            )
        } else if (x_min < center.num & center.num < x_max) {
            breaks.num <- c(
                seq(
                    x_min, center.num,
                    length.out = n.num%/%2 + 1
                ),
                seq(
                    center.num, x_max,
                    length.out = n.num%/%2 + 1
                )
            )
        } else {
            center.num <- stats::median(x.num, na.rm = TRUE)
            breaks.num <- c(
                seq(
                    x_min, center.num,
                    length.out = n.num%/%2 + 1
                ),
                seq(
                    center.num, x_max,
                    length.out = n.num%/%2 + 1
                )
            )
        }
    } else if (method == "density") {
        breaks.num <- stats::quantile(
            x.num, na.rm = TRUE,
            probs = seq(0, 1, length.out = n.num)
        )
    } else {
        stop("Method.chr muste be one of 'linear' or 'density'.\n")
    }
    return(breaks.num[!duplicated(breaks.num)])
}
