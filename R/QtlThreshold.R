#' Find threshold for outliers based on quantiles.
#'
#' QtlThreshold
#' @keywords internal
#' @description Find threshold for outliers triming based on quantiles.
#' @param x <numeric>: Numeric vector.
#' @param prctThr <numeric>: Percentage (0-100) threshold. (Default 5)
#' @param tails <character>: Bounds to return, "lower", "upper" or "both". (Default "both")
#' @return Numerical vector of thresholds values for outliers triming.
#' @examples
#' set.seed(1111)
#' x <- 0:100
#' x <- sort(x)
#' x
#' QtlThreshold(x, prctThr = 5, tails = "lower")
#' QtlThreshold(x, prctThr = 5, tails = "both")
#' QtlThreshold(x, prctThr = 5, tails = "upper")
#'
QtlThreshold <- function(
    x = NULL, prctThr = 5, tails = "both"
) {
    probs.num <- dplyr::case_when(
        tails == "both" ~ c(prctThr / 200, 1 - (prctThr / 200)),
        tails == "upper" ~ c(NA, 1 - (prctThr / 100)),
        tails == "lower" ~ c(prctThr / 100, NA)
    )
    return(stats::quantile(x, na.rm = TRUE, probs.num))
}
