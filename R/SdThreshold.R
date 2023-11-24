#' Find threshold for outliers based on sd.
#'
#' SdThreshold
#' @keywords internal
#' @description Find threshold for outliers triming based on standard deviation.
#' @param x <numeric>: numeric vector.
#' @param sdThr <numeric>: number of standard deviation. (Default 3)
#' @param tails <character>: bounds to return, "lower", "upper" or "both". (Default "both")
#' @return numerical vector of thresholds values for outliers triming
#' @examples
#' set.seed(1111)
#' x <- rnorm(1000)
#' x <- sort(x)
#' x
#' SdThreshold(x, sdThr = 2, tails = "lower")
#' SdThreshold(x, sdThr = 2, tails = "both")
#' SdThreshold(x, sdThr = 2, tails = "upper")
#'
SdThreshold <- function(
    x = NULL, sdThr = 3, tails = "both"
) {
    mu.num <- mean(x, na.rm = TRUE)
    sdev.num <- stats::sd(x, na.rm = TRUE)
    upper.num <- mu.num + (sdThr * sdev.num)
    lower.num <- mu.num - (sdThr * sdev.num)
    thr <- dplyr::case_when(
        tails == "both" ~ c(lower.num, upper.num),
        tails == "upper" ~ c(NA, upper.num),
        tails == "lower" ~ c(lower.num, NA)
    )
    return(thr)
}
