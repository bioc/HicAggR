#' Remove outliers.
#'
#' TrimOutliers
#' @keywords internal
#' @description Replace values of a numerical vector that are below a minimal
#'  thresholds and/or above maximal thresholds.
#' @param x <numeric>: Numeric vector.
#' @param thr <numeric>: Numeric vector of length 2. first value is minimal
#'  threshold, second value maximal threshold (Default find threshold based
#'  on standarrd deviation. see `SdThreshold` function)
#' @param clip <logical>: If TRUE the value out of bounds are replace with
#'  threshodls values. If FALSE the Values out of bound are replace with NA
#'  (Default FALSE).
#' @return Trimed Numerical vector.
#' @examples
#' set.seed(1111)
#' x <- rnorm(1000)
#' x <- sort(x)
#' x[990:1000]
#' SdThreshold(x)
#' TrimOutliers(x)[990:1000]
#' TrimOutliers(x, clip = TRUE)[990:1000]
#'
TrimOutliers <- function(
    x,
    thr = SdThreshold(x),
    clip = FALSE
) {
    if (clip) {
        x[which(x > thr[2])] <- thr[2]
        x[which(x < thr[1])] <- thr[1]
    } else {
        x[which(x > thr[2] | thr[1] > x)] <- NA
    }
    return(x)
}
