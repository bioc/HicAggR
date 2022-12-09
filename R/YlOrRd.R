#' YlOrRd palette.
#'
#' YlOrRd
#' @description Create a YlOrRd palette.
#' @param paletteLength.num <numeric>: Color number.
#' @param space <numeric>: A character string; interpolation in RGB or CIE Lab color spaces. See ?grDevices::colorRamp for more details. (Default "rgb")
#' @param interpolationMethod <numeric>: Use spline or linear interpolation. See ?grDevices::colorRamp for more details. (Default "linear")
#' @param bias <numeric>: A positive number. Higher values give more widely spaced colors at the high end. See ?grDevices::colorRamp for more details. (Default 1)
#' @return A vector of color.
#' @examples
#' YlOrRd(9)
#'
YlOrRd <- function(
    paletteLength.num = NULL, space = "rgb",
    interpolationMethod = "linear", bias = 1
) {
    (grDevices::colorRampPalette(
        colors = c(
            "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C",
            "#FC4E2A", "#E31A1C", "#BD0026", "#800026"
        ),
        space = space,
        interpolate = interpolationMethod,
        bias = bias
    ))(paletteLength.num)
}
