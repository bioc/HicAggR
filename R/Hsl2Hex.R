#' Convert HSL to Hex.
#'
#' Hsl2Hex
#' @keywords internal
#' @description Convert a color in HSL (Hue,Saturation,Light) to hexadecimal format.
#' @param hslColor <charcater>: A vector of the color's HSL code.
#' @param alpha <logical>: Whether the alpha layer should be returned. (Default FALSE)
#' @return A character of the color's hexadecimal code.
#' @examples
#' Hsl2Hex(c(43.8, 0.873, 0.492, 0.498), alpha = TRUE)
#'
Hsl2Hex <- function(
    hslColor = NULL, alpha = FALSE
) {
    color.hex <- Hsl2Rgb(hslColor, alpha) |>
        Rgb2Hex(alpha)
    return(color.hex)
}
