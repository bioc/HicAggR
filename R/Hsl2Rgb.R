#' Convert HSL to RGB.
#'
#' Hsl2Rgb
#' @keywords internal
#' @description Convert a color in HSl (Hue,Saturation,Light) format to RGB format.
#' @param hslColor <charcater>: A vector of the color's HSL code.
#' @param alpha <logical>: Whether the alpha layer should be returned. (Default FALSE)
#' @return An integer vector of the color's RGB code.
#' @examples
#' Hsl2Rgb(c(43.8, 0.873, 0.492, 0.498), alpha = TRUE)
#'
Hsl2Rgb <- function(
    hslColor = NULL, alpha = FALSE
) {
    if (3 > length(hslColor) | length(hslColor) > 4) {
        err.chr <- paste0(
            "Need 3 or 4 values beetween 0 and 255, first value for hue, ",
            "second for saturation, third for light and last for alpha"
        )
        stop(err.chr)
    } else {
        if (IsHsl(hslColor)) {
            if (length(hslColor) == 3) {
                alphaValue <- 255
            } else {
                alphaValue <- hslColor[4] * 255
                hslColor <- hslColor[seq_len(3)]
            }
            C <- hslColor[2] * (1 - abs(2 * hslColor[3] - 1))
            X <- C * (1 - abs((hslColor[1] / 60) %% 2 - 1))
            m <- hslColor[3] - C / 2
            if (0 <= hslColor[1] & hslColor[1] < 60) {
                rbgColor <- c(
                    red = (C + m) * 255,
                    green = (X + m) * 255,
                    blue = (0 + m) * 255
                )
            } else if (60 <= hslColor[1] & hslColor[1] < 120) {
                rbgColor <- c(
                    red = (X + m) * 255,
                    green = (C + m) * 255,
                    blue = (0 + m) * 255
                )
            } else if (120 <= hslColor[1] & hslColor[1] < 180) {
                rbgColor <- c(
                    red = (0 + m) * 255,
                    green = (C + m) * 255,
                    blue = (X + m) * 255
                )
            } else if (180 <= hslColor[1] & hslColor[1] < 240) {
                rbgColor <- c(
                    red = (0 + m) * 255,
                    green = (X + m) * 255,
                    blue = (C + m) * 255
                )
            } else if (240 <= hslColor[1] & hslColor[1] < 300) {
                rbgColor <- c(
                    red = (X + m) * 255,
                    green = (0 + m) * 255,
                    blue = (C + m) * 255
                )
            } else if (300 <= hslColor[1] & hslColor[1] < 360) {
                rbgColor <- c(
                    red = (C + m) * 255,
                    green = (0 + m) * 255,
                    blue = (X + m) * 255
                )
            }
            if (alpha) {
                rbgColor <- c(rbgColor, alpha = alphaValue)
            }
        } else {
            err.chr <- paste0(
                "Need 3 or 4 values beetween 0 and 255, first value for hue, ",
                "second for saturation, third for light and last for alpha"
            )
            stop(err.chr)
        }
    }
    return(round(rbgColor))
}
