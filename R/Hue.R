#' Hue palette.
#'
#' Hue
#' @description Create an Hue palette.
#' @param paletteLength <numeric>: Color number.
#' @param rotation <numeric>: If positive, rotates clockwise in the color space, reversing if the number is negative. If is NULL compute rotation according to hueRange parameter. (Default NULL)
#' @param hueRange <numeric>: Degree range in color space between 0 and 360. (Default c(0,360))
#' @param saturation <numeric>: Saturation value between 0 and 1. (Default 0.65)
#' @param lightness <numeric>: Lightness value between 0 and 1. (Default 0.65)
#' @param alphaValue <numeric>: Opacity value between 0 and 1. (Default 1)
#' @param alpha <logical>: Whether the alpha layer should be returned. (Default FALSE)
#' @return A vector of color.
#' @examples
#' Hue(paletteLength = 9)
#'
Hue <- function(
    paletteLength = 9, rotation = NULL, hueRange = c(0, 360),
    saturation = 0.65, lightness = 0.65, alphaValue = 1,
    alpha = FALSE
) {
    if (paletteLength > 2) {
        if (is.null(rotation)) {
            if (abs(diff(hueRange)) >=
                180) {
                rotation <- sign(diff(hueRange))
            } else {
                rotation <- -sign(diff(hueRange))
            }
        }
        if (sign(diff(hueRange)) !=
            sign(rotation)) {
            hueRange[which.min(hueRange)] <- 360 +
                hueRange[which.min(hueRange)]
                
        }
        distGap.num <- seq(
            hueRange[1],
            hueRange[2],
            length.out = paletteLength) |>
            diff() |>
            mean() |>
            abs()
        gap.num <- min(
            abs(diff(hueRange)%%360),
            abs(360 - abs(diff(hueRange)%%360))
        )
        if (paletteLength > 1 && gap.num < distGap.num) {
            adjust.num <- (distGap.num - gap.num)/(paletteLength +1) *
                paletteLength
            hueRange[which.max(hueRange)] <- 
                hueRange[which.max(hueRange)] - adjust.num
        }
        hue.lst <- seq(
            hueRange[1],
            hueRange[2],
            length.out = paletteLength)%%360
    } else if (paletteLength == 2) {
        hue.lst <- Hue(
            paletteLength = 5, rotation = rotation,
            hueRange = hueRange, saturation = saturation,
            lightness = lightness, alphaValue = alphaValue,
            alpha = alpha)[c(2, 3)]
    } else if (paletteLength == 1) {
        hue.lst <- Hue(
            paletteLength = 3, rotation = rotation,
            hueRange = hueRange, saturation = saturation,
            lightness = lightness, alphaValue = alphaValue,
            alpha = alpha)[2]
    }
    if (is.numeric(hue.lst)) {
        hue.lst <- lapply(
            hue.lst, function(hue.num) {
                Hsl2Hex(
                    c(
                        hue = hue.num, saturation = saturation,
                        light = lightness, alpha = alphaValue
                    ),
                    alpha = alpha
                )
            }
        ) |>
            unlist()
        return(hue.lst)
    } else {
        return(hue.lst)
    }
}
