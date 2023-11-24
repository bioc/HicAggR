#' Check HSL color format.
#'
#' IsHsl
#' @keywords internal
#' @description Check if a color is in HSL color format.
#' @param colour <character or numeric>: A color.
#' @return A logical.
#' @examples
#' IsHsl("red")
#' IsHsl("#FFFFFF")
#' IsHsl(c(125, 125, 125))
#' IsHsl(c(43.8, 0.873, 0.492))
#'
IsHsl <- function(
    colour = NULL
) {
    logical.bln <- lapply(
            colour[-1],
            function(val) {0 <= val & val <= 1 }
        ) |>
        unlist() |>
        sort()
    logical.bln <- logical.bln[[1]]
    return(
        class(colour) %in% c("list", "numeric", "integer") &&
        0 <= colour[[1]] &&
        colour[[1]] < 360 &&
        logical.bln
    )
}
