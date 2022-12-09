#' Check RGB color format.
#'
#' IsRgb
#' @keywords internal
#' @description Check if a color is in RGB color format.
#' @param colour <character or numeric>: A color.
#' @return A logical.
#' @examples
#' IsRgb("red")
#' IsRgb("#FFFFFF")
#' IsRgb(c(125, 125, 125))
#' IsRgb(c(43.8, 0.873, 0.492))
#'
IsRgb <- function(
    colour = NULL
) {
    logical.bln <- lapply(
            colour,
            function(val) {0 <= val & val <= 255}
        ) |>
        unlist() |>
        sort()
    logical.bln <- logical.bln[[1]]
    return(
        (!IsHsl(colour)) &&
        (class(colour) %in% c("list", "numeric", "integer") &&
        logical.bln)
    )
}
