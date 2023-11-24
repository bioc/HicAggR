#' Add list as attributes.
#'
#' AddAttr
#' @keywords internal
#' @description Add list as attributes to any object with or without overwrite.
#' @param x <any>: An object to which attributes are to be added.
#' @param attrs <list>: A named list of new attributes.
#' @param overwrite <logical>: Whether an overwrite is required on attributes with the same name.(Default FALSE)
#' @return The object with new attributes.
#' @examples
#' x <- seq_len(10)
#' x <- AddAttr(x, list(dim = c(2, 5)))
#' x
#' x <- AddAttr(x, list(dim = c(5, 2)))
#' x
#' x <- AddAttr(x, list(dim = c(5, 2)), overwrite = TRUE)
#' x
#'
AddAttr <- function(
    x = NULL, attrs = NULL,
    overwrite = FALSE
) {
    intersectAttr <- intersect(
        names(attributes(x)),
        names(attrs)
    )
    if (overwrite & length(intersectAttr)) {
        attrs <- c(
            attributes(x)[
                which(names(attributes(x)) != intersectAttr)
            ],
            attrs
        )
    } else if (length(intersectAttr)) {
        attrs <- c(
            attributes(x),
            attrs[
                which(names(attrs) != intersectAttr)
            ]
        )
    } else {
        attrs <- c(
            attributes(x),
            attrs
        )
    }
    attributes(x) <- attrs
    return(x)
}