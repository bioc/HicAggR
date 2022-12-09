#' Turns a nested list "inside-out".
#'
#' TransposeList
#' @keywords internal
#' @description Turns a nested list "inside-out".
#' @param nestedList <list[list]>: A nested list to transpose.
#' @return The tranposed nested list.
#' @examples
#' my_lst <- list(
#'     first = list("A1", "B1", "C1"),
#'     second = list("A2", "B2"),
#'     third = list(NULL, "B3")
#' )
#' TransposeList(my_lst)
#'
TransposeList <- function(
    nestedList
) {
    nestedList |>
        lapply(length) |>
        unlist() |>
        max() |>
        seq_len() |>
        lapply(function(newLst.ndx) {
            new.lst <- nestedList |>
                lapply(function(ele.lst) {
                    if (length(ele.lst) >= newLst.ndx) {
                        return(ele.lst[[newLst.ndx]])
                    } else {
                        return(NA)
                    }
                }) |>
                unlist()
            return(new.lst[!is.na(new.lst)])
        })
}
