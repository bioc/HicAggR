#' Add a value around a matrix.
#'
#' PadMtx
#' @keywords internal
#' @description Add a value around a matrix.
#' @param mtx <matrix>: Numerical matrix.
#' @param padSize <numeric>: Number of columns or rows to add. (Default 1)
#' @param val <numeric>: Value to add. If Null create mirror of choosen sides. (Default 0)
#' @param side <character>: Side to pad, must be one or some of 'top','bot','right' or 'left'. (Default c('top','bot','right','left') )
#' @return A matrix.
#' @examples
#' mtx <- matrix(seq_len(25), 5, 5)
#' PadMtx(
#'     mtx = mtx,
#'     padSize = 1,
#'     val = 0,
#'     side = c("top", "bot", "right", "left")
#' )
#' PadMtx(
#'     mtx = mtx,
#'     padSize = 1,
#'     val = NULL,
#'     side = c("top", "bot", "right", "left")
#' )
#' PadMtx(
#'     mtx = mtx,
#'     padSize = 1,
#'     val = 0,
#'     side = c("right", "left")
#' )
#' PadMtx(
#'     mtx = mtx,
#'     padSize = 1,
#'     val = 0,
#'     side = c("top")
#' )
#'
PadMtx <- function(
    mtx = NULL, padSize = 1, val = 0,
    side = c("top", "bot", "right", "left")
) {
    if ("top" %in% side) {
        if (!is.null(val)) {
            row.lst <- rep(list(rep(val, dim(mtx)[2])), padSize)
            row.pad <- do.call(rbind, row.lst)
        } else {
            row.pad <- mtx[padSize:1, ]
        }
        mtx <- rbind(row.pad, mtx)
    }
    if ("bot" %in% side) {
        if (!is.null(val)) {
            row.lst <- rep(list(rep(val, dim(mtx)[2])), padSize)
            row.pad <- do.call(rbind, row.lst)
        } else {
            row.pad <- mtx[
                (nrow(mtx) - padSize + 1):nrow(mtx),
            ]
        }
        mtx <- rbind(mtx, row.pad)
    }
    if ("left" %in% side) {
        if (!is.null(val)) {
            col.lst <- rep(list(rep(val, dim(mtx)[1])), padSize)
            col.pad <- do.call(cbind, col.lst)
        } else {
            col.pad <- mtx[, padSize:1]
        }
        mtx <- cbind(col.pad, mtx)
    }
    if ("right" %in% side) {
        if (!is.null(val)) {
            col.lst <- rep(list(rep(val, dim(mtx)[1])), padSize)
            col.pad <- do.call(cbind, col.lst)
        } else {
            col.pad <- mtx[,
                (ncol(mtx) - padSize + 1):ncol(mtx)
            ]
        }
        mtx <- cbind(mtx, col.pad)
    }
    return(mtx)
}
