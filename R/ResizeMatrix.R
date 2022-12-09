#' Resize a matrix
#'
#' ResizeMatrix
#' @keywords internal
#' @description Resize a numericam matrix in new dimension.
#' @param mtx <matrix>: A numerical matrix to resize.
#' @param newDim <integer>: The number of rows and cols in resized matrix.
#' @return Resized matrix.
#' @examples
#' mtx <- matrix(0, 11, 11)
#' mtx[which(as.logical(seq_len(11 * 11) %% 2))] <- seq_len(ceiling((11 * 11) / 2))
#' mtx[2, ] <- 100
#' mtx[, 7] <- 200
#' mtx
#' ResizeMatrix(mtx = mtx, newDim = c(7, 7))
#' ResizeMatrix(mtx = mtx, newDim = c(13, 13))
ResizeMatrix <- function(
    mtx, newDim = dim(mtx)
) {
    # Rescaling
    newCoord.mtx <- as.matrix(
        expand.grid(seq_len(newDim[1]),
        seq_len(newDim[2]))
    )
    rescaleCol.ndx <- MinMaxScale(newCoord.mtx[, 1], 1, dim(mtx)[1])
    rescaleRow.ndx <- MinMaxScale(newCoord.mtx[, 2], 1, dim(mtx)[2])
    # Interpolation
    col.ndx <- floor(rescaleCol.ndx)
    row.ndx <- floor(rescaleRow.ndx)
    xGap.num <- rescaleCol.ndx - col.ndx
    yGap.num <- rescaleRow.ndx - row.ndx
    xGap.num[col.ndx == dim(mtx)[1]] <- 1
    yGap.num[row.ndx == dim(mtx)[2]] <- 1
    col.ndx[col.ndx == dim(mtx)[1]] <- dim(mtx)[1] - 1
    row.ndx[row.ndx == dim(mtx)[2]] <- dim(mtx)[2] - 1
    # Output
    resizedMatrice.mtx <- matrix(
        NA,
        nrow = newDim[1],
        ncol = newDim[2])
    resizedMatrice.mtx[newCoord.mtx] <- mtx[cbind(col.ndx, row.ndx)] *
        (1 - xGap.num) *
        (1 - yGap.num) + mtx[cbind(col.ndx + 1, row.ndx)] *
        xGap.num *
        (1 - yGap.num) + mtx[cbind(col.ndx, row.ndx + 1)] *
        (1 - xGap.num) *
        yGap.num + mtx[cbind(col.ndx + 1, row.ndx + 1)] *
        xGap.num *
        yGap.num
    return(resizedMatrice.mtx)
}
