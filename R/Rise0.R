#' Explicit zeros in sparse matrix.
#'
#' Rise0
#' @keywords internal
#' @description Explicit some implicit zeros in sparse matrix.
#' @param spMtx <dgCMatrix or dgCMatrix coercible>: A sparse matrix.
#' @param indices <numeric>: Vector of positions of the zeros to be explicits (column driven). If NULL and coords NULL all zeros are explicits. (Default NULL)
#' @param coords <data.frame>: A coordinate data frame for zeros to explicit Row index in fisrt column, columns index in second columns. If NULL the indices parameter is used (Default NULL)
#' @return Sparse matrix with some explicit zeros.
#' @examples
#' set.seed(123)
#' spMtx <- as(matrix(floor(runif(7 * 13, 0, 2)), 7, 13), "dgCMatrix")
#' spMtx
#' Rise0(spMtx = spMtx, indices = c(1, 3, 6, 10, 12))
#' Rise0(spMtx = spMtx, coords = data.frame(i = c(1, 5, 3), j = c(1, 2, 3)))
#' Rise0(spMtx = spMtx)
#'
Rise0 <- function(
    spMtx = NULL, indices = NULL, coords = NULL
) {
    if (is.null(coords)) {
        if (is.null(indices)) {
            indices <- which(as.vector(spMtx) == 0)
        }
        coords <- data.frame(
            i = (indices - 1) %% dim(spMtx)[1] + 1,
            j = ((indices - 1) %/% dim(spMtx)[1]) + 1,
            x = 0
        )
    }
    coords$x <- 0
    names(coords) <- c("i", "j", "x")
    mat.dtf <- rbind(MeltSpm(spMtx), coords)
    mat.dtf <- dplyr::arrange(mat.dtf, "j", "i")
    return(Matrix::sparseMatrix(
        i = mat.dtf$i,
        j = mat.dtf$j,
        x = mat.dtf$x,
        dims = dim(spMtx)
    ))
}
