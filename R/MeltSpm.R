
#' Coerce matrix in tibble.
#'
#' MeltSpm
#' @keywords internal
#' @description Coerce a sparse matrix M in tibble where columns: i is row index, j is column index and x the value M`[`i,j`]`.
#' @param spMtx <dgCMatrix or dgCMatrix coercible>: A matrix.
#' @return A tibble.
#' @examples
#' i <- c(1, 1, 2, 2, 3, 3, 4, 4, 4, 4)
#' j <- c(1, 4, 2, 5, 1, 4, 2, 3, 4, 5)
#' x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' spMtx <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(5, 5))
#' spMtx
#' meltedMat.tbl <- MeltSpm(spMtx)
#' meltedMat.tbl[order(meltedMat.tbl$i), ]
#'
MeltSpm <- function(
    spMtx = NULL
) {
    if (NotIn("dgCMatrix", class(spMtx))) {
        spMtx <- methods::as(spMtx, "dgCMatrix")
    }
    dp.num <- diff(spMtx@p)
    mat.tbl <- tibble::tibble(
        i = as.integer(spMtx@i + 1),
        j = seq_len(spMtx@Dim[2]) |>
            lapply(function(j.ndx) {
                    rep.num <- dp.num[j.ndx]
                    return(rep(j.ndx, rep.num))
            }) |>
            unlist(),
        x = spMtx@x
    ) |>
        AddAttr(attrs = list(
            matrice.attr = attributes(spMtx)[which(
                NotIn(
                    names(attributes(spMtx)),
                    c("i", "p", "Dimnames", "x", "factors", "class")
                )
            )]
        ))
    return(mat.tbl)
}
