#' Apply a function over two RLE
#'
#' ReduceRun
#' @keywords internal
#' @description Apply a function on the values over two RLE and return one RLE.
#' @param firstRle <rle or Rle>: First rle.
#' @param secondRle <rle or Rle>>: Second rle.
#' @param reduceMethod <character>: Name of a function to apply e.g paste, sum, mean.
#' @param ... <...>: Other parameter for the reduce function.
#' @return Reduced Rle
#' @examples
#' firstRle <- rle(c("A", "A", "B"))
#' secondRle <- rle(c("A", "B", "B"))
#' ReduceRun(firstRle = firstRle, secondRle = secondRle, reduceMethod = "paste", sep = "_")
#' firstRle <- S4Vectors::Rle(c(1, 2, 3))
#' secondRle <- S4Vectors::Rle(c(5, 5, 5))
#' ReduceRun(firstRle = firstRle, secondRle = secondRle, reduceMethod = "sum")
#'
ReduceRun <- function(
    firstRle, secondRle, reduceMethod = "paste", ...
) {
    if (methods::is(firstRle, "rle")) {
        firstLen.num <- firstRle$length
        firstVal.vec <- firstRle$values
    } else if (methods::is(firstRle, "Rle")) {
        firstLen.num <- S4Vectors::runLength(firstRle)
        firstVal.vec <- S4Vectors::runValue(firstRle)
    }
    if (methods::is(secondRle, "rle")) {
        secondLen.num <- secondRle$length
        secondVal.vec <- secondRle$values
    } else if (methods::is(secondRle, "Rle")) {
        secondLen.num <- S4Vectors::runLength(secondRle)
        secondVal.vec <- S4Vectors::runValue(secondRle)
    }
    newLen.num <- NULL
    newVal.vec <- NULL
    lens.lst <- list(
        firstLen.num = firstLen.num,
        secondLen.num = secondLen.num)
    vals.lst <- list(
        firstVal.vec = firstVal.vec,
        secondVal.vec = secondVal.vec)
    while (
        !is.na(lens.lst[[1]][1]) &&
        !is.na(lens.lst[[2]][1]) &&
        !is.na(vals.lst[[1]][1]) &&
        !is.na(vals.lst[[2]][1])
    ) {
        A.len <- lens.lst[[which.min(c(lens.lst[[1]][1],lens.lst[[2]][1]))]][1]
        A.val <- eval(parse(text = reduceMethod))(
            vals.lst[[1]][1],
            vals.lst[[2]][1],
            ...
        )
        newVal.vec <- c(newVal.vec, A.val)
        newLen.num <- c(newLen.num, A.len)
        lapply(seq_along(lens.lst), function(i) {
            if ((lens.lst[[i]][1] - A.len) > 0) {
                lens.lst[[i]][1] <<- (lens.lst[[i]][1] - A.len)
            } else {
                lens.lst[[i]] <<- lens.lst[[i]][-1]
                vals.lst[[i]] <<- vals.lst[[i]][-1]
            }
        }) |>
            invisible()
    }
    return(S4Vectors::Rle(values = newVal.vec, lengths = newLen.num))
}
