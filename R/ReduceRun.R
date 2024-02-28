#' Apply a function over two RLE
#'
#' ReduceRun
#' @keywords internal
#' @description Apply a function on the values over two RLE and return one RLE.
#' @param firstRle <rle or Rle>: First rle.
#' @param secondRle <rle or Rle>>: Second rle.
#' @param reduceMethod <character>: Name of a function to apply
#'  e.g paste, sum, mean.
#' @param ... <...>: Other parameter for the reduce function.
#' @return Reduced Rle
#' @examples
#' firstRle <- rle(c("A", "A", "B"))
#' secondRle <- rle(c("A", "B", "B"))
#' ReduceRun(firstRle = firstRle, secondRle = secondRle,
#'  reduceMethod = "paste", sep = "_")
#' firstRle <- S4Vectors::Rle(c(1, 2, 3))
#' secondRle <- S4Vectors::Rle(c(5, 5, 5))
#' ReduceRun(firstRle = firstRle, secondRle = secondRle, reduceMethod = "sum")
#'
ReduceRun <- function(
    firstRle, secondRle, reduceMethod = "paste", ...
) {
    if (methods::is(firstRle, "rle")) {
        firstRle <- S4Vectors::Rle(
        values = firstRle$values,
        lengths = firstRle$lengths)
    }
    if(is.numeric(S4Vectors::runValue(firstRle))) {
        firstVal.vec <- as.numeric(firstRle)
    } else if (is.character(S4Vectors::runValue(firstRle))) {
        firstVal.vec <- as.character(firstRle)
    }
    if (methods::is(secondRle, "rle")) {
        secondRle <- S4Vectors::Rle(
        values = secondRle$values,
        lengths = secondRle$lengths)
    }
    if (is.numeric(S4Vectors::runValue(secondRle))) {
        secondVal.vec <- as.numeric(secondRle)
    }else if (is.character(S4Vectors::runValue(secondRle))) {
        secondVal.vec <- as.character(secondRle)
    }
    vals.lst <- list(
        firstVal.vec = firstVal.vec,
        secondVal.vec = secondVal.vec
    )

    # These changes were made to simplify the code and remove super-assignment
    # Trials with list2env were consistent failures...

    if(is.character(firstVal.vec)){
        # walking across the first values vector is risky,
        # but in both uses of this function (in ExtractSubmatrix),
        # first and second values should have the same length
        newVal.vec <- vapply(seq_along(firstVal.vec), 
            FUN = function(x){eval(parse(text = reduceMethod))(
                vals.lst[[1]][x],
                vals.lst[[2]][x],
        ...
        )}, FUN.VALUE = character(1))
    } else if(is.numeric(firstVal.vec)){
        newVal.vec <- vapply(seq_along(firstVal.vec), 
            FUN = function(x){eval(parse(text = reduceMethod))(
                vals.lst[[1]][x],
                vals.lst[[2]][x],
        ...
        )}, FUN.VALUE = numeric(1))
    }
    return(S4Vectors::Rle(newVal.vec))
}
