#' Get file name.
#'
#' GetFileName
#' @keywords internal
#' @description Function as `basename()` with the option to not return the file extension.
#' @param file <character>: The path to the file.
#' @param extension <logical>: Whether the file extension should be returned with the file name. (Default FALSE)
#' @return A character string.
#' @examples
#' filePath.pth ="my/path/to/my/file.txt"
#' GetFileName(file=filePath.pth, extension=FALSE)
#' GetFileName(file=filePath.pth, extension=TRUE)

GetFileName <- function(
    file = NULL, extension = FALSE
) {
    ifelse(extension,
        fileName.str <- basename(file),
        fileName.str <- sub(pattern = "(.*)\\..*$",
        replacement = "\\1", basename(file))
    )
    return(fileName.str)
}