#' Get file extension
#'
#' GetFileExtension
#' @keywords internal
#' @description Give the extension of a file from the path.
#' @param file <character>: The path to the file.
#' @return A character string
#' @examples
#' filePath.pth ="my/path/to/my/file.txt"
#' GetFileExtension(file=filePath.pth)

GetFileExtension <- function(
    file = NULL
) {
    fileName.chr <- GetFileName(file, extension = TRUE) |>
        strsplit(".", fixed = TRUE) |>
        unlist()
    return(fileName.chr[length(fileName.chr)])
}