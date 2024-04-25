#' .validHicMatrices
#' check the validity of hic matrices
#'
#' @param matrices list of hic matrices
#' @keywords internal
#' @return raises error or returns nothing
#' @noRd
#' @importFrom checkmate assertList assertCount assertDataFrame
#' @importFrom checkmate makeAssertCollection reportAssertions
.validHicMatrices <- function(matrices){
    matricesFails <- checkmate::makeAssertCollection()
    checkmate::assertList(
        x = matrices,
        types = "ContactMatrix",
        add = matricesFails)
    checkmate::assertCount(
        x = attr(matrices,"resolution"),
        na.ok = FALSE,
        add = matricesFails
    )
    checkmate::assertDataFrame(
        x = attr(matrices,"chromSize"),
        add = matricesFails
    )
    if(!matricesFails$isEmpty()){
        matricesFails$push("Unsupported HiC matrices list object")
        checkmate::reportAssertions(matricesFails)
    }
}

#' .validGranges
#' check the validity of GRanges objects
#' @keywords internal
#' @param gRanges object to test
#' @param testForList are GRanges in a list?
#' @param nullValid is null value permitted?
#'
#' @return raises error or returns nothing
#' @noRd
#' @importFrom checkmate assertClass assert checkClass checkList
.validGranges <- function(
    gRanges, testForList = FALSE,
    nullValid = FALSE){
    if(!testForList){
        checkmate::assertClass(
            x = gRanges,
            classes = "GRanges",
            null.ok = nullValid
        )
    }else{
        checkmate::assert(
            checkmate::checkClass(
                x = gRanges,
                classes = "GRanges",
                null.ok = nullValid
            ),
            checkmate::checkList(
                x = gRanges,
                types = "GRanges",
                null.ok = nullValid
            ),
            checkmate::checkClass(
                x = gRanges,
                classes = "GRangesList",
                null.ok = nullValid
            ), combine = "or"
        )
    }
}

#' .validSubmatrices
#' check the validity of submatrices
#' @keywords internal
#' @param submatrices list of submatrices
#'
#' @return raises error or returns nothing
#' @noRd
#' @importFrom checkmate checkClass assert checkList
.validSubmatrices <- function(submatrices){
    checkmate::assert(
        checkmate::checkList(
            x = submatrices,
            types = "matrix",
            null.ok = FALSE
        ),
        checkmate::checkClass(
            x = attr(submatrices, "interactions"),
            classes = "GInteractions",
            null.ok = FALSE
        )
    )
}
#' Add list as attributes.
#'
#' AddAttr
#' @keywords internal
#' @description Add list as attributes to any object with or without overwrite.
#' @param x <any>: An object to which attributes are to be added.
#' @param attrs <list>: A named list of new attributes.
#' @param overwrite <logical>: Whether an overwrite is required on attributes
#'  with the same name.(Default FALSE)
#' @return The object with new attributes.
#' @noRd
AddAttr <- function(
    x = NULL, attrs = NULL,
    overwrite = FALSE
) {
    intersectAttr <- intersect(
        names(attributes(x)),
        names(attrs)
    )
    if (overwrite && length(intersectAttr)) {
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

#' D.melanogaster Beaf-32 ChIP-seq.
#'
#' @description Drosophila Melanogaster Beaf32 peaks on 2L and 2R chromosomes.
#' @docType data
#' @usage data(Beaf32_Peaks.gnr)
#' @format An object of class GRanges.
#' @keywords datasets
#' @examples
#' data(Beaf32_Peaks.gnr)
#' Beaf32_Peaks.gnr
#'
"Beaf32_Peaks.gnr"

#' Data frames row binding.
#'
#' BindFillRows
#' @keywords internal
#' @description Bind data frames by rows after filling missing columns with NA.
#' @param df_Lst <data.frames or list[data.frame]>: Data frames to bind or list
#'  of data.frames. If is a data.frame create a list with arguments
#'  `df_Lst` and `...`, else `...` are ignored.
#' @param ... <data.frames or list[data.frame]>: Data frames to bind or list of
#'  data.frames.
#' @return The binded data frame
#' @noRd
BindFillRows <- function(
    df_Lst, ...
) {
    if (is.data.frame(df_Lst)) {
        df_Lst <- list(df_Lst, ...)
    }
    df_Lst <- lapply(
        seq_along(df_Lst),
        function(data.ndx) {
            data.df <- df_Lst[[data.ndx]]
            dataNames.chr <- lapply(df_Lst[-data.ndx], names) |>
                unlist() |>
                unique()
            data.df[setdiff(dataNames.chr, names(data.df))] <- NA
            return(data.df)
        }
    )
    return(do.call(rbind, df_Lst))
}

#' Blur a matrix.
#'
#' BoxBlur
#' @keywords internal
#' @description Blur a matrix with a one dimensional kernel.
#' @param mtx <matrix>: Numerical matrix.
#' @param boxKernel <numeric>: The numerical vector for kernel.
#'  If NULL apply a GaussBox (see 'GaussBox' function) (Default NULL)
#' @param kernSize <numeric>: If boxKernel is NULL, size of kernel for
#'  'GaussBox' function. (Default NULL)
#' @param stdev <numeric>: If boxKernel is NULL, standard deviation
#'  parameter for 'GaussBox' function. (Default NULL)
#' @return Blurred matrix.
#' @noRd
BoxBlur <- function(
    mtx, boxKernel = NULL, kernSize = NULL,
    stdev = 1
) {
    if (is.null(boxKernel)) {
        boxKernel <- GaussBox(
            stdev = stdev, kernScale = "1",
            kernSize = kernSize
        )
    }
    pad.num <- (length(boxKernel) - 1)/2
    mtx <- PadMtx(
        mtx = mtx, padSize = pad.num,
        val = NULL, side = c("top", "bot", "right", "left")
    )
    matVsmth.mtx2 <- vapply(
        ((1 + pad.num):(dim(mtx)[2] - pad.num)), FUN = 
        function(j) {
            (t(mtx[, (j - pad.num):(j + pad.num)]) *
                boxKernel) |>
                apply(2, Plus)
        }, FUN.VALUE = array(1,dim = dim(mtx)[2])
    )
    matHsmth.mtx2 <- t(vapply(
        ((1 + pad.num):(dim(matVsmth.mtx2)[1] - pad.num)), FUN = 
        function(i) {
            (matVsmth.mtx2[(i - pad.num):(i + pad.num), ] *
                boxKernel) |>
                apply(2, Plus) |>
                t()
        }, FUN.VALUE = array(1,dim = dim(matVsmth.mtx2)[2])
))
    indices <- which(matHsmth.mtx2 == 0)
    if (length(indices)) {
        matHsmth.mtx2 <- Rise0(matHsmth.mtx2, indices = indices)
    }
    return(matHsmth.mtx2)
}

#' Cut a vector.
#'
#' BreakVector
#' @keywords internal
#' @description Compute the n+1 breaks of a vector in a linear or density based
#'  way with the possibility to fix minimal, center and maximal values.
#' @param x <numeric>: Numerical vector.
#' @param x_min <numeric>: Minimal fixed value.
#' @param center <numeric>: Center fixed value.
#' @param x_max <numeric>: Maximal fixed value.
#' @param n <numeric>: Number of tile (return n+1 breaks).
#' @param method <character>: Kind of breaking. "linear" or "density".
#'  (Default "linear")
#' @return Numerical vector of breaks.
#' @noRd
BreakVector <- function(
    x = NULL, x_min = NULL,
    center = NULL, x_max = NULL,
    n = 10, method = "linear"
) {
    n <- n + 1
    if (method == "linear") {
        x <- x[which(!is.na(x))]
        if (is.null(x_min)) {
            x_min <- min(x, na.rm = TRUE)
        }
        if (is.null(x_max)) {
            x_max <- max(x, na.rm = TRUE)
        }
        if (is.null(center)) {
            breaks.num <- seq(
                x_min, x_max,
                length.out = n
            )
        } else if (x_min < center && center < x_max) {
            breaks.num <- c(
                seq(
                    x_min, center,
                    length.out = n%/%2 + 1
                ),
                seq(
                    center, x_max,
                    length.out = n%/%2 + 1
                )
            )
        } else {
            center <- stats::median(x, na.rm = TRUE)
            breaks.num <- c(
                seq(
                    x_min, center,
                    length.out = n%/%2 + 1
                ),
                seq(
                    center, x_max,
                    length.out = n%/%2 + 1
                )
            )
        }
    } else if (method == "density") {
        breaks.num <- stats::quantile(
            x, na.rm = TRUE,
            probs = seq(0, 1, length.out = n)
        )
    } else {
        stop("Method.chr muste be one of 'linear' or 'density'.\n")
    }
    return(breaks.num[!duplicated(breaks.num)])
}

#' Gaussian formula.
#'
#' Gauss
#' @keywords internal
#' @description Gaussian formula in 1 or 2 dimension.
#' @param x <numeric>: x value.
#' @param y <numeric>: y value for 2 dimensional gaussian.
#' @param stdev <numeric>: Standard deviation parameter of the gaussian.
#'  (Default 1)
#' @param mu <numeric>: Mean deviation parameter of the gaussian. (Default 0)
#' @return Result of Gaussian formula
#' @noRd
Gauss <- function(
    x = NULL, y = NULL, stdev = 1, mu = 0
) {
    x <- x[1]
    y <- y[1]
    if (is.null(y)) {
        return(1 / (stdev * sqrt(2*pi)) * exp(-((x - mu)^2) / (2*stdev^2)))
    } else {
        return(1 / (2 *pi*stdev^2) * exp(-((x^2 + y^2) / (2*stdev^2))))
    }
}

#' One dimension Gaussian kernel.
#'
#' GaussBox
#' @keywords internal
#' @description One dimension Gaussian kernel.
#' @param stdev <numeric>: Standard deviation parameter of the gaussian.
#'  (Default 1)
#' @param kernSize <numeric>: Kernel size. If NULL size is 1+4*stdev.
#'  (Default NULL)
#' @param kernScale <character>: Scaling kind of box. If "1" sum of kernel
#'  equals 1. If "int" Minimal value of kernel is 1 and all entry are integer.
#'  If "none", kernel is not scale. (Default "1")
#' @return numerical vector.
#' @noRd
GaussBox <- function(
    stdev = 1, kernSize = NULL, kernScale = "1"
) {
    if (is.null(kernSize)) {
        kernSize <- 1 + 4 * stdev
    }
    x <- as.vector(
        scale(seq_len(kernSize),
        scale = FALSE,
        center = TRUE)
    )
    box <- lapply(x, function(x) {
        xInterval.num <- seq((x - 0.5), (x + 0.5), by = 0.01) |>
            lapply(function(xi) {
                Gauss(x = xi, stdev = stdev)
            }) |>
            unlist() |>
            mean()
        return(xInterval.num)
    }) |>
        unlist()
    if (kernScale == "1") {
        box <- {
            box / sum(abs(box))
        }
    } else if (kernScale == "int") {
        box <- ceiling(box / box[1])
    }
    return(box)
}

#' Base pairs convertion.
#'
#' GenomicSystem
#' @description Convert numbers of base into string with order of magnitude
#'  (Kbp, Mbp, Gbp) and vice versa.
#' @param x <character or numeric>: The number
#' to convert or string to convert.
#' @param digits <integer>: The number of significant
#' digits to be used. See [signif()] for more informations.
#' (Default 3)
#' @return The converted number or string.
#' @export
#' @examples
#' GenomicSystem(1540, 3)
#' GenomicSystem(1540, 2)
#' GenomicSystem("1Mbp")
#' GenomicSystem("1Kbp")
#' GenomicSystem("1k")
#'
GenomicSystem <- function(x, digits = 3) {
    if (is.numeric(x)) {
        dplyr::case_when(
            x >= 1e+09 ~ paste0(signif(x * 10^(-9), digits),"Gbp"),
            x >= 1e+06 ~ paste0(signif(x * 10^(-6), digits),"Mbp"),
            x >= 1000 ~ paste0(signif(x * 10^(-3), digits),"Kbp"),
            x >= 0 ~ paste0(signif(x * 10^(0), digits), "Bp"))
    } else if (is.character(x)) {
        x <- toupper(x)
        tryCatch(as.numeric(x),
            error = function(e){
                stop("An error occurred on GenomicSystem conversion:\n", e)
            },
            warning = function(w){
                dplyr::case_when(
                grepl(x = x, pattern = "G") ~ (10^9) * 
                    as.numeric(gsub(
                    x = x,
                    pattern = "[a-z]",
                    ignore.case = TRUE,
                    replacement = ""
                    )),
                grepl(x = x, pattern = "M") ~ (10^6) * 
                    as.numeric(gsub(
                    x = x,
                    pattern = "[a-z]",
                    ignore.case = TRUE,
                    replacement = ""
                    )),
                grepl(x = x, pattern = "K") ~ (10^3) * 
                    as.numeric(gsub(
                    x = x,
                    pattern = "[a-z]",
                    ignore.case = TRUE,
                    replacement = ""
                    ))
                )
            }
        )
    }
}

#' Get file extension
#'
#' GetFileExtension
#' @keywords internal
#' @description Give the extension of a file from the path.
#' @param file <character>: The path to the file.
#' @return A character string
#' @noRd

GetFileExtension <- function(
    file = NULL
) {
    fileName.chr <- GetFileName(file, extension = TRUE) |>
        strsplit(".", fixed = TRUE) |>
        unlist()
    return(fileName.chr[length(fileName.chr)])
}

#' Get file name.
#'
#' GetFileName
#' @keywords internal
#' @description Function as `basename()` with the option to not return the file
#'  extension.
#' @param file <character>: The path to the file.
#' @param extension <logical>: Whether the file extension should be returned
#'  with the file name. (Default FALSE)
#' @return A character string.
#' @noRd
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

#' In situ Hi-C control.
#'
#' @description In situ Hi-C on non-heat treated S2 cells
#'  (Drosophila Melanogaster) with MboI on chromosome 2R and 2L download from
#'  [4DN](https://data.4dnucleome.org/experiment-set-replicates/4DNESFOADERB/)
#'  portal (Ray J, Munn PR, et al., 2019). This data is the result of the
#'  [HicAggR::ImportHiC()] function.
#' @docType data
#' @usage data(HiC_Ctrl.cmx_lst)
#' @format A a list of ContactMatrix objects. Each element correspond to the
#'  interaction matrix of two chromosomes.
#' @keywords datasets
#' @examples
#' data(HiC_Ctrl.cmx_lst)
#' HiC_Ctrl.cmx_lst
#'
"HiC_Ctrl.cmx_lst"

#' In situ Hi-C heat treated.
#'
#' @description In situ Hi-C on heat treated S2 cells (Drosophila Melanogaster)
#'  with MboI on chromosome 2R and 2L download from
#'  [4DN](https://data.4dnucleome.org/experiment-set-replicates/4DNESFI64TG3/)
#'  portal (Ray J, Munn PR, et al., 2019). This data is the result of the
#'  [HicAggR::ImportHiC()] function.
#' @docType data
#' @usage data(HiC_HS.cmx_lst)
#' @format A a list of ContactMatrix objects. Each element correspond to the
#'  interaction matrix of two chromosomes.
#' @keywords datasets
#' @examples
#' data(HiC_HS.cmx_lst)
#' HiC_HS.cmx_lst
#'
"HiC_HS.cmx_lst"

#' Convert HSL to Hex.
#'
#' Hsl2Hex
#' @keywords internal
#' @description Convert a color in HSL (Hue,Saturation,Light) to
#'  hexadecimal format.
#' @param hslColor <charcater>: A vector of the color's HSL code.
#' @param alpha <logical>: Whether the alpha layer should be returned.
#'  (Default FALSE)
#' @return A character of the color's hexadecimal code.
#' @noRd
#'
Hsl2Hex <- function(
    hslColor = NULL, alpha = FALSE
) {
    color.hex <- Hsl2Rgb(hslColor, alpha) |>
        Rgb2Hex(alpha)
    return(color.hex)
}

#' Convert HSL to RGB.
#'
#' Hsl2Rgb
#' @keywords internal
#' @description Convert a color in HSl (Hue,Saturation,Light) format to RGB.
#' @param hslColor <charcater>: A vector of the color's HSL code.
#' @param alpha <logical>: Whether the alpha layer should be returned.
#'  (Default FALSE)
#' @return An integer vector of the color's RGB code.
#' @noRd
#'
Hsl2Rgb <- function(
    hslColor = NULL, alpha = FALSE
) {
    if (3 > length(hslColor) || length(hslColor) > 4) {
        err.chr <- paste0(
            "Need 3 or 4 values beetween 0 and 255, first value for hue, ",
            "second for saturation, third for light and last for alpha"
        )
        stop(err.chr)
    } else {
        if (IsHsl(hslColor)) {
            if (length(hslColor) == 3) {
                alphaValue <- 255
            } else {
                alphaValue <- hslColor[4] * 255
                hslColor <- hslColor[seq_len(3)]
            }
            C <- hslColor[2] * (1 - abs(2 * hslColor[3] - 1))
            X <- C * (1 - abs((hslColor[1] / 60) %% 2 - 1))
            m <- hslColor[3] - C / 2
            if (0 <= hslColor[1] && hslColor[1] < 60) {
                rbgColor <- c(
                    red = (C + m) * 255,
                    green = (X + m) * 255,
                    blue = (0 + m) * 255
                )
            } else if (60 <= hslColor[1] && hslColor[1] < 120) {
                rbgColor <- c(
                    red = (X + m) * 255,
                    green = (C + m) * 255,
                    blue = (0 + m) * 255
                )
            } else if (120 <= hslColor[1] && hslColor[1] < 180) {
                rbgColor <- c(
                    red = (0 + m) * 255,
                    green = (C + m) * 255,
                    blue = (X + m) * 255
                )
            } else if (180 <= hslColor[1] && hslColor[1] < 240) {
                rbgColor <- c(
                    red = (0 + m) * 255,
                    green = (X + m) * 255,
                    blue = (C + m) * 255
                )
            } else if (240 <= hslColor[1] && hslColor[1] < 300) {
                rbgColor <- c(
                    red = (X + m) * 255,
                    green = (0 + m) * 255,
                    blue = (C + m) * 255
                )
            } else if (300 <= hslColor[1] && hslColor[1] < 360) {
                rbgColor <- c(
                    red = (C + m) * 255,
                    green = (0 + m) * 255,
                    blue = (X + m) * 255
                )
            }
            if (alpha) {
                rbgColor <- c(rbgColor, alpha = alphaValue)
            }
        } else {
            err.chr <- paste0(
                "Need 3 or 4 values beetween 0 and 255, first value for hue, ",
                "second for saturation, third for light and last for alpha"
            )
            stop(err.chr)
        }
    }
    return(round(rbgColor))
}

#' Hue palette.
#'
#' Hue
#' @description Create an Hue palette.
#' @param paletteLength <numeric>: Color number.
#' @param rotation <numeric>: If positive, rotates clockwise
#' in the color space, reversing if the number is negative. If is NULL,
#' compute rotation according to hueRange parameter. (Default NULL)
#' @param hueRange <numeric>: Degree range in color space
#' between 0 and 360. (Default c(0,360))
#' @param saturation <numeric>: Saturation value between 0 and 1.
#' (Default 0.65)
#' @param lightness <numeric>: Lightness value between 0 and 1.
#' (Default 0.65)
#' @param alphaValue <numeric>: Opacity value between 0 and 1.
#' (Default 1)
#' @param alpha <logical>: Whether the alpha layer should be
#' returned. (Default FALSE)
#' @return A vector of color.
#' @export
#' @examples
#' Hue(paletteLength = 9)
#'
Hue <- function(
    paletteLength = 9, rotation = NULL, hueRange = c(0, 360),
    saturation = 0.65, lightness = 0.65, alphaValue = 1,
    alpha = FALSE
) {
    if (paletteLength > 2) {
        if (is.null(rotation)) {
            if (abs(diff(hueRange)) >=
                180) {
                rotation <- sign(diff(hueRange))
            } else {
                rotation <- -sign(diff(hueRange))
            }
        }
        if (sign(diff(hueRange)) !=
            sign(rotation)) {
            hueRange[which.min(hueRange)] <- 360 +
                hueRange[which.min(hueRange)]
                
        }
        distGap.num <- seq(
            hueRange[1],
            hueRange[2],
            length.out = paletteLength) |>
            diff() |>
            mean() |>
            abs()
        gap.num <- min(
            abs(diff(hueRange)%%360),
            abs(360 - abs(diff(hueRange)%%360))
        )
        if (paletteLength > 1 && gap.num < distGap.num) {
            adjust.num <- (distGap.num - gap.num)/(paletteLength +1) *
                paletteLength
            hueRange[which.max(hueRange)] <- 
                hueRange[which.max(hueRange)] - adjust.num
        }
        hue.lst <- seq(
            hueRange[1],
            hueRange[2],
            length.out = paletteLength)%%360
    } else if (paletteLength == 2) {
        hue.lst <- Hue(
            paletteLength = 5, rotation = rotation,
            hueRange = hueRange, saturation = saturation,
            lightness = lightness, alphaValue = alphaValue,
            alpha = alpha)[c(2, 3)]
    } else if (paletteLength == 1) {
        hue.lst <- Hue(
            paletteLength = 3, rotation = rotation,
            hueRange = hueRange, saturation = saturation,
            lightness = lightness, alphaValue = alphaValue,
            alpha = alpha)[2]
    }
    if (is.numeric(hue.lst)) {
        hue.lst <- lapply(
            hue.lst, function(hue.num) {
                Hsl2Hex(
                    c(
                        hue = hue.num, saturation = saturation,
                        light = lightness, alpha = alphaValue
                    ),
                    alpha = alpha
                )
            }
        ) |>
            unlist()
        return(hue.lst)
    } else {
        return(hue.lst)
    }
}

#' Check HSL color format.
#'
#' IsHsl
#' @keywords internal
#' @description Check if a color is in HSL color format.
#' @param colour <character or numeric>: A color.
#' @return A logical.
#' @noRd
IsHsl <- function(
    colour = NULL
) {
    logical.bln <- lapply(
            colour[-1],
            function(val) {0 <= val & val <= 1 }
        ) |>
        unlist() |>
        sort()
    logical.bln <- logical.bln[[1]]
    return(
        class(colour) %in% c("list", "numeric", "integer") &&
        0 <= colour[[1]] &&
        colour[[1]] < 360 &&
        logical.bln
    )
}

#' Check RGB color format.
#'
#' IsRgb
#' @keywords internal
#' @description Check if a color is in RGB color format.
#' @param colour <character or numeric>: A color.
#' @return A logical.
#' @noRd
#'
IsRgb <- function(
    colour = NULL
) {
    logical.bln <- lapply(
            colour,
            function(val) {0 <= val & val <= 255}
        ) |>
        unlist() |>
        sort()
    logical.bln <- logical.bln[[1]]
    return(
        (!IsHsl(colour)) &&
        (class(colour) %in% c("list", "numeric", "integer") &&
        logical.bln)
    )
}

#' Configure parallel parameters.
#'
#' MakeParallelParam
#' @keywords internal
#' @description Create BiocParallel parameter according to OS.
#' @param cores <numerical> : An integer to specify the number of cores.
#'  (Default 1)
#' @param verbose <logical>: A logical value. If TRUE show the progression in
#'  console. (Default TRUE)
#' @return Parrallel parameter according number of cores and OS to use with
#'  BiocParallel package.
#' @noRd
#'
MakeParallelParam <- function(
    cores = 1, verbose = FALSE
) {
    if (!is.numeric(cores) ||
        cores < 2 ||
        .Platform$OS.type == "windows"
    ) {
        return(BiocParallel::SerialParam(progressbar = verbose))
    } else {
        return(
            BiocParallel::MulticoreParam(
                workers = cores,
                progressbar = verbose
            )
        )
    }
}

#' Scale values by mean.
#'
#' MeanScale
#' @keywords internal
#' @description Scale values with mean.
#' @param x <numeric>: Numerical vector.
#' @return Scaled numeric vector.
#' @noRd
MeanScale <- function(
    x
) {
    (x - mean(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x,na.rm = TRUE))
}

#' Coerce matrix in tibble.
#'
#' MeltSpm
#' @keywords internal
#' @description Coerce a sparse matrix M in tibble where columns:
#'  i is row index, j is column index and x the value M`[`i,j`]`.
#' @param spMtx <dgCMatrix or dgCMatrix coercible>: A matrix.
#' @return A tibble.
#' @noRd
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

#' Scales values on min-max range.
#'
#' MinMaxScale
#' @keywords internal
#' @description Scale values on min-max range.
#' @param x <numeric>: Numerical vector.
#' @param x_min <numeric>: Minimal value after scaling.
#' @param x_max <numeric>: Maximal value after scaling.
#' @return Scaled numeric vector.
#' @noRd
MinMaxScale <- function(
    x, x_min = (0), x_max = 1
) {
    x_min +
    (
        ((x - min(x, na.rm = TRUE)) * (x_max - x_min)) /
        (max(x, na.rm = TRUE ) - min(x, na.rm = TRUE))
    )
}

#' Exclusion binary operator.
#'
#' NotIn
#' @keywords internal
#' @description Binary operator, inverse to \%in\%.
#' @param lhs <vector or NULL>: Values to be compared against rhs
#' @param rhs <vector or NULL>: Values to be compared against lhs
#' @return A boolean.
#' @noRd
"NotIn" <- function(
    lhs, rhs
) {
    return(match(lhs, rhs, nomatch = 0L) == 0L)
}

#' Add a value around a matrix.
#'
#' PadMtx
#' @keywords internal
#' @description Add a value around a matrix.
#' @param mtx <matrix>: Numerical matrix.
#' @param padSize <numeric>: Number of columns or rows to add. (Default 1)
#' @param val <numeric>: Value to add. If Null create mirror of choosen sides.
#'  (Default 0)
#' @param side <character>: Side to pad, must be one or some of 'top','bot',
#' 'right' or 'left'. (Default c('top','bot','right','left') )
#' @return A matrix.
#' @noRd
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

#' Sum by removing NA.
#'
#' Plus
#' @keywords internal
#' @description Perform sum by removing NA. If all values are NA return NA
#'  instead 0.
#' @param x <numerical>: A numerical vector
#' @return  The sum of numbers.
#' @noRd
Plus <- function(
    x
) {
    if (all(is.na(x))) {
        c(x[0], NA)
    } else {
        sum(x, na.rm = TRUE)
    }
}

#' Find threshold for outliers based on quantiles.
#'
#' QtlThreshold
#' @keywords internal
#' @description Find threshold for outliers triming based on quantiles.
#' @param x <numeric>: Numeric vector.
#' @param prctThr <numeric>: Percentage (0-100) threshold.
#' (Default 5)
#' @param tails <character>: Bounds to return, "lower", "upper"
#' or "both". (Default "both")
#' @return Numerical vector of thresholds values for outliers triming.
#' @importFrom checkmate assertChoice assertNumeric
#' @export
#' @examples
#' set.seed(1111)
#' x <- 0:100
#' x <- sort(x)
#' x
#' QtlThreshold(x, prctThr = 5, tails = "lower")
#' QtlThreshold(x, prctThr = 5, tails = "both")
#' QtlThreshold(x, prctThr = 5, tails = "upper")
QtlThreshold <- function(
    x = NULL, prctThr = 5, tails = "both"
) {
    checkmate::assertChoice(
        x = tails,
        choices = c('lower', 'upper', 'both'),
        null.ok = FALSE
    )
    checkmate::assertNumeric(
        x = prctThr,
        lower = 0, upper = 100,
        null.ok = FALSE
    )
    checkmate::assertNumeric(
        x = x,
        null.ok = FALSE
    )
    probs.num <- dplyr::case_when(
        tails == "both" ~ c(prctThr / 200, 1 - (prctThr / 200)),
        tails == "upper" ~ c(NA, 1 - (prctThr / 100)),
        tails == "lower" ~ c(prctThr / 100, NA)
    )
    return(stats::quantile(x, na.rm = TRUE, probs.num))
}

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
#' @importFrom S4Vectors Rle runValue
ReduceRun <- function(
    firstRle, secondRle, reduceMethod = "paste", ...
) {
    if (methods::is(firstRle, "rle")) {
        firstRle <- S4Vectors::Rle(
        values = firstRle$values,
        lengths = firstRle$lengths)
    }
    if(is.numeric(S4Vectors::runValue(firstRle))) {
        firstValues <- as.numeric(firstRle)
    } else if (is.character(S4Vectors::runValue(firstRle))) {
        firstValues <- as.character(firstRle)
    } else if (is.factor(S4Vectors::runValue(firstRle))) {
        if(is.numeric(levels(S4Vectors::runValue(firstRle)))) {
            firstValues <- as.numeric(firstRle)
        } else if (is.character(levels(S4Vectors::runValue(firstRle)))) {
            firstValues <- as.character(firstRle)
        }
    }
    if (methods::is(secondRle, "rle")) {
        secondRle <- S4Vectors::Rle(
        values = secondRle$values,
        lengths = secondRle$lengths)
    }
    if (is.numeric(S4Vectors::runValue(secondRle))) {
        secondValues <- as.numeric(secondRle)
    }else if (is.character(S4Vectors::runValue(secondRle))) {
        secondValues <- as.character(secondRle)
    } else if (is.factor(S4Vectors::runValue(secondRle))) {
        if(is.numeric(levels(S4Vectors::runValue(secondRle)))) {
            secondValues <- as.numeric(secondRle)
        } else if (is.character(levels(S4Vectors::runValue(secondRle)))) {
            secondValues <- as.character(secondRle)
        }
    }
    vals.lst <- list(
        firstVal.vec = firstValues,
        secondVal.vec = secondValues
    )

    # These changes were made to simplify the code and remove super-assignment
    # Trials with list2env were consistent failures...

    if(is.character(firstValues)){
        # walking across the first values vector is risky,
        # but in both uses of this function (in ExtractSubmatrix),
        # first and second values should have the same length
        newVal.vec <- vapply(seq_along(firstValues), 
            FUN = function(x){eval(parse(text = reduceMethod))(
                vals.lst[[1]][x],
                vals.lst[[2]][x],
        ...
        )}, FUN.VALUE = character(1))
    } else if(is.numeric(firstValues)){
        newVal.vec <- vapply(seq_along(firstValues), 
            FUN = function(x){eval(parse(text = reduceMethod))(
                vals.lst[[1]][x],
                vals.lst[[2]][x],
        ...
        )}, FUN.VALUE = numeric(1))
    }
    return(S4Vectors::Rle(newVal.vec))
}

#' Resize a matrix
#'
#' ResizeMatrix
#' @keywords internal
#' @description Resize a numericam matrix in new dimension.
#' @param mtx <matrix>: A numerical matrix to resize.
#' @param newDim <integer>: The number of rows and cols in
#' resized matrix.
#' @return Resized matrix.
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


#' Convert RGB to Hex.
#'
#' Rgb2Hex
#' @keywords internal
#' @description Convert a color in RGB format to hexadecimal format.
#' @param rbgColor <integer>: An integer of the color's RGB code.
#' @param alpha <logical>: Whether the alpha layer should be returned.
#'  (Default FALSE)
#' @return A character of the color's hexadecimal code.
#' @noRd
Rgb2Hex <- function(
    rbgColor = NULL, alpha = FALSE
) {
    if (3 > length(rbgColor) || length(rbgColor) > 4) {
        err.chr <- paste0(
            "Need 3 or 4 values beetween 0 and 255, first value for red, ",
            "second for green, third for blue and last for alpha"
        )
        stop(err.chr)
    } else {
        if (IsRgb(rbgColor)) {
            if (length(rbgColor) == 3) {
                rbgColor <- c(rbgColor, 255)
            }
            hex.col <- lapply(rbgColor, function(val) {
                fisrtBit <- val %/% 16
                if (fisrtBit > 9) {
                    fisrtBit <- letters[fisrtBit - 9]
                }
                secondBit <- val %% 16
                if (secondBit > 9) {
                    secondBit <- letters[secondBit - 9]
                }
                return(c(fisrtBit, secondBit))
            }) |>
                unlist()
            hex.col <- paste0(c("#", hex.col), collapse = "")
        } else {
            err.chr <- paste0(
                "Need 3 or 4 values beetween 0 and 255, first value for red, ",
                "second for green, third for blue and last for alpha"
            )
            stop(err.chr)
        }
    }
    if (!alpha) {
        hex.col <- substr(hex.col, 1, nchar(hex.col) - 2)
    }
    return(hex.col)
}

#' Explicit zeros in sparse matrix.
#'
#' Rise0
#' @keywords internal
#' @description Explicit some implicit zeros in sparse matrix.
#' @param spMtx <dgCMatrix or dgCMatrix coercible>: A sparse matrix.
#' @param indices <numeric>: Vector of positions of the zeros to be explicits
#'  (column driven). If NULL and coords NULL all zeros are explicits.
#'  (Default NULL)
#' @param coords <data.frame>: A coordinate data frame for zeros to explicit
#'  Row index in fisrt column, columns index in second columns. If NULL the
#'  indices parameter is used (Default NULL)
#' @return Sparse matrix with some explicit zeros.
#' @noRd
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

#' Find threshold for outliers based on sd.
#'
#' SdThreshold
#' @keywords internal
#' @description Find threshold to trim outliers based on standard deviation.
#' @param x <numeric>: numeric vector.
#' @param sdThr <numeric>: number of standard deviation.
#' (Default 3)
#' @param tails <character>: bounds to return, "lower", "upper"
#' or "both". (Default "both")
#' @return numerical vector of thresholds values for outliers triming
#' @noRd 
SdThreshold <- function(
    x = NULL, sdThr = 3, tails = "both"
) {
    thr <- dplyr::case_when(
        tails == "both" ~ c(
            (mean(x, na.rm = TRUE) - (sdThr * stats::sd(x, na.rm = TRUE))),
            (mean(x, na.rm = TRUE) + (sdThr * stats::sd(x, na.rm = TRUE)))),
        tails == "upper" ~ c(NA, 
            (mean(x, na.rm = TRUE) + (sdThr * stats::sd(x, na.rm = TRUE)))),
        tails == "lower" ~ c((
            mean(x, na.rm = TRUE) - (sdThr * stats::sd(x, na.rm = TRUE))), NA)
    )
    return(thr)
}

#' Get all sequences lengths.
#'
#' SeqEnds
#' @description Get all sequences lengths for each ranges of a GRanges object.
#' @param gRanges <GRanges>: A GRanges object.
#' @return An integer vector.
#' @export
#' @examples
#' GRange.grn <- GenomicRanges::GRanges(
#'     seqnames = S4Vectors::Rle(c("chr1", "chr2", "chr1"), c(1, 3, 1)),
#'     ranges = IRanges::IRanges(101:105, end = 111:115,
#'         names = letters[seq_len(5)]),
#'     strand = S4Vectors::Rle(BiocGenerics::strand(c("-", "+", "*", "+")),
#'         c(1, 1, 2, 1)),
#'     seqinfo = c(chr1 = 200, chr2 = 300),
#'     score = seq_len(5)
#' )
#' SeqEnds(GRange.grn)
#'
SeqEnds <- function(
    gRanges
) {
    GenomeInfoDb::seqlengths(gRanges)[
        as.character(GenomeInfoDb::seqnames(gRanges))
    ]
}

#' Convert String to GRanges.
#'
#' StrToGRanges
#' @description Convert ranges describe with string 
#' (i.e seqname:start-end:strand) in GRanges object.
#' @param stringRanges <character>: Strings to convert on
#' GRanges.
#' @return A GRanges object.
#' @export
#' @importFrom S4Vectors mcols
#' @examples
#' StrToGRanges("chr1:1-100:+")
#' StrToGRanges(c("chr1:1-100:+", "chr2:400-500:-", "chr1:10-50:*"))
#'
StrToGRanges <- function(
    stringRanges
) {
    x.gnr <- lapply(stringRanges, function(x.chr) {
        x.chr <- unlist(strsplit(x.chr, ":"))
        seqnames.chr <- x.chr[1]
        ranges.num <- strsplit(x.chr[2], "-") |>
            unlist() |>
            as.numeric()
        start.num <- ranges.num[1]
        end.num <- ifelse(
            is.na(ranges.num[2]),
            yes = start.num,
            no = ranges.num[2]
        )
        strand.chr <- ifelse(
            is.na(x.chr[3]),
            yes = "*",
            no = x.chr[3]
        )
        x.gnr <- GenomicRanges::GRanges(
            seqnames = seqnames.chr,
            ranges = IRanges::IRanges(
                start = start.num,
                end = end.num
            ),
            strand = strand.chr
        )
        return(x.gnr)
    }) |>
        MergeGRanges()
    S4Vectors::mcols(x.gnr)$names <- names(stringRanges)
    return(x.gnr)
}

#' D.melanogaster TADs.
#'
#' @description Drosophila Melanogaster TADs on chromosome 2R and 2L
#'  ([F.Ramirez, 2018](https://doi.org/10.1038/s41467-017-02525-w))
#' @docType data
#' @usage data(TADs_Domains.gnr)
#' @format An object of class GRanges.
#' @keywords datasets
#' @examples
#' data(TADs_Domains.gnr)
#' TADs_Domains.gnr
#'
"TADs_Domains.gnr"

#' D.melanogaster Transcription starting sites.
#'
#' @description Data from a CHip Seq experiment
#' @docType data
#' @usage data(TSS_Peaks.gnr)
#' @format An object of class GRanges.
#' @keywords datasets
#' @examples
#' data(TSS_Peaks.gnr)
#' TSS_Peaks.gnr
#'
"TSS_Peaks.gnr"

#' Turns a nested list "inside-out".
#'
#' TransposeList
#' @keywords internal
#' @description Turns a nested list "inside-out".
#' @param nestedList <list[list]>: A nested list to transpose.
#' @return The tranposed nested list.
#' @noRd
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

#' Remove outliers.
#'
#' TrimOutliers
#' @keywords internal
#' @description Replace values of a numerical vector that are below a minimal
#'  thresholds and/or above maximal thresholds.
#' @param x <numeric>: Numeric vector.
#' @param thr <numeric>: Numeric vector of length 2. first
#' value is minimal threshold, second value maximal threshold
#' (Default find threshold based on standarrd deviation.
#' see `SdThreshold` function)
#' @param clip <logical>: If TRUE the values out of bounds are
#' replaced with thresholds values. If FALSE the Values out of bound are
#' replaced with NA (Default FALSE).
#' @return Trimed Numerical vector.
TrimOutliers <- function(
    x,
    thr = SdThreshold(x),
    clip = FALSE
) {
    if (clip) {
        x[which(x > thr[2])] <- thr[2]
        x[which(x < thr[1])] <- thr[1]
    } else {
        x[which(x > thr[2] | thr[1] > x)] <- NA
    }
    return(x)
}


#' viridis palette.
#'
#' viridis
#' @description Create a viridis palette.
#' @param paletteLength <numeric>: Color number.
#' @param space <numeric>: A character string; interpolation
#' in RGB or CIE Lab color spaces. See ?grDevices::colorRamp for more details.
#'  (Default "rgb")
#' @param interpolationMethod <numeric>: Use spline or linear
#' interpolation. See ?grDevices::colorRamp for more details.
#' (Default "linear")
#' @param bias <numeric>: A positive number. Higher values give
#' more widely spaced colors at the high end. See ?grDevices::colorRamp
#' for more details. (Default 1)
#' @return A vector of color.
#' @export
#' @examples
#' viridis(9)
#'
viridis <- function(
    paletteLength = NULL, space = "rgb",
    interpolationMethod = "linear", bias = 1
) {
    (grDevices::colorRampPalette(
        colors = c(
            "#440154", "#440256", "#450457", "#450559", "#46075A", "#46085C",
            "#460A5D", "#460B5E", "#470D60", "#470E61", "#471063", "#471164",
            "#471365", "#481467", "#481668", "#481769", "#48186A", "#481A6C",
            "#481B6D", "#481C6E", "#481D6F", "#481F70", "#482071", "#482173",
            "#482374", "#482475", "#482576", "#482677", "#482878", "#482979",
            "#472A7A", "#472C7A", "#472D7B", "#472E7C", "#472F7D", "#46307E",
            "#46327E", "#46337F", "#463480", "#453581", "#453781", "#453882",
            "#443983", "#443A83", "#443B84", "#433D84", "#433E85", "#423F85",
            "#424086", "#424186", "#414287", "#414487", "#404588", "#404688",
            "#3F4788", "#3F4889", "#3E4989", "#3E4A89", "#3E4C8A", "#3D4D8A",
            "#3D4E8A", "#3C4F8A", "#3C508B", "#3B518B", "#3B528B", "#3A538B",
            "#3A548C", "#39558C", "#39568C", "#38588C", "#38598C", "#375A8C",
            "#375B8D", "#365C8D", "#365D8D", "#355E8D", "#355F8D", "#34608D",
            "#34618D", "#33628D", "#33638D", "#32648E", "#32658E", "#31668E",
            "#31678E", "#31688E", "#30698E", "#306A8E", "#2F6B8E", "#2F6C8E",
            "#2E6D8E", "#2E6E8E", "#2E6F8E", "#2D708E", "#2D718E", "#2C718E",
            "#2C728E", "#2C738E", "#2B748E", "#2B758E", "#2A768E", "#2A778E",
            "#2A788E", "#29798E", "#297A8E", "#297B8E", "#287C8E", "#287D8E",
            "#277E8E", "#277F8E", "#27808E", "#26818E", "#26828E", "#26828E",
            "#25838E", "#25848E", "#25858E", "#24868E", "#24878E", "#23888E",
            "#23898E", "#238A8D", "#228B8D", "#228C8D", "#228D8D", "#218E8D",
            "#218F8D", "#21908D", "#21918C", "#20928C", "#20928C", "#20938C",
            "#1F948C", "#1F958B", "#1F968B", "#1F978B", "#1F988B", "#1F998A",
            "#1F9A8A", "#1E9B8A", "#1E9C89", "#1E9D89", "#1F9E89", "#1F9F88",
            "#1FA088", "#1FA188", "#1FA187", "#1FA287", "#20A386", "#20A486",
            "#21A585", "#21A685", "#22A785", "#22A884", "#23A983", "#24AA83",
            "#25AB82", "#25AC82", "#26AD81", "#27AD81", "#28AE80", "#29AF7F",
            "#2AB07F", "#2CB17E", "#2DB27D", "#2EB37C", "#2FB47C", "#31B57B",
            "#32B67A", "#34B679", "#35B779", "#37B878", "#38B977", "#3ABA76",
            "#3BBB75", "#3DBC74", "#3FBC73", "#40BD72", "#42BE71", "#44BF70",
            "#46C06F", "#48C16E", "#4AC16D", "#4CC26C", "#4EC36B", "#50C46A",
            "#52C569", "#54C568", "#56C667", "#58C765", "#5AC864", "#5CC863",
            "#5EC962", "#60CA60", "#63CB5F", "#65CB5E", "#67CC5C", "#69CD5B",
            "#6CCD5A", "#6ECE58", "#70CF57", "#73D056", "#75D054", "#77D153",
            "#7AD151", "#7CD250", "#7FD34E", "#81D34D", "#84D44B", "#86D549",
            "#89D548", "#8BD646", "#8ED645", "#90D743", "#93D741", "#95D840",
            "#98D83E", "#9BD93C", "#9DD93B", "#A0DA39", "#A2DA37", "#A5DB36",
            "#A8DB34", "#AADC32", "#ADDC30", "#B0DD2F", "#B2DD2D", "#B5DE2B",
            "#B8DE29", "#BADE28", "#BDDF26", "#C0DF25", "#C2DF23", "#C5E021",
            "#C8E020", "#CAE11F", "#CDE11D", "#D0E11C", "#D2E21B", "#D5E21A",
            "#D8E219", "#DAE319", "#DDE318", "#DFE318", "#E2E418", "#E5E419",
            "#E7E419", "#EAE51A", "#ECE51B", "#EFE51C", "#F1E51D", "#F4E61E",
            "#F6E620", "#F8E621", "#FBE723", "#FDE725"
        ),
        space = space,
        interpolate = interpolationMethod,
        bias = bias
    ))(paletteLength)
}

#' Convert string to function.
#'
#' WrapFunction
#' @keywords internal
#' @description Wrap a string into a function.
#' @param ... <character>: A string that could be parse and
#' eval as a function.
#' @return The result of the function or a function.
WrapFunction <- function(
    ...
) {
    eval(parse(text = paste(..., collapse = " ")))
}

#' YlGnBu palette.
#'
#' YlGnBu
#' @description Create a YlGnBu palette.
#' @param paletteLength <numeric>: Color number.
#' @param space <numeric>: A character string; interpolation
#' in RGB or CIE Lab color spaces. See ?grDevices::colorRamp for more details.
#'  (Default "rgb")
#' @param interpolationMethod <numeric>: Use spline or linear
#' interpolation. See ?grDevices::colorRamp for more details.
#' (Default "linear")
#' @param bias <numeric>: A positive number. Higher values give
#' more widely spaced colors at the high end. See `?grDevices::colorRamp`
#' for more details.
#' (Default 1)
#' @return A vector of color.
#' @export
#' @examples
#' YlGnBu(9)
#'
YlGnBu <- function(
    paletteLength = NULL, space = "rgb",
    interpolationMethod = "linear", bias = 1
) {
    (grDevices::colorRampPalette(
        colors = c(
            "#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4",
            "#1D91C0", "#225EA8", "#253494", "#081D58"
        ),
        space = space,
        interpolate = interpolationMethod,
        bias = bias
    ))(paletteLength)
}

#' YlOrRd palette.
#'
#' YlOrRd
#' @description Create a YlOrRd palette.
#' @param paletteLength <numeric>: Color number.
#' @param space <numeric>: A character string; interpolation
#' in RGB or CIE Lab color spaces. See `?grDevices::colorRamp`
#' for more details. (Default "rgb")
#' @param interpolationMethod <numeric>: Use spline or linear
#' interpolation. See `?grDevices::colorRamp` for more details.
#' (Default "linear")
#' @param bias <numeric>: A positive number. Higher values give
#' more widely spaced colors at the high end. See `?grDevices::colorRamp`
#' for more details. (Default 1)
#' @return A vector of color.
#' @export
#' @examples
#' YlOrRd(9)
#'
YlOrRd <- function(
    paletteLength = NULL, space = "rgb",
    interpolationMethod = "linear", bias = 1
) {
    (grDevices::colorRampPalette(
        colors = c(
            "#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C",
            "#FC4E2A", "#E31A1C", "#BD0026", "#800026"
        ),
        space = space,
        interpolate = interpolationMethod,
        bias = bias
    ))(paletteLength)
}
