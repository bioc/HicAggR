#' Compute quantification on extracted submatrices.
#'
#' GetQuantif
#' @description Function that computes quantification of contact frequencies in a given area and returns it in a named vector.
#' @param matrices <List[matrix]>: A matrices list.
#' @param areaFun <character or function>: A character or function that allows to extract an area from each matrix that composes the matrices list (Default "center").
#' \itemize{
#' \item "C" or "CENTER": 3x3 pixels at the intersection between anchor and bait.
#' \item "UL" or "UPPER_LEFT": pixels in the uppper left square
#' \item "UR" or "UPPER_RIGHT": pixels in the uppper right square
#' \item "BL" or "BOTTOM_LEFT": pixels in the bottom left square
#' \item "BR" or "BOTTOM_RIGHT": pixels in the bottom right square
#' \item "U" or "UPPER": pixels above the center area
#' \item "B" or "BOTTOM": pixels below the center area
#' \item "L" or "LEFT": pixels in the left of the center area
#' \item "R" or "RIGHT": pixels in the right of the center area
#' \item "D" or "DONUT": pixels that surrounds the center area
#' }
#' @param operationFun <character or function>: A character or function specifying the operation applied to the selected area (Default "mean_rm0").
#' \itemize{
#' \item "mean_rm0": apply a mean after replacing 0 with NA
#' \item "median_rm0": apply a median after replacing 0 with NA
#' \item "sum_rm0": apply a sum after replacing 0 with NA
#' \item "median": apply a median
#' \item "sum": apply a sum
#' \item "mean" or other character: apply a mean
#' }
#' @param varName <character>: The name of a column in GInteraction attributes of matrices used as named in the output vector (Default NULL). By default, sub-matrices IDs are used.
#' @return A GRange object.
#' @examples
#' # Data
#' data(Beaf32_Peaks.gnr)
#' data(HiC_Ctrl.cmx_lst)
#'
#' # Index Beaf32
#' Beaf32_Index.gnr <- IndexFeatures(
#'     gRangeList = list(Beaf = Beaf32_Peaks.gnr),
#'     chromSizes = data.frame(
#'         seqnames = c("2L", "2R"),
#'         seqlengths = c(23513712, 25286936)
#'     ),
#'     binSize = 100000
#' )
#'
#' # Beaf32 <-> Beaf32 Pairing
#' Beaf_Beaf.gni <- SearchPairs(indexAnchor = Beaf32_Index.gnr)
#' Beaf_Beaf.gni <- Beaf_Beaf.gni[seq_len(2000)] # subset 2000 first for exemple
#'
#' # Matrices extractions center on Beaf32 <-> Beaf32 point interaction
#' interactions_PF.mtx_lst <- ExtractSubmatrix(
#'     genomicFeature = Beaf_Beaf.gni,
#'     hicLst = HiC_Ctrl.cmx_lst,
#'     referencePoint = "pf"
#' )
#' GetQuantif(
#'     matrices = interactions_PF.mtx_lst,
#'     areaFun = "center",
#'     operationFun = "mean"
#' ) |> head()
#'
GetQuantif <- function(
    matrices, areaFun = "center",
    operationFun = "mean_rm0", varName = NULL
) {
    # Define operationFunction
    if (is.null(operationFun)) {
        operationFun <- function(area) {
            c(area[which(!is.na(area))])
        }
    } else if (!is.function(operationFun)) {
        operationFun <- dplyr::case_when(
            operationFun == "mean_rm0" ~
                "function(x){x[which(x==0)]<-NA;mean(x,na.rm=TRUE)}",
            operationFun == "median_rm0" ~
                "function(x){x[which(x==0)]<-NA;stats::median(x,na.rm=TRUE)}",
            operationFun == "sum_rm0" ~
                "function(x){x[which(x==0)]<-NA;sum(x,na.rm=TRUE)}",
            operationFun == "median" ~
                "function(x){stats::median(x,na.rm=TRUE)}",
            operationFun == "sum" ~
                "function(x){sum(x,na.rm=TRUE)}",
            operationFun == "mean" ~
                "function(x){mean(x,na.rm=TRUE)}",
            TRUE ~
                "function(x){mean(x,na.rm=TRUE)}"
        )
        operationFun <- WrapFunction(operationFun)
    }
    # Define extraction function
    matriceDim <- attributes(matrices)$matriceDim
    if (!is.function(areaFun) &&
        attributes(matrices)$referencePoint == "pf"
    ) {
        # Compute rows and cols index Center
        center.num <- c(
            floor((matriceDim + 1)/2) +
                ifelse(matriceDim >= 9, yes = 1, no = 0),
            ceiling((matriceDim + 1)/2) -
                ifelse(matriceDim >= 9, yes = 1, no = 0)
        )
        centerStart.num <- min(center.num)
        centerEnd.num <- max(center.num)
        center.chr <- ifelse(
            centerStart.num == centerEnd.num,
            yes = centerStart.num,
            no = paste0(centerStart.num, ":", centerEnd.num)
        )
        # Rest
        first.chr  <- paste0("1:", centerStart.num - 2)
        second.chr <- paste0(centerEnd.num + 2, ":", matriceDim)
        # Compute rows and cols index for Donut Thick
        donutThick.num <- ifelse(matriceDim >= 9, yes = 2, no = 1)
        # Top
        topDonutStart.num <- (centerStart.num - 1 - donutThick.num)
        topDonutEnd.num   <- centerStart.num - 2
        topDonut.num      <- topDonutStart.num:topDonutEnd.num
        # Left
        leftDonut.num <- topDonut.num
        # Bot
        botDonutEnd.num   <- centerEnd.num + 1 + donutThick.num
        botDonutStart.num <- centerEnd.num + 2
        botDonut.num      <- botDonutStart.num:botDonutEnd.num
        # Right
        rightDonut.num <- botDonut.num
        # Width
        widthDonut.num <- (centerStart.num - 1):(centerEnd.num + 1)
        # Height
        heigthDonut.num <- widthDonut.num
        # Coords
        donutCoord.dtf <- rbind(
            expand.grid(topDonut.num,    leftDonut.num),
            expand.grid(topDonut.num,    widthDonut.num),
            expand.grid(topDonut.num,    rightDonut.num),
            expand.grid(heigthDonut.num, leftDonut.num),
            expand.grid(heigthDonut.num, rightDonut.num),
            expand.grid(botDonut.num,    leftDonut.num),
            expand.grid(botDonut.num,    widthDonut.num),
            expand.grid(botDonut.num,    rightDonut.num)
        )
        donutRows.chr <- paste(donutCoord.dtf[, 1], collapse = ",")
        donutCols.chr <- paste(donutCoord.dtf[, 2], collapse = ",")
        donut.chr <- paste0(
            "cbind(c(",
            donutRows.chr,
            "), c(",
            donutCols.chr, ") )"
        )
        # Wrap index in a function
        areaFun <- dplyr::case_when(
            toupper(areaFun) %in% c("C", "CENTER") ~
                list(center.chr, center.chr),
            toupper(areaFun) %in% c("UL", "UPPER_LEFT") ~
                list(first.chr, first.chr),
            toupper(areaFun) %in% c("UR", "UPPER_RIGHT") ~
                list(first.chr, second.chr),
            toupper(areaFun) %in% c("BL", "BOTTOM_LEFT") ~
                list(second.chr, first.chr),
            toupper(areaFun) %in% c("BR", "BOTTOM_RIGHT") ~
                list(second.chr, second.chr),
            toupper(areaFun) %in% c("U", "UPPER") ~
                list(first.chr, center.chr),
            toupper(areaFun) %in% c("B", "BOTTOM") ~
                list(second.chr, center.chr),
            toupper(areaFun) %in% c("L", "LEFT") ~
                list(center.chr, first.chr),
            toupper(areaFun) %in% c("R", "RIGHT") ~
                list(center.chr, second.chr),
            toupper(areaFun) %in% c("D", "DONUT") ~
                list(donut.chr),
            TRUE ~
                list(center.chr, center.chr)
        ) |>
            paste(collapse = ",")
        areaFun <- WrapFunction(
            paste0("function(matrice.mtx){ matrice.mtx[", areaFun, "] }")
        )
    } else if (!is.function(areaFun) &&
        attributes(matrices)$referencePoint == "rf"
    ) {
        shift <- attributes(matrices)$shift
        # Compute rows and cols index Anchor
        anchorStart.num <- max(1,
            floor((matriceDim-2)*shift/(1+2*shift))+
                1 - ifelse(matriceDim >= 9, yes = 1, no = 0)
        )
        anchorEnd.num <- max(1,
            min(matriceDim,
                floor((matriceDim-2)*shift/(1+2*shift))+
                1 + ifelse(matriceDim >= 9, yes = 1, no = 0)
            )
        )
        anchor.chr <- ifelse(
            anchorStart.num == anchorEnd.num,
            yes = anchorStart.num,
            no = paste0(anchorStart.num, ":", anchorEnd.num)
        )
        # Bait
        baitStart.num <- max(1,
            ceiling((matriceDim-2)*(1+shift)/(1+2*shift)) +
                2 - ifelse(matriceDim >= 9, yes = 1, no = 0)
        )
        baitEnd.num <- max(1,
            min(matriceDim,
                ceiling((matriceDim-2)*(1+shift)/(1+2*shift))+
                2 + ifelse(matriceDim >= 9, yes = 1, no = 0)
            )
        )
        bait.chr <- ifelse(
            baitStart.num == baitEnd.num,
            yes = baitStart.num,
            no = paste0(baitStart.num, ":", baitEnd.num)
        )
        # Rest
        first.chr <- paste0("1:", anchorStart.num - 2)
        second.chr <- paste0(baitEnd.num + 2, ":", matriceDim)
        ULwidth.chr <- paste0("1:", baitStart.num - 2)
        inner.chr <- paste0(anchorEnd.num + 2, ":", baitStart.num - 2)
        BRheight.chr <- paste0(anchorEnd.num + 2, ":", matriceDim)
        donut.chr <- NULL
        # Computability
        U.lgk <- anchorStart.num >= 3
        R.lgk <- matriceDim >= baitEnd.num + 2
        B.lgk <- (((matriceDim + 1)/2) - anchorEnd.num) >= 1
        L.lgk <- (baitStart.num - ((matriceDim + 1)/2)) >= 1
        UL.lgk <- U.lgk && baitStart.num >= 3
        UR.lgk <- U.lgk && R.lgk
        BR.lgk <- matriceDim >= anchorEnd.num + 2 &&
            R.lgk
        BL.lgk <- (((matriceDim + 1)/2) - anchorEnd.num) >= 2 &&
            (baitStart.num -((matriceDim + 1)/2)) >= 2
        D.lgk <- sum(
            U.lgk, R.lgk, B.lgk, L.lgk,
            UL.lgk, UR.lgk, BR.lgk, BL.lgk) >= 1
        # Compute rows and cols index for Donut
        if (D.lgk) {
            # Thick
            donutThick.num <- ifelse(matriceDim >= 9, yes = 2, no = 1)
            # Top
            topDonutStart.num <- (anchorStart.num - 1 - donutThick.num)
            topDonutEnd.num <- anchorStart.num - 2
            topDonut.num <- topDonutStart.num:topDonutEnd.num
            # Left
            leftDonutStart.num <- (baitStart.num - 1 - donutThick.num)
            leftDonutEnd.num <- baitStart.num - 2
            leftDonut.num <- leftDonutStart.num:leftDonutEnd.num
            # Bot
            botDonutEnd.num <- anchorEnd.num + 1 + donutThick.num
            botDonutStart.num <- anchorEnd.num + 2
            botDonut.num <- botDonutStart.num:botDonutEnd.num
            # Right
            rightDonutEnd.num <- baitEnd.num + 1 + donutThick.num
            rightDonutStart.num <- baitEnd.num + 2
            rightDonut.num <- rightDonutStart.num:rightDonutEnd.num
            # Widht
            widthDonut.num <- (baitStart.num - 1):(baitEnd.num + 1)
            # Height
            heigthDonut.num <- (anchorStart.num - 1):(anchorEnd.num + 1)
            # Coord
            donutCoord.dtf <- rbind(
                expand.grid(topDonut.num, leftDonut.num),
                expand.grid(topDonut.num, widthDonut.num),
                expand.grid(topDonut.num, rightDonut.num),
                expand.grid(heigthDonut.num, leftDonut.num),
                expand.grid(heigthDonut.num, rightDonut.num),
                expand.grid(botDonut.num, leftDonut.num),
                expand.grid(botDonut.num, widthDonut.num),
                expand.grid(botDonut.num, rightDonut.num)
            )
            donutCoord.dtf <- dplyr::filter(
                donutCoord.dtf,
                donutCoord.dtf$Var1 >= 1 &
                donutCoord.dtf$Var2 <= matriceDim &
                donutCoord.dtf$Var1 <= donutCoord.dtf$Var2
            )
            donutRows.chr <- paste(donutCoord.dtf[, 1], collapse = ",")
            donutCols.chr <- paste(donutCoord.dtf[, 2], collapse = ",")
            donut.chr <- paste0(
                "cbind( c(",
                donutRows.chr,
                "), c(",
                donutCols.chr,
                ") )"
            )
        }
        # Wrap index in a function
        areaFun <- dplyr::case_when(
            toupper(areaFun) %in% c("C", "CENTER") ~
                list(anchor.chr, bait.chr),
            toupper(areaFun) %in% c("UL", "UPPER_LEFT") && UL.lgk ~ 
                list(first.chr, ULwidth.chr),
            toupper(areaFun) %in% c("UR", "UPPER_RIGHT") && UR.lgk ~ 
                list(first.chr, second.chr),
            toupper(areaFun) %in% c("BL", "BOTTOM_LEFT") && BL.lgk ~ 
                list(inner.chr, inner.chr),
            toupper(areaFun) %in% c("BR", "BOTTOM_RIGHT") && BR.lgk ~ 
                list(BRheight.chr, second.chr),
            toupper(areaFun) %in% c("U", "UPPER") && U.lgk ~ 
                list(first.chr, bait.chr),
            toupper(areaFun) %in% c("B", "BOTTOM") && B.lgk ~ 
                list(second.chr, bait.chr),
            toupper(areaFun) %in% c("L", "LEFT") && L.lgk ~ 
                list(anchor.chr, first.chr),
            toupper(areaFun) %in% c("R", "RIGHT") && R.lgk ~ 
                list(anchor.chr, second.chr),
            toupper(areaFun) %in% c("D", "DONUT") && D.lgk ~ 
                list(donut.chr),
            TRUE ~ 
                list(anchor.chr, bait.chr)
        ) |>
            paste(collapse = ",")
        areaFun <- WrapFunction(
            paste0("function(matrice.mtx){ matrice.mtx[", areaFun, "] }")
        )
    }
    # Compute quantif
    quantif.num <- lapply(
        matrices, function(mtx) {
            mtxQuantif.num <- operationFun(areaFun(mtx)) |>
                stats::setNames(NULL)
            rownames(mtxQuantif.num) <- NULL
            colnames(mtxQuantif.num) <- NULL
            return(mtxQuantif.num)
        }
    )
    # Get Names
    if (!is.null(varName)) {
        interactions.dtf <- data.frame(
            S4Vectors::mcols(attributes(matrices)$interactions)
            )
        names.chr_lst <- dplyr::arrange(
            interactions.dtf,
            factor(interactions.dtf$submatrix.name,
            levels = names(quantif.num))
        ) |>
            dplyr::pull(varName)
    } else {
        names.chr_lst <- names(quantif.num)
    }
    # Repeted Index if names.chr_lst is a nested List
    lengths.num <- lapply(names.chr_lst, length)
    lengths.num <- lengths.num |>
        stats::setNames(seq_along(lengths.num))
    repeted.ndx <- rep(
        names(lengths.num),
        lengths.num
    ) |> as.numeric()
    # Add attributes
    quantif.num <- unlist(quantif.num[repeted.ndx]) |>
        stats::setNames(unlist(names.chr_lst)) |>
        AddAttr(
            attrs = c(
                attributes(matrices),
                interactions =
                    attributes(matrices)$interactions[repeted.ndx],
                operation = operationFun,
                area = areaFun,
                duplicated = which(duplicated(repeted.ndx))
            )
        )
    return(quantif.num)
}
