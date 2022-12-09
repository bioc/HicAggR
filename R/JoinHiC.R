
#' Merge HiC chunk.
#'
#' JoinHiC
#' @description Create mega contactMatrix from a list of contactMatrix.
#' @param hicLst <List[contactMatrix]>: The HiC maps list.
#' @return A contactMatrix.
#' @examples
#' data(HiC_Ctrl.cmx_lst)
#' Mega_Ctrl.cmx <- JoinHiC(HiC_Ctrl.cmx_lst)
#'
JoinHiC <- function(
    hicLst
) {
    chromSizes <- attributes(hicLst)$chromSize
    mapNames <- names(hicLst)
    mapPosition <- data.frame(apply(
        mapNames |>
            strsplit("_") |>
            simplify2array() |>
            t(),
        seq_len(2),
        function(ele) {which(chromSizes$name == ele)}
    ))
    colnames(mapPosition) <- c("i", "j")
    megaDim <- sum(chromSizes$dimension)
    mega.dtf <- data.frame(
        i = NULL, j = NULL,
        x = NULL
    )
    for (i in seq_along(chromSizes$name)) {
        addI <- 0
        if (i > 1) { addI <- sum(chromSizes$dimension[seq_len(i) -1]) }
        for (j in seq(i, length(chromSizes$name))) {
            addJ <- 0
            if (j > 1) { addJ <- sum(chromSizes$dimension[seq_len(j) -1]) }
            if (length(which(mapPosition$i == i & mapPosition$j == j))) {
                map.dtf <- MeltSpm(
                    hicLst[[
                        which(mapPosition$i == i & mapPosition$j == j )
                    ]]@matrix
                )
                map.dtf$i <- map.dtf$i + addI
                map.dtf$j <- map.dtf$j + addJ
                mega.dtf <- rbind(mega.dtf, map.dtf)
            }
        }
    }
    mega.spm <- Matrix::sparseMatrix(
        i = mega.dtf$i,
        j = mega.dtf$j,
        x = mega.dtf$x,
        dims = c(megaDim, megaDim)
    )
    seqlengths.lst <- chromSizes$length |>
        stats::setNames(chromSizes$name)
    binnedGenome.gnr <- GenomicRanges::tileGenome(
        seqlengths.lst,
        tilewidth = attributes(hicLst)$resolution,
        cut.last.tile.in.chrom = TRUE
    )
    megaHic <- InteractionSet::ContactMatrix(
        mega.spm,
        binnedGenome.gnr,
        binnedGenome.gnr
    )
    megaHic@metadata <- attributes(hicLst)[-which(
        names(attributes(hicLst)) == "names"
    )]
    megaHic@metadata$kind <- "U"
    megaHic@metadata$symmetric <- TRUE
    return(megaHic)
}
