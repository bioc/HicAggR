#' Genomic distance bias correction.
#'
#' OverExpectedHiC
#' @description Function that normalises HiC matrices by expected values
#'  computed per genomic distance.
#' @param hicLst <List[ContactMatrix][InteractionSet::ContactMatrix()]>:
#'  The HiC maps list.
#' @param method <character> Options are "mean_non_zero",
#' "mean_total", or "lieberman". Look at details for more.
#' (Default: "mean_non_zero")
#' @param cores <numerical> : Number of cores to be used.
#' (Default 1)
#' @param verbose <logical>: Show the progression in console?
#'  (Default FALSE)
#' @param plot_contact_vs_dist Whether to plot
#' contact vs distance curve per chromosome ("per_seq"),
#' all chromosomes ("total") or not (NULL). (Default "per_seq")
#' @return A matrices list.
#' @export
#' @importFrom rlang .data
#' @import ggplot2
#' @importFrom reshape melt
#' @importFrom checkmate checkChoice
#' @details
#' Methods to calculate expected values per distance:
#' \itemize{
#' \item "mean_non_zero": for each distance, average contact value is
#' calculated using only non-zero values.
#' \item "mean_total": for each distance, average contact value is
#' calculated using all values at this distance.
#' \item "lieberman": for each distance, contact values are summed and
#' divided by chromsome length minus distance. Only for cis contacts.
#' }
#' @examples
#' # Note: run HicAggR::BalanceHiC before OverExpectedHiC calculation.
#' data(HiC_Ctrl.cmx_lst)
#' OverExpectedHiC(HiC_Ctrl.cmx_lst)
#'
OverExpectedHiC <- function(
    hicLst, method = "mean_non_zero",
    verbose = FALSE, cores = 1, 
    plot_contact_vs_dist="per_seq"
) {
    .validHicMatrices(matrices = hicLst)
    checkmate::checkChoice(
        x = method,
        choices = c("mean_non_zero", "mean_total", "lieberman"),
        null.ok = FALSE
    )
    resolution.num <- attributes(hicLst)$resolution
    matricesNames.chr <- names(hicLst)
    multicoreParam <- MakeParallelParam(
        cores = cores,
        verbose = verbose
    )
    expected.lst <- BiocParallel::bplapply(
        BPPARAM = multicoreParam, seq_along(matricesNames.chr),
        function(matrixName.ndx) {
            matrixName.chr <- matricesNames.chr[matrixName.ndx]
            if (hicLst[[matrixName.chr]]@metadata$type == "cis") {
                hic.dtf <- MeltSpm(hicLst[[matrixName.chr]]@matrix)
                if(method != "mean_non_zero"){
                    hic_dtf_total <- 
                    as.matrix(hicLst[[matrixName.chr]]@matrix)|>
                    reshape::melt()|>
                    `colnames<-`(c("i","j","x"))|>
                    dplyr::filter(.data$i <= .data$j) |>
                    dplyr::mutate(distance =
                        1 + (.data$j - .data$i) * resolution.num)|>
                        dplyr::select(c("x", "distance"))
                }
                hic.dtf <- dplyr::mutate(
                    hic.dtf,
                    distance = 
                        (1 + (hic.dtf$j - hic.dtf$i) * resolution.num)) |>
                    dplyr::select(c("x", "distance"))
                if(method == "lieberman"){
                    chr_length <- attributes(hicLst)$chromSize$length[which(
                        attributes(hicLst)$chromSize$name == 
                        unlist(strsplit(matrixName.chr,split = "_"))[1])]

                    expected.dtf <- dplyr::group_by(
                        hic_dtf_total, distance = hic_dtf_total$distance) |>
                        # get distance, chr_length and contacts!!!
                        dplyr::summarise_at(.vars="x",
                            .funs = list(total=sum)) |>
                        dplyr::mutate(factor=chr_length - .data$distance) |>
                        dplyr::mutate(expected=.data$total/.data$factor) |>
                        dplyr::select(-c("factor","total"))

                    }else if(method == "mean_total"){
                    expected.dtf <- dplyr::group_by(
                        hic_dtf_total, distance = hic_dtf_total$distance) |>
                        dplyr::summarise_at(
                        .vars = "x", .funs = list(expected = mean)
                        )
                    }else{
                expected.dtf <- dplyr::group_by(
                    hic.dtf, distance = hic.dtf$distance) |>
                    dplyr::summarise_at(
                        .vars = "x", .funs = list(expected = mean)
                    )
                    }
                expected.num <- dplyr::left_join(
                    hic.dtf,
                    expected.dtf,
                    by = "distance"
                )
            } else {
                expected.num <- list(
                    expected = sum(hicLst[[matrixName.chr]]@matrix@x)/
                        Reduce(`*`, hicLst[[matrixName.chr]]@matrix@Dim)
                )
            }
            expected.lst <- list(expected.num) |>
                stats::setNames(matrixName.chr)
            return(expected.lst)
        }
    ) |> do.call(what = "c")
    cisNames.chr <- NULL
    for (matrixName.chr in matricesNames.chr) {
        hicLst[[matrixName.chr]]@metadata$expected <-
            expected.lst[[matrixName.chr]]$expected
        hicLst[[matrixName.chr]]@matrix@x <-
            hicLst[[matrixName.chr]]@matrix@x/
            expected.lst[[matrixName.chr]]$expected
        if (hicLst[[matrixName.chr]]@metadata$type == "cis") {
            cisNames.chr <- c(cisNames.chr, matrixName.chr)
        }
    }
    attributes(hicLst)$mtx <- "o/e"
    expected.dtf <- expected.lst[cisNames.chr] |>
        do.call(what = rbind)
    expected.lst <- expected.lst[cisNames.chr]

    if(plot_contact_vs_dist=="per_seq" && length(expected.lst)>0){
        data <- do.call(rbind,
            lapply(
                seq_along(expected.lst),
                function(i){return(expected.lst[[i]] |>
                dplyr::filter(.data$distance>1) |>
                dplyr::filter(
                !duplicated(.data$distance)|!duplicated(.data$expected)) |>
                        dplyr::mutate(seqnames=names(expected.lst)[i]))
            })
        )
        p <- ggplot2::ggplot(data,
            ggplot2::aes(x=log10(.data$distance),
                y=log10(.data$expected),col=.data$seqnames))+
            ggplot2::geom_line()+
            ggplot2::ggtitle(label = "Contact vs distance (per chromosome)")+
            ggplot2::theme_bw()
            plot(p)
        }
    attributes(hicLst)$expected <- dplyr::group_by(
        expected.dtf,
        distance = expected.dtf$distance) |>
        dplyr::summarise_at(.vars = "expected", .funs = list(expected = mean))
        if(plot_contact_vs_dist=="total" && length(expected.lst)>0){
        data <- attributes(hicLst)$expected |>
            dplyr::filter(.data$distance>1)

        p <- ggplot2::ggplot(data,
            ggplot2::aes(x=log10(.data$distance),
                y=log10(.data$expected),col="#619CFF"))+
            ggplot2::geom_line(show.legend = FALSE)+
            ggplot2::ggtitle(label = "Contact vs distance")+
            ggplot2::theme_bw()
        plot(p)
        }
    return(hicLst)
}