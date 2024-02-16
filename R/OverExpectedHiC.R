#' Genomic distance bias correction.
#'
#' OverExpectedHiC
#' @description Function that normalises HiC matrices by expected values
#'  computed per genomic distance.
#' @param hicLst <List[ContactMatrix][InteractionSet::ContactMatrix()]>:
#'  The HiC maps list.
#' @param cores <numerical> : Number of cores to be used. (Default 1)
#' @param verbose <logical>: If TRUE show the progression in console.
#'  (Default FALSE)
#'  @param plot_contact_vs_dist Whether to plot contact vs distance curve
#'  per chromosome ("per_seq"), all chromosomes ("total") or not (NULL). 
#'  (Default "per_seq")
#' @return A matrices list.
#' @import ggplot2
#' @examples
#' # Note: run HicAggR::BalanceHiC before OverExpectedHiC calculation.
#' data(HiC_Ctrl.cmx_lst)
#' OverExpectedHiC(HiC_Ctrl.cmx_lst)
#'
OverExpectedHiC <- function(
    hicLst, verbose = FALSE, cores = 1, plot_contact_vs_dist="per_seq"
) {
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
                hic.dtf <- dplyr::mutate(
                    hic.dtf,
                    distance = 1 + (hic.dtf$j - hic.dtf$i) * resolution.num) |>
                    dplyr::select(c("x", "distance"))
                expected.dtf <- dplyr::group_by(
                    hic.dtf, distance = hic.dtf$distance) |>
                    dplyr::summarise_at(
                        .vars = "x", .funs = list(expected = mean)
                    )
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
				dplyr::filter(distance>1) |>
				dplyr::filter(
				!duplicated(distance)|!duplicated(expected)) |>
						dplyr::mutate(seqnames=names(expected.lst)[i]))
			})
		)
    p <- ggplot2::ggplot(data,
		ggplot2::aes(x=log10(distance),
			y=log10(expected),col=seqnames))+
		ggplot2::geom_line()+
		ggplot2::ggtitle(label = "Contact vs distance (per chromosome)")+
		ggplot2::theme_bw()
    print(p)
  	}
    attributes(hicLst)$expected <- dplyr::group_by(
        expected.dtf,
        distance = expected.dtf$distance) |>
        dplyr::summarise_at(.vars = "expected", .funs = list(expected = mean))
  	if(plot_contact_vs_dist=="total" && length(expected.lst)>0){
    	data <- attributes(hicLst)$expected |>
      	dplyr::filter(distance>1)
                      
    p <- ggplot2::ggplot(data,
		ggplot2::aes(x=log10(distance),
			y=log10(expected),col="#619CFF"))+
		ggplot2::geom_line(show.legend = F)+
		ggplot2::ggtitle(label = "Contact vs distance")+
		ggplot2::theme_bw()
    print(p)
  	}
    return(hicLst)
}