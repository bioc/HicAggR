#' Compute HiC matrix-balancing.
#'
#' BalanceHiC
#' @description Apply a matrix-balancing normalization method to a
#'  list of contacts matrix.
#' @param hicLst <List[ContactMatrix][InteractionSet::ContactMatrix()]>:
#'  The HiC maps list.
#' @param method <character> : The kind of normalization method.
#'  One of "ICE", "VC" or "VC_SQRT" (Default "ICE")
#' @param interactionType <character> : "cis", "trans", c("cis", "trans"),
#'  "all".
#'  If NULL normalization is apply on cis contactMatrix then trans
#'  contactMatrix (equivalent to c("cis", "trans")). If is "all", normalization
#'  is apply on all contactMatrix at once. (Default NULL)
#' @param maxIter <numerical>: The maximum iteration number.
#' @param qtlTh <numerical>: The threshold quantile below which the bins will
#'  be ignored. (Default 0.15)
#' @param cores <numerical> : Number of cores to be used. (Default 1)
#' @param verbose <logical>: If TRUE show the progression in console.
#'  (Default FALSE)
#' @return A matrices list.
#' @examples
#' data(HiC_Ctrl.cmx_lst)
#'
#' HiC_Ctrl_ICE.cmx_lst <- BalanceHiC(HiC_Ctrl.cmx_lst,
#'     interactionType = "cis",
#'     method = "ICE"
#' )
#'
#' HiC_Ctrl_VC.cmx_lst <- BalanceHiC(HiC_Ctrl.cmx_lst,
#'     interactionType = c("cis", "trans"),
#'     method = "VC"
#' )
#'
#' HiC_Ctrl_VC_SQRT.cmx_lst <- BalanceHiC(HiC_Ctrl.cmx_lst,
#'     interactionType = "all",
#'     method = "VC_SQRT"
#' )
#'
BalanceHiC <- function(
    hicLst, method = "ICE", interactionType = NULL, maxIter = 50,
    qtlTh = 0.15, cores = 1, verbose = FALSE
) {
    if (!is.null(interactionType) &&
        "all" %in% interactionType
    ) {
        megaHic <- JoinHiC(hicLst)
        if (method == "VC") {
            megaHic <- VCnorm(
                megaHic,
                qtlTh = qtlTh,
                vcsqrt = FALSE
            )
        } else if (method == "VC_SQRT") {
            megaHic <- VCnorm(
                megaHic,
                qtlTh = qtlTh,
                vcsqrt = TRUE
            )
        } else if (method == "ICE") {
            megaHic <- ICEnorm(
                megaHic,
                qtlTh = qtlTh,
                maxIter = maxIter
            )
        }
        hicLst <- CutHiC(megaHic, verbose = verbose)
    } else if (!is.null(interactionType) &&
        "cis" %in% interactionType &&
        NotIn("trans", interactionType)
    ) {
        matricesKind.tbl <- attributes(hicLst)$matricesKind
        if ("cis" %in% matricesKind.tbl$type) {
            cisMatricesNames.chr <- dplyr::filter(
                    matricesKind.tbl,
                    matricesKind.tbl$type == "cis"
                ) |>
                dplyr::pull("name")
            if (length(cisMatricesNames.chr)) {
                multicoreParam <- MakeParallelParam(
                    cores = cores,
                    verbose = verbose
                )
                hicLst[cisMatricesNames.chr] <- BiocParallel::bplapply(
                    BPPARAM = multicoreParam, seq_along(cisMatricesNames.chr),
                    function(ele.ndx) {
                        matrixName.chr <- cisMatricesNames.chr[[ele.ndx]]
                        if (method == "VC") {
                            hic <- VCnorm(
                                hicLst[[matrixName.chr]],
                                qtlTh = qtlTh,
                                vcsqrt = FALSE
                            )
                        } else if (method == "VC_SQRT") {
                            hic <- VCnorm(
                                hicLst[[matrixName.chr]],
                                qtlTh = qtlTh,
                                vcsqrt = TRUE
                            )
                        } else if (method == "ICE") {
                            hic <- ICEnorm(
                                hicLst[[matrixName.chr]],
                                qtlTh = qtlTh,
                                maxIter = maxIter
                            )
                        }
                        return(hic)
                    }
                )
            }
            transMatricesNames.chr <- dplyr::filter(
                matricesKind.tbl,
                matricesKind.tbl$type == "trans") |>
                dplyr::pull("name")
            mess.chr <- paste0(
                paste(transMatricesNames.chr, collapse = ", "),
                " removed from output."
            )
            if(verbose){
                message(mess.chr)
            }
            if (length(transMatricesNames.chr)) {
                attr.lst <- attributes(hicLst)
                attr.lst$matricesKind <- dplyr::filter(
                    attr.lst$matricesKind,
                    NotIn(attr.lst$matricesKind$name, transMatricesNames.chr)
                )
                chroms.chr <- attr.lst$matricesKind$name |>
                    strsplit("_") |>
                    unlist() |>
                    unique()
                attr.lst$chromSize <- dplyr::filter(
                    attr.lst$chromSize,
                    attr.lst$chromSize$name == chroms.chr
                )
                hicLst <- hicLst[-which(
                    names(hicLst) %in% transMatricesNames.chr
                )] |>
                    AddAttr(attrs=attr.lst)
            }
        } else {
            mess.chr <-paste0(
                "No cis matrix, Normalisation ",
                "won't be applied on cis matrices"
            )
            if(verbose){
                message(mess.chr)
            }
        }
    } else if (!is.null(interactionType) &&
        "trans" %in% interactionType && NotIn("cis", interactionType)
    ) {
        matricesKind.tbl <- attributes(hicLst)$matricesKind
        if ("trans" %in% matricesKind.tbl$type) {
            transMatricesNames.chr <- dplyr::filter(
                    matricesKind.tbl,
                    matricesKind.tbl$type == "trans") |>
                dplyr::pull("name")
            chromNames.chr <- transMatricesNames.chr |>
                strsplit("_") |>
                unlist() |>
                unique()
            chromSize.tbl <- attributes(hicLst)$chromSize
            chromSize.tbl <- dplyr::filter(
                chromSize.tbl,
                chromSize.tbl$name %in% chromNames.chr
            )
            trans.cmx_lst <- hicLst[transMatricesNames.chr] |>
                AddAttr(attrs = list(
                    resolution = attributes(hicLst)$resolution,
                    chromSize = chromSize.tbl,
                    matricesKind = matricesKind.tbl
                ))
            megaHic <- JoinHiC(trans.cmx_lst)
            if (method == "VC") {
                megaHic <- VCnorm(
                    megaHic,
                    qtlTh = qtlTh,
                    vcsqrt = FALSE
                )
            } else if (method == "VC_SQRT") {
                megaHic <- VCnorm(
                    megaHic,
                    qtlTh = qtlTh,
                    vcsqrt = TRUE
                )
            } else if (method == "ICE") {
                megaHic <- ICEnorm(
                    megaHic,
                    qtlTh = qtlTh,
                    maxIter = maxIter
                )
            }
            trans.cmx_lst <- CutHiC(megaHic, verbose = verbose)
            hicLst[transMatricesNames.chr] <-
                trans.cmx_lst[transMatricesNames.chr]
            cisMatricesNames.chr <- dplyr::filter(
                    matricesKind.tbl,
                    matricesKind.tbl$type == "cis"
                ) |>
                dplyr::pull("name")
            mess.chr <- paste0(
                paste(cisMatricesNames.chr, collapse = ", "),
                " removed from output."
            )
            if(verbose){
                message(mess.chr)
            }
            if (length(cisMatricesNames.chr)) {
                attr.lst <- attributes(hicLst)
                attr.lst$matricesKind <- dplyr::filter(
                    attr.lst$matricesKind,
                    NotIn(attr.lst$matricesKind$name, cisMatricesNames.chr)
                )
                chroms.chr <- attr.lst$matricesKind$name |>
                    strsplit("_") |>
                    unlist() |>
                    unique()
                attr.lst$chromSize <- dplyr::filter(
                    attr.lst$chromSize,
                    attr.lst$chromSize$name == chroms.chr
                )
                hicLst <- hicLst[-which(
                    names(hicLst) %in% cisMatricesNames.chr
                )] |>
                    AddAttr(attrs = attr.lst)
            }
        } else {
            mess.chr <- paste0(
                "No trans matrix, Normalisation ",
                "won't be applied on trans matrices"
            )
            if(verbose){
                message(mess.chr)
            }
        }
    } else {
        matricesKind.tbl <- attributes(hicLst)$matricesKind
        if (is.null(interactionType) ||
            "cis" %in% interactionType
        ) {
            if ("cis" %in% matricesKind.tbl$type) {
                cisMatricesNames.chr <- dplyr::filter(
                        matricesKind.tbl,
                        matricesKind.tbl$type == "cis"
                    ) |>
                    dplyr::pull("name")
                if (length(cisMatricesNames.chr)) {
                    multicoreParam <- MakeParallelParam(
                        cores = cores,
                        verbose = verbose
                    )
                    hicLst[cisMatricesNames.chr] <-
                        BiocParallel::bplapply(
                            BPPARAM = multicoreParam,
                            seq_along(cisMatricesNames.chr),
                            function(ele.ndx) {
                            matrixName.chr <- cisMatricesNames.chr[[ele.ndx]]
                            if (method == "VC") {
                                hic <- VCnorm(
                                    hicLst[[matrixName.chr]],
                                    qtlTh = qtlTh,
                                    vcsqrt = FALSE
                                )
                            } else if (method == "VC_SQRT") {
                                hic <- VCnorm(
                                    hicLst[[matrixName.chr]],
                                    qtlTh = qtlTh,
                                    vcsqrt = TRUE
                                )
                            } else if (method == "ICE") {
                                hic <- ICEnorm(
                                    hicLst[[matrixName.chr]],
                                    qtlTh = qtlTh,
                                    maxIter = maxIter
                                )
                            }
                            return(hic)
                            }
                    )
                }
            } else {
                mess.chr <- paste0(
                    "No cis matrix, Normalisation ",
                    "won't be applied on cis matrices"
                )
                if(verbose){
                    message(mess.chr)
                }
            }
        }
        if (is.null(interactionType) ||
            "trans" %in% interactionType
        ) {
            if ("trans" %in% matricesKind.tbl$type) {
                transMatricesNames.chr <- dplyr::filter(
                        matricesKind.tbl,
                        matricesKind.tbl$type == "trans"
                    ) |>
                    dplyr::pull("name")
                chromNames.chr <- transMatricesNames.chr |>
                    strsplit("_") |>
                    unlist() |>
                    unique()
                chromSize.tbl <- attributes(hicLst)$chromSize
                chromSize.tbl <- dplyr::filter(
                    chromSize.tbl,
                    chromSize.tbl$name %in% chromNames.chr
                )
                trans.cmx_lst <- hicLst[transMatricesNames.chr] |>
                    AddAttr(attrs = list(
                        resolution = attributes(hicLst)$resolution,
                        chromSize = chromSize.tbl,
                        matricesKind = matricesKind.tbl
                    ))
                megaHic <- JoinHiC(trans.cmx_lst)
                if (method == "VC") {
                    megaHic <- VCnorm(
                        megaHic,
                        qtlTh = qtlTh,
                        vcsqrt = FALSE
                    )
                } else if (method == "VC_SQRT") {
                    megaHic <- VCnorm(
                        megaHic,
                        qtlTh = qtlTh,
                        vcsqrt = TRUE
                    )
                } else if (method == "ICE") {
                    megaHic <- ICEnorm(
                        megaHic,
                        qtlTh = qtlTh,
                        maxIter = maxIter
                    )
                }
                trans.cmx_lst <- CutHiC(megaHic, verbose = verbose)
                hicLst[transMatricesNames.chr] <-
                    trans.cmx_lst[transMatricesNames.chr]
            } else {
                mess.chr <- paste0(
                    "No trans matrix, Normalisation ",
                    "won't be applied on trans matrices"
                )
                if(verbose){
                    message(mess.chr)
                }
            }
        }
    }
    return(hicLst)
}
