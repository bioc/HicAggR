#' Indexes GRanges on genome.
#'
#' IndexFeatures
#' @description Function that indexes a GRanges object on binned genome and constraints. Needed prior HicAggR::SearchPairs() function.
#' @param gRangeList <GRanges or GRangesList or list[GRanges]>: GRanges object, list of GRanges or GRangesList containing coordinates to index.
#' @param genomicConstraint <GRanges>: GRanges object of constraint regions. Note that bins in the same constraint region only will be paired in HicAggR::SearchPairs(). If NULL chromosomes in chromSizes are used as constraints (Default NULL)
#' @param chromSizes <data.frame>: A data.frame containing chromosomes names and lengths in base pairs (see example).
#' @param binSize <integer>: Bin size in bp - corresponds to HiC matrix resolution.
#' @param metadataColName <character> : A character vector that specify the metadata columns of GRanges on which apply the summary method if multiple ranges are indexed in the same bin.
#' @param method <character>: A string defining which summary method is used on metadata columns defined in metadataColName if multiple ranges are indexed in the same bin. Use 'mean', 'median', 'sum', 'max' or 'min'. (Default 'mean'')
#' @param cores <integer> : Number of cores used. (Default 1)
#' @param verbose <logical>: If TRUE show the progression in console. (Default FALSE)
#' @return A GRanges object.
#' @examples
#' data(Beaf32_Peaks.gnr)
#' Beaf32_Index.gnr <- IndexFeatures(
#'     gRangeList = list(Beaf = Beaf32_Peaks.gnr),
#'     chromSizes = data.frame(
#'         seqnames = c("2L", "2R"),
#'         seqlengths = c(23513712, 25286936)
#'     ),
#'     binSize = 100000
#' )
#'
IndexFeatures <- function(
    gRangeList = NULL, genomicConstraint = NULL, chromSizes = NULL,
    binSize = NULL, method = "mean", metadataColName = NULL,
    cores = 1, verbose = FALSE
) {
    # Constraint Informations
    if (is.null(genomicConstraint)) {
        genomicConstraint <- GenomicRanges::GRanges(
            seqnames = chromSizes[, 1],
            ranges = IRanges::IRanges(
                start = rep(1, length(chromSizes[, 2])),
                end = chromSizes[, 2]
            ),
            strand = "*", name = chromSizes[, 1]
        )
    } else {
        if (is.null(genomicConstraint$name) |
            length(which(!is.na(genomicConstraint$name))) == 0) {
            genomicConstraint$name <- paste0(
                "Constraint_",
                seq_along(genomicConstraint)
            )
        }
    }
    seqLevelsStyle.chr <- GenomeInfoDb::seqlevelsStyle(genomicConstraint)
    if (length(seqLevelsStyle.chr) > 1) {
        seqLevelsStyle.chr <- seqLevelsStyle.chr[[1]]
        GenomeInfoDb::seqlevelsStyle(genomicConstraint) <- seqLevelsStyle.chr
    }
    binnedConstraint.gnr <- BinGRanges(
        gRange = genomicConstraint,
        chromSizes = chromSizes,
        binSize = binSize,
        verbose = verbose,
        reduceRanges = FALSE,
        cores = cores
    )
    # Feature Names
    if (inherits(gRangeList, "GRanges")) {
        gRangeList <- list(Features = gRangeList)
    } else if (inherits(gRangeList, "GRangesList")) {
        gRangeList <- as.list(gRangeList)
    }
    if (is.null(names(gRangeList))) {
        gRangeList <- stats::setNames(
            gRangeList,
            paste0("Feature_", seq_along(gRangeList))
        )
    }
    gRangeOrder.ndx <- lapply(gRangeList, length) |>
        unlist() |>
        order(decreasing = TRUE)
    gRangeList <- gRangeList[gRangeOrder.ndx]
    feature.chr_vec <- names(gRangeList)
    # GRanges Binning
    multicoreParam <- MakeParallelParam(cores = cores, verbose = FALSE) # DD 230310 correct how parallezation was done...
    binnedFeature.lst <- BiocParallel::bplapply(BPPARAM = multicoreParam,
        seq_along(gRangeList),function(feature.ndx) {
            feature.chr <- feature.chr_vec[[feature.ndx]]
            feature.gnr <- IRanges::subsetByOverlaps(
                gRangeList[[feature.chr]],
                genomicConstraint
            )
            # adding name metadata if the GRanges do not have name
            if(is.null(feature.gnr$name)){
                if(!is.null(names(gRangeList)[feature.ndx])){
                    feature.gnr$name = paste0(names(gRangeList)[feature.ndx],'_',seq(1,length(feature.gnr)))
                }else{
                    feature.gnr$name = paste0("peak_feature",feature.ndx,'_',seq(1,length(feature.gnr)))
                }
            }
            GenomeInfoDb::seqlevelsStyle(feature.gnr) <- seqLevelsStyle.chr
            binnedFeature.gnr <- BinGRanges(
                gRange = feature.gnr, chromSizes = chromSizes,
                binSize = binSize, method = method,
                metadataColName = metadataColName,
                verbose = verbose, reduceRanges = TRUE,
                cores = cores
            )
            binnedFeat.tbl <- tibble::tibble(
                BinnedFeature.ndx = seq_along(binnedFeature.gnr),
                Feature.name = binnedFeature.gnr$name
            ) |>
                tidyr::unnest(cols = "Feature.name")
            binnedFeat.tbl <- dplyr::group_by(
                binnedFeat.tbl,
                Feature.name = binnedFeat.tbl$Feature.name
            )
            binnedFeat.tbl <- tidyr::nest(binnedFeat.tbl) |>
                stats::setNames(c("Feature.name", "BinnedFeature.ndx"))
            binnedConstraint.tbl <- tibble::tibble(
                BinnedConstraint.ndx = seq_along(binnedConstraint.gnr),
                Constraint.name = binnedConstraint.gnr$name
            )
            binnedConstraint.tbl <- dplyr::group_by(
                binnedConstraint.tbl,
                Constraint.name = binnedConstraint.tbl$Constraint.name
            )
            binnedConstraint.tbl <- tidyr::nest(binnedConstraint.tbl) |>
                stats::setNames(c("Constraint.name", "BinnedConstraint.ndx"))
            featConstOvlp.ovlp <- GenomicRanges::findOverlaps(
                feature.gnr,
                genomicConstraint
            )
            featConstOvlp.tbl <- tibble::tibble(
                Feature.name = feature.gnr$name[featConstOvlp.ovlp@from],
                Constraint.name = genomicConstraint$name[featConstOvlp.ovlp@to]
            ) |>
                dplyr::left_join(binnedFeat.tbl, by = "Feature.name") |>
                dplyr::select(-"Feature.name") |>
                tidyr::unnest(cols = "BinnedFeature.ndx") |>
                unique()
            featConstOvlp.tbl <- dplyr::group_by(
                featConstOvlp.tbl,
                Constraint.name = featConstOvlp.tbl$Constraint.name
            )
            featConstOvlp.tbl <- tidyr::nest(featConstOvlp.tbl) |>
                stats::setNames(c("Constraint.name", "BinnedFeature.ndx")) |>
                dplyr::left_join(binnedConstraint.tbl, by = "Constraint.name")

            binnedFeature.gnr_lst <- lapply(seq_len(nrow(featConstOvlp.tbl)),
                function(row.ndx) {
                    ranges.ndx <-
                        featConstOvlp.tbl$BinnedFeature.ndx[row.ndx] |>
                        unlist(use.names = FALSE)
                    constraint.ndx <-
                        featConstOvlp.tbl$BinnedConstraint.ndx[row.ndx] |>
                        unlist(use.names = FALSE)
                    subBinnedFeature.gnr <- IRanges::subsetByOverlaps(
                        binnedFeature.gnr[ranges.ndx],
                        binnedConstraint.gnr[constraint.ndx]
                    )
                    subBinnedFeature.gnr$constraint <-
                        featConstOvlp.tbl$Constraint.name[row.ndx]
                    return(subBinnedFeature.gnr)
                }
            )
            binnedFeature.gnr <- MergeGRanges(
                binnedFeature.gnr_lst,
                sortRanges = FALSE,
                reduceRanges = FALSE
            )
            binnedFeature.gnr$bln <- TRUE
            names(S4Vectors::mcols(binnedFeature.gnr)) <- paste0(
                feature.chr,
                ".",
                names(S4Vectors::mcols(binnedFeature.gnr))
            )
            names(S4Vectors::mcols(binnedFeature.gnr))[which(
                names(S4Vectors::mcols(binnedFeature.gnr)) ==
                    paste0(feature.chr, ".bin")
            )] <- "bin"
            names(S4Vectors::mcols(binnedFeature.gnr))[which(
                names(S4Vectors::mcols(binnedFeature.gnr)) ==
                    paste0(feature.chr, ".constraint")
            )] <- "constraint"
            binnedFeature.gnr$name <- paste0(binnedFeature.gnr$bin,
                ":",
                binnedFeature.gnr$constraint
            )
            metadataBinnedFeature.dtf <- data.frame(
                S4Vectors::mcols(binnedFeature.gnr)
            )
            S4Vectors::mcols(binnedFeature.gnr) <- NULL
            return(list(
                binnedFeature.gnr = binnedFeature.gnr,
                featureMetadata.dtf = metadataBinnedFeature.dtf
            ))
        }
    )
    binnedIndex.gnr <- binnedFeature.lst |>
        lapply("[[", "binnedFeature.gnr") |>
        MergeGRanges(sortRanges = FALSE, reduceRanges = FALSE)
    featureMetadata.lst_dtf <- binnedFeature.lst |>
        lapply("[[", "featureMetadata.dtf")
    S4Vectors::mcols(binnedIndex.gnr) <- BindFillRows(featureMetadata.lst_dtf)
    ids.lst <- binnedIndex.gnr$name
    dupplicatedIds.lst <- unique(ids.lst[duplicated(ids.lst)])
    idDuplicated.ndx <- which(ids.lst %in% dupplicatedIds.lst)
    # Merge GRanges into one index and duplicated Bin handle
    if (length(idDuplicated.ndx)) {
        binnedIndexDuplicated.dtf <- data.frame(
            binnedIndex.gnr[idDuplicated.ndx]
        )
        binnedIndexDuplicated.tbl <- tibble::tibble(binnedIndexDuplicated.dtf)
        binnedIndexDuplicated.tbl <- dplyr::group_by(
            binnedIndexDuplicated.tbl,
            name = binnedIndexDuplicated.tbl$name
        )
        binnedIndexDuplicated.tbl <- tidyr::nest(binnedIndexDuplicated.tbl)
        binnedIndexNoDuplicated.dtf <- data.frame(
            binnedIndex.gnr[-idDuplicated.ndx]
        )
        binnedIndexNoDuplicated.tbl <- tibble::tibble(
            binnedIndexNoDuplicated.dtf
        )
        multicoreParam <- MakeParallelParam(
            cores = cores,
            verbose = verbose
        )
        binnedIndexDuplicated.lst <- BiocParallel::bplapply(
            BPPARAM = multicoreParam, seq_len(nrow(binnedIndexDuplicated.tbl)),
            function(row.ndx) {
                rowName.chr <- binnedIndexDuplicated.tbl$name[[row.ndx]]
                row <- binnedIndexDuplicated.tbl$data[[row.ndx]]
                col.lst <- lapply(
                    seq_along(row),
                    function(col.ndx) {
                        col <- dplyr::pull(row, col.ndx)
                        if (length(unique(stats::na.omit(col))) == 0) {
                            return(NA)
                        } else if (length(unique(stats::na.omit(col))) == 1) {
                            return(unique(stats::na.omit(col)))
                        } else {
                        return(list(unlist(col)))
                        }
                    }
                ) |>
                    stats::setNames(names(row))
                binnedIndexDuplicated.tbl <- tibble::as_tibble(col.lst) |>
                    tibble::add_column(name = rowName.chr)
                return(binnedIndexDuplicated.tbl)
            }
        )
        binnedIndexDuplicated.tbl <- dplyr::bind_rows(
            binnedIndexDuplicated.lst
        )
        binnedIndex.gnr <- rbind(
                binnedIndexDuplicated.tbl,
                binnedIndexNoDuplicated.tbl
            ) |>
            data.frame() |>
            methods::as("GRanges")
    }
    for (featureName.chr in feature.chr_vec) {
        colname.chr <- paste0(featureName.chr, ".bln")
        S4Vectors::mcols(binnedIndex.gnr)[
            which(is.na(S4Vectors::mcols(binnedIndex.gnr)[,colname.chr])),
            colname.chr
        ] <- FALSE
        S4Vectors::mcols(binnedIndex.gnr)[,colname.chr] <- methods::as(
            S4Vectors::mcols(binnedIndex.gnr)[,colname.chr], "Rle"
        )
    }
    columOrder.chr <- c(
            "name",
            "bin",
            "constraint",
            names(S4Vectors::mcols(binnedIndex.gnr))
        ) |>
        unique()
    S4Vectors::mcols(binnedIndex.gnr) <- as.data.frame(
            S4Vectors::mcols(binnedIndex.gnr)
        ) |>
        dplyr::select(dplyr::all_of(columOrder.chr)) # should avoid warning tidyselect
    # When indexing a GRangeList seqinfo gets lost somehow
    chromSizes = chromSizes[which(chromSizes[,1]!="All"),]
    GenomeInfoDb::seqinfo(binnedIndex.gnr) = GenomeInfoDb::Seqinfo(seqnames = chromSizes[,1],seqlengths = chromSizes[,2],)
    return(sort(binnedIndex.gnr))
}
