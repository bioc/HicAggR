#' removes from background couples those found
#' in target couples
#' @keywords internal
#' @param bg_couples GInteractions object with
#' background couples
#' @param matrices hic matrices list with target couples
#' @return GInteractions object with clean background couples
#' @noRd
.cleanCouples <- function(
    bg_couples,
    matrices){
    # discard couples that are in target couples
    bg_couples <- bg_couples[
    which(!attr(bg_couples,"NAMES")%in%attr(matrices,"names"))]
    return(bg_couples)
}

#' generates background couples between
#' randomly picked bins and baits from
#' target couples
#' @keywords internal
#' @param targetCouples GInteractions object with
#' target couples
#' @param genomicConstraint GRanges object with
#' constraint coordinates
#' @param indexAnchor indexed anchors
#' @param indexBait indexed baits
#' @param chromSizes data.frame with chromosome
#' names and lengths
#' @param resolution numeric of hic resolution
#' @param cores number of cores
#' @param verbose verbosity
#' @importFrom GenomicRanges GRanges width trim
#' @importFrom IRanges IRanges subsetByOverlaps
#' @return GInteractions object with background couples
#' @noRd
.findRandomBins <- function(
    targetCouples,
    genomicConstraint=NULL,
    indexAnchor=NULL,
    indexBait=NULL,
    chromSizes,
    resolution,
    cores = cores,
    verbose = FALSE){
    
    if(verbose){
        message("random bins")
    }
    if(!is.null(genomicConstraint)){
    dist_const <- c(min(targetCouples$distance),
        max(GenomicRanges::width(genomicConstraint)))
    }else{
    dist_const <- c(min(targetCouples$distance),max(targetCouples$distance))
    }
    if (is.null(genomicConstraint)) {
    genomicConstraint <- GenomicRanges::GRanges(
        seqnames = as.character(chromSizes[, 1]),
        ranges = IRanges::IRanges(
        start = rep(1, length(chromSizes[, 2])),
        end = as.numeric(chromSizes[, 2])
        ),
        strand = "*", name = as.character(chromSizes[, 1])
    )
    }
    binned_constraint <- BinGRanges(genomicConstraint,
                                chromSizes,
                                binSize = resolution)
    binned_constraint$name <- paste0("random_",
        seq(1,length(binned_constraint)))
    binned_constraint <- IRanges::subsetByOverlaps(
    binned_constraint,
    (GenomicRanges::trim(unique(indexBait)+(resolution*2))),
    invert = TRUE)
    
    binned_constraint <- IRanges::subsetByOverlaps(
    binned_constraint,
    (GenomicRanges::trim(unique(indexAnchor)+(resolution*2))),
    invert = TRUE)
    
    binned_constraint.idx <- IndexFeatures(binned_constraint,
                genomicConstraint,
                chromSizes,
                binSize = resolution,cores=cores)
    bait.idx <- IndexFeatures(unique(anchors(targetCouples,type="second")),
                                        genomicConstraint,
                                        chromSizes,
                                        binSize = resolution,cores=cores)
    if(dist_const[1]==0){dist_const[1]<-resolution+1}
    background_pairs <- SearchPairs(binned_constraint.idx,
                                bait.idx,
                                minDist = dist_const[1],
                                maxDist = dist_const[2],cores=cores)
    if(verbose){
        message("Number of background pairs: ",
        length(background_pairs))
    }
    return(background_pairs)
}

#' generate inter-TAD couples
#' 
#' @keywords internal
#' @param targetCouples GInteractions object with
#' target couples
#' @param genomicConstraint GRanges object with
#' TAD coordinates
#' @param indexAnchor indexed anchors
#' @param indexBait indexed baits
#' @param chromSizes data.frame with chromosome
#' names and lengths
#' @param resolution numeric of hic resolution
#' @param secondaryConst.var column name
#' of secondary variable, such as compartiment
#' info
#' @param cores number of cores
#' @param verbose verbosity
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges width GRanges
#' @importFrom dplyr mutate filter
#' @importFrom S4Vectors mcols
#' @return GInteractions object with background couples
#' @noRd
.interTads <- function(
    targetCouples,
    genomicConstraint = NULL,
    indexAnchor = NULL,
    indexBait = NULL,
    chromSizes,
    resolution,
    secondaryConst.var = NULL,
    cores = cores,
    verbose = FALSE){
    if(verbose){
        message("inter-TADs")
    }
    if(!is.null(genomicConstraint) && 
    !all(unique(indexAnchor$constraint)%in%
        GenomeInfoDb::seqlevels(indexAnchor)) && 
    !all(unique(indexBait$constraint)%in%
        GenomeInfoDb::seqlevels(indexBait))){
    dist_const <- c(min(targetCouples$distance),
        max(GenomicRanges::width(genomicConstraint)))
    }else{
    stop("Inter-TAD background couples require TADs info!
        Anchor and Bait coordinates should be indexed with TADs info!")
    }
    indexAnchor.tmp <- indexAnchor |>
    as.data.frame() |>
    dplyr::mutate(tad_name = .data$constraint) |>
    dplyr::mutate(constraint = .data$seqnames)|>
    GenomicRanges::GRanges()
    
    indexBait.tmp <- indexBait |>
    as.data.frame() |>
    dplyr::mutate(tad_name = .data$constraint) |>
    dplyr::mutate(constraint = .data$seqnames)|>
    GenomicRanges::GRanges()
    
    bg_couples <- SearchPairs(
    indexAnchor = indexAnchor.tmp,
    indexBait = indexBait.tmp,
    minDist = dist_const[1],
    maxDist = dist_const[2],
    verbose = FALSE,
    cores = cores)
    
    diff.tad.couples <- as.data.frame(S4Vectors::mcols(bg_couples)) |>
    dplyr::filter(.data$anchor.tad_name!=.data$bait.tad_name) |>
    rownames()
    bg_couples <- bg_couples[
    which(attr(bg_couples,"NAMES")%in%diff.tad.couples)
    ]
    if(verbose){
        message("Number of inter-TAD couples ",length(bg_couples))
    }
    if(!is.null(secondaryConst.var)){
    ## discard couples in the same compartment
    ## names of compartment info in couples
    compart.var1 <- paste0("anchor.",secondaryConst.var)
    compart.var2 <- paste0("bait.",secondaryConst.var)
    diff.compart.couples <- as.data.frame(S4Vectors::mcols(bg_couples)) |>
        dplyr::filter(compart.var1!=compart.var2) |>
        rownames()
    bg_couples <- bg_couples[
        which(attr(bg_couples,"names")%in%diff.compart.couples)]
    }
    return(bg_couples)
}
#' Fit a polynomial model with 2 degrees on the background couples
#' log(counts)~distance. Get residuals mean and sd for the background,
#' predict target couples, calculate residuals for targets and compute
#' z.scores to perform z.test. Adjust p.values.
#'
#' @param df_bg data.frame of distance
#' to counts for the background couples.
#' @param df_target data.frame of distance
#' to counts for the target couples.
#' @param method_adjust method used to adjust
#' p.values. (Default: "BH")
#'
#' @return data.frame with interactions id, z.score, p.value and adjusted
#' p.values for each target couple
#' @importFrom stats residuals lm poly sd predict pnorm p.adjust as.formula
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @noRd
.computeZscore <- function(
    df_bg = NULL,
    df_target = NULL,
    method_adjust = "BH"
){
    model_poly_df2 <- stats::lm(formula = stats::as.formula(
        "log(counts)~stats::poly(distance,df=2)"),
        data = df_bg,
        subset = (df_bg$counts > 0))
    avg_residuals <- mean(stats::residuals(model_poly_df2))
    sd_residuals <- stats::sd(stats::residuals(model_poly_df2))
    target_predictions <- stats::predict(
        object = model_poly_df2,
        newdata = data.frame(
            distance = df_target$distance
        )
    )
    target_residuals <- df_target$counts - exp(target_predictions)
    z_output <- data.frame(
        name = df_target$names,
        z.score = (
            (target_residuals - avg_residuals)/sd_residuals)
        )|> dplyr::mutate(
            p.value = stats::pnorm(
                q = .data$z.score,
                lower.tail = FALSE
            )
        )|>
        dplyr::mutate(
            adj.p_val = stats::p.adjust(
                p = .data$p.value,
                method = method_adjust
            )
        )
    return(z_output)
}

#' CompareToBackground
#'
#' @description Computes z.test for each target couple over background couples.
#'
#' @param hicList <List[ContactMatrix][InteractionSet::ContactMatrix()]>:
#' The HiC maps list.
#' @param matrices <list[matrix]>: The matrices list.
#' @param genomicConstraint <GRanges>: GRanges object of
#' constraint regions. If `NULL`, chromosomes in chromSizes are used as
#' constraints (Default NULL)
#' @param secondaryConst.var <character>: A string defining
#' column name containing compartment information in the
#' metadata of anchor and bait <GRanges> objects. (Default NULL)
#' @param chromSizes <data.frame>: A data.frame containing chromosomes
#' names and lengths in base pairs.
#' @param n_background <integer> : Number of background
#' couples to keep. We recommend to use `set.seed` prior to launching
#' this function with non null n_background. (Default NULL)
#' @param indexAnchor <GRanges>: A first indexed GRanges object
#' used as pairs anchor (must be indexed using `IndexFeatures()`).
#' @param indexBait <GRanges>: A second indexed GRanges object
#' used as pairs bait (must be indexed using `IndexFeatures()`).
#' @param cores <integer> : Number of cores used. (Default 1)
#' @param areaFun <character or function>: A character
#' or function that allows to extract an area from each matrix that
#' composes the matrices list (Default "center").
#' Look at [GetQuantif] for more info.
#' @param operationFun <character or function>: A character
#' or function specifying the operation applied to the selected area.
#' (Default "mean").
#' Look at [GetQuantif] for more info.
#' @param bg_type <character>: Type of background couples
#' to generate. Possible choices: "random_anchors", "inter_TAD",
#' "inter_compartment", NULL (Defaults to "random_anchors").
#' More information in details...
#' @param verbose <logical> details on progress? (Default: FALSE)
#' @param p_adj_method <character> method used
#' to adjust p.values. More from `stats::p.adjust()`. (Default: "BH")
#' @param ... arguments to pass to [PrepareMtxList], inorder to treat
#' background matrices.
#' @return returns a <list> object with the z.test output for each
#' target couple, values for the target couples and values for the
#' background couples.
#' 
#' @export
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom IRanges subsetByOverlaps
#' @importFrom rlang .data
#' @importFrom dplyr group_by mutate summarise left_join
#' @importFrom S4Vectors mcols
#' @details
#' Types of background couples possible:
#' \itemize{
#' \item "random_anchors": picks random bins as anchors and forms couples 
#' with bait bins.
#' If genomicConstraint is supplied, only intra-TAD random-bait couples are
#'  kept.
#' Else intra-TAD random-bait couples within a distance constraint 
#' corresponding to the minimal and maximal distances of target couples.
#' \item "inter_TAD": If target couples were formed using TAD information 
#' with non `NULL`
#' genomicConstraint argument, then inter-TAD anchor-bait couples are used 
#' as background.
#' Distance constraint applied correspond to the minimal distance of target
#'  couples and
#' maximal width of supplied TADs.
#' \item "inter_compartment": If `secondaryConst.var` is not
#' `NULL` and both indexAnchor and indexBait objects contain the
#' provided variable name, then background couples are formed between anchors
#' and baits located in different compartments.
#' \item "NULL": If `NULL`, `random_anchors` are set by default.
#' }
#'
#' Notes on the comparison between bg and target couples:
#' We noticed that o/e values tend to be
#' skewed towards very long distance interactions. As a result, long
#' distance background couples tend to influence strongly mean and sd,
#' resulting in more long distance target couples being significant.
#' So rather than computing z.score over all background couples,
#' we've chosen to fit a polynomial with 2 degrees on the log(counts)
#' vs distance data of the background couples. Z.scores are then
#' computed per target couple by comparing residuals of the target counts
#' as predicted by the model and the residuals of the background couples.
#'
#' @examples
#' h5_path <- system.file("extdata",
#'     "Control_HIC_10k_2L.h5",
#'     package = "HicAggR", mustWork = TRUE
#' )
#' binSize=10000
#' data(Beaf32_Peaks.gnr)
#' data(TADs_Domains.gnr)
#' hicLst <- ImportHiC(
#'   file      = h5_path,
#'   hicResolution       = binSize,
#'   chromSizes = data.frame(seqnames = c("2L"), 
#'   seqlengths = c(23513712)),
#'   chrom_1   = c("2L")
#' )
#' hicLst <- BalanceHiC(hicLst)
#' hicLst <- OverExpectedHiC(hicLst)
#' # Index Beaf32
#' Beaf32_Index.gnr <- IndexFeatures(
#'   gRangeList = list(Beaf = Beaf32_Peaks.gnr),
#'   chromSizes = data.frame(seqnames = c("2L"), 
#'    seqlengths = c(23513712)),
#'   genomicConstraint = TADs_Domains.gnr,
#'   binSize = binSize
#' )
#' Beaf_Beaf.gni <- SearchPairs(indexAnchor = Beaf32_Index.gnr)
#' 
#' interactions_Ctrl.mtx_lst <- ExtractSubmatrix(
#'  genomicFeature = Beaf_Beaf.gni,
#'  hicLst = hicLst,
#'  referencePoint = "pf"
#' )
#' interactions_Ctrl.mtx_lst <- PrepareMtxList(
#'  matrices = interactions_Ctrl.mtx_lst
#' )
#' output_bgInterTAD = CompareToBackground(hicList = hicLst,
#'  matrices = interactions_Ctrl.mtx_lst,
#'  indexAnchor = Beaf32_Index.gnr,
#'  indexBait = Beaf32_Index.gnr,
#'  genomicConstraint = TADs_Domains.gnr,
#'  chromSizes = data.frame(seqnames = c("2L"), 
#'    seqlengths = c(23513712)),
#'  bg_type="inter_TAD"
#' )
CompareToBackground <- function(
    hicList = NULL,
    matrices = NULL,
    indexAnchor = NULL,
    indexBait = NULL,
    genomicConstraint = NULL,
    secondaryConst.var = NULL,
    chromSizes = NULL,
    n_background = NULL,
    areaFun="center",
    operationFun="mean",
    bg_type = NULL,
    cores = 1,
    verbose = FALSE,
    p_adj_method = 'BH',
    ...){

    targetCouples <- attributes(matrices)$interactions
    resolution <- attributes(matrices)$resolution

    if(bg_type == "random_anchors" || is.null(bg_type)){
        bg_couples <- .findRandomBins(targetCouples,
                        genomicConstraint = genomicConstraint,
                        indexAnchor=indexAnchor,
                        indexBait=indexBait,
                        chromSizes=chromSizes,
                        resolution=resolution,
                        cores=cores)
    }else{
        if(bg_type == "inter_TAD"){
        bg_couples <- .interTads(targetCouples,
                                genomicConstraint = genomicConstraint,
                                indexAnchor=indexAnchor,
                                indexBait=indexBait,
                                chromSizes=chromSizes,
                                resolution=resolution,
                                secondaryConst.var = NULL,
                                cores=cores)

        bg_couples <- .cleanCouples(bg_couples,matrices)
        }
        if(bg_type == "inter_compartment"){
        bg_couples <- .interTads(targetCouples,
                                genomicConstraint = genomicConstraint,
                                indexAnchor=indexAnchor,
                                indexBait=indexBait,
                                chromSizes=chromSizes,
                                resolution=resolution,
                                secondaryConst.var = secondaryConst.var,
                                cores=cores)
        bg_couples <- .cleanCouples(bg_couples,matrices)
        }
    }

    bg_couples <- .cleanCouples(bg_couples,matrices)
    if(length(bg_couples) < 500){
        ## If no tad is given just use random bins
        ## using indexAnchor with the same distance constraints
        ## would result in having the same couples as target logically
        bg_couples <- .findRandomBins(targetCouples,
                                    genomicConstraint=NULL,
                                    indexAnchor,
                                    indexBait,
                                    chromSizes,
                                    resolution)
        bg_couples <- .cleanCouples(bg_couples,matrices)
        if(verbose){
            message(length(bg_couples))
        }
    }

    if(!is.null(n_background) && n_background < length(bg_couples)){
        to_keep <- sample(seq(1,length(bg_couples)),
            size = n_background)
        bg_couples <- bg_couples[to_keep]
    }

    GenomeInfoDb::seqinfo(bg_couples) <- 
        GenomeInfoDb::Seqinfo(seqnames = as.character(chromSizes[[1]]),
                                    seqlengths = as.numeric(chromSizes[[2]]))
    bg_counts <- ExtractSubmatrix(genomicFeature = bg_couples,
                    hicLst = hicList,
                    hicResolution = resolution,
                    cores = cores,
                    referencePoint = "pf")
    bg_counts <- PrepareMtxList(bg_counts, ...)
    bg_quantifs <- GetQuantif(
        matrices = bg_counts,
        areaFun = areaFun,
        operationFun = operationFun
    ) |> as.numeric() |> 
        `names<-`(attr(bg_counts,"names"))

    target_quantifs <- GetQuantif(
        matrices = matrices,
        areaFun = areaFun,
        operationFun = operationFun
    ) |> as.numeric() |> 
        `names<-`(attr(matrices,"names"))
    z_output <- .computeZscore(
        df_bg = data.frame(
            names = attributes(bg_counts)$interactions$name,
            distance = attributes(bg_counts)$interactions$distance,
            counts = as.numeric(bg_quantifs)
        ),
        df_target = data.frame(
            names = attributes(matrices)$interactions$name,
            distance = attributes(matrices)$interactions$distance,
            counts = as.numeric(target_quantifs)
        ),
        method_adjust = p_adj_method
    )

    return(list(z.test = z_output, target_quantifs = target_quantifs,
        bg_quantifs = bg_quantifs))
}
