
#' Computes z.test for each target couple over background couples.
#'
#' @param hicList <List[ContactMatrix][InteractionSet::ContactMatrix()]>: 
#' The HiC maps list.
#' @param matrices <list[matrix]>: The matrices list.
#' @param genomicConstraint <GRanges>: GRanges object of constraint regions. 
#'  If NULL chromosomes in chromSizes are used as constraints (Default NULL)
#' @param secondaryConst.var <character>: A string defining column name 
#' containing compartment information in the metadata of anchor and 
#' bait <GRanges> objects. (Default NULL)
#' @param chromSizes <data.frame>: A data.frame containing chromosomes 
#' names and lengths in base pairs.
#' @param n_background <integer> : Number of background couples to keep. 
#' (Default NULL)
#' @param indexAnchor <GRanges>: A first indexed GRanges object used as 
#' pairs anchor (must be indexed using `IndexFeatures()`).
#' @param indexBait <GRanges>: A second indexed GRanges object used as 
#' pairs bait (must be indexed using `IndexFeatures()`).
#' @param cores <integer> : Number of cores used. (Default 1)
#' @param areaFun <character or function>: A character or function that allows
#'  to extract an area from each matrix that composes the matrices list 
#' (Default "center"). Look at [GetQuantif] for more info.
#' @param operationFun <character or function>: A character or function 
#' specifying the operation applied to the selected area (Default "mean"). 
#' Look at [GetQuantif] for more info.
#' @param bg_type <character>: Type of background couples to generate.
#'\itemize{
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
#' \item "inter_compartment": If `secondaryConst.var` is not `NULL` and both 
#' indexAnchor and indexBait objects contain the provided variable name, 
#' then background couples are formed between anchors and baits located in 
#' different compartments.
#' \item "NULL": If `NULL`, `random_anchors` are set by default.
#' }
#' @param verbose <logical> details on progress? (Default: FALSE)
#' @param ... arguments to pass to [PrepareMtxList], inorder to treat 
#' background matrices.
#' @return returns a <list> object with the z.test output for each 
#' target couple, values for the target couples and values for the 
#' background couples.
#' 
#' @export
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom IRanges subsetByOverlaps
#' @importFrom BSDA z.test
#' @importFrom stats p.adjust sd
#' @importFrom broom tidy
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
#' output_bgInterTAD = compare_to_background(hicList = hicLst,
#'  matrices = interactions_Ctrl.mtx_lst,
#'  indexAnchor = Beaf32_Index.gnr,
#'  indexBait = Beaf32_Index.gnr,
#'  genomicConstraint = TADs_Domains.gnr,
#'  chromSizes = data.frame(seqnames = c("2L"), 
#'    seqlengths = c(23513712)),
#'  bg_type="inter_TAD"
#' )
#' 
compare_to_background <- function(hicList = NULL,
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
                    ...){

    targetCouples <- attributes(matrices)$interactions
    resolution <- attributes(matrices)$resolution
    .findRandomBins <- function(
        targetCouples,
        genomicConstraint=NULL,
        indexAnchor=NULL,
        indexBait=NULL,
        chromSizes,
        resolution,
        N=n_background,
        cores=cores){
        
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
            seqnames = chromSizes[, 1],
            ranges = IRanges::IRanges(
            start = rep(1, length(chromSizes[, 2])),
            end = chromSizes[, 2]
            ),
            strand = "*", name = chromSizes[, 1]
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
        if(!is.null(n_background)){
        background_pairs <- background_pairs[
            sample(seq(1,length(background_pairs)),size = N)]
        }
        if(verbose){
            message("Number of background pairs: ",
            length(background_pairs))
        }
        return(background_pairs)
    }
    
    .cleanCouples <- function(bg_couples,
                            matrices){
        # discard couples that are in target couples
        bg_couples <- bg_couples[
        which(!attr(bg_couples,"NAMES")%in%attr(matrices,"names"))]
        return(bg_couples)
    }
    
    
    .interTads <- function(targetCouples,
                            genomicConstraint=NULL,
                            indexAnchor=NULL,
                            indexBait=NULL,
                            chromSizes,
                            resolution,
                            N=n_background,
                            secondaryConst.var = NULL,
                            cores=cores){
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
        dplyr::mutate(tad_name = "constraint") |>
        dplyr::mutate(constraint = "seqnames")|>
        GenomicRanges::GRanges()
        
        indexBait.tmp <- indexBait |>
        as.data.frame() |>
        dplyr::mutate(tad_name = "constraint") |>
        dplyr::mutate(constraint = "seqnames")|>
        GenomicRanges::GRanges()
        
        bg_couples <- SearchPairs(
        indexAnchor = indexAnchor.tmp,
        indexBait = indexBait.tmp,
        minDist = dist_const[1],
        maxDist = dist_const[2],
        verbose = FALSE,
        cores = cores)
        
        diff.tad.couples <- as.data.frame(S4Vectors::mcols(bg_couples)) |>
        dplyr::filter("anchor.tad_name"!="bait.tad_name") |>
        rownames()
        bg_couples <- bg_couples[
        which(attr(bg_couples,"NAMES")%in%diff.tad.couples)
        ]
        if(verbose){
            message("Number of inter-TAD couples ",length(bg_couples))
        }
        if(!is.null(secondaryConst.var)){
        # discard couples in the same compartment
        # names of compartment info in couples
        compart.var1 <- paste0("anchor.",secondaryConst.var)
        compart.var2 <- paste0("bait.",secondaryConst.var)
        diff.compart.couples <- as.data.frame(S4Vectors::mcols(bg_couples)) |>
            dplyr::filter(compart.var1!=compart.var2) |>
            rownames()
        bg_couples <- bg_couples[
            which(attr(bg_couples,"names")%in%diff.compart.couples)]
        }
        
        bg_couples <- .cleanCouples(bg_couples,matrices)
        return(bg_couples)
    }
    
    
    if(bg_type == "random_anchors" || is.null(bg_type)){
        bg_couples <- .findRandomBins(targetCouples,
                        genomicConstraint = genomicConstraint,
                        indexAnchor=indexAnchor,
                        indexBait=indexBait,
                        chromSizes=chromSizes,
                        resolution=resolution,
                        N=n_background,
                        cores=cores)
    }else{
        if(bg_type == "inter_TAD"){
        bg_couples <- .interTads(targetCouples,
                                genomicConstraint = genomicConstraint,
                                indexAnchor=indexAnchor,
                                indexBait=indexBait,
                                chromSizes=chromSizes,
                                resolution=resolution,
                                N=n_background,
                                secondaryConst.var = NULL,
                                cores=cores)
        }
        if(bg_type == "inter_compartment"){
        bg_couples <- .interTads(targetCouples,
                                genomicConstraint = genomicConstraint,
                                indexAnchor=indexAnchor,
                                indexBait=indexBait,
                                chromSizes=chromSizes,
                                resolution=resolution,
                                N=n_background,
                                secondaryConst.var = secondaryConst.var,
                                cores=cores)
        }
    }
    
    bg_couples <- .cleanCouples(bg_couples,matrices)
    if(length(bg_couples) < 500){
        # If no tad is given just use random bins
        # using indexAnchor with the same distance constraints
        # would result in having the same couples as target logically
        # couples with random
        bg_couples <- .findRandomBins(targetCouples,
                                    genomicConstraint=NULL,
                                    indexAnchor,
                                    indexBait,
                                    chromSizes,
                                    resolution,
                                    n_background)
        bg_couples <- .cleanCouples(bg_couples,matrices)
        if(verbose){
            message(length(bg_couples))
        }
    }
    
    GenomeInfoDb::seqinfo(bg_couples) <- 
        GenomeInfoDb::Seqinfo(seqnames = chromSizes$seqnames,
                                    seqlengths = chromSizes$seqlengths)
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
    ) |> as.numeric()
    sigma.bg <- stats::sd(bg_quantifs)
    
    target_quantifs <- GetQuantif(
        matrices = matrices,
        areaFun = areaFun,
        operationFun = operationFun
    ) |> as.numeric() |> 
        `names<-`(attr(matrices,"names"))
    z.output <- lapply(target_quantifs, 
            function(i){broom::tidy(BSDA::z.test(x=bg_quantifs,
                                            sigma.x = sigma.bg,
                                            mu=i,
                                            alternative = "less"))}
            ) |> dplyr::bind_rows() |> 
        dplyr::mutate(names = names(target_quantifs)) |>
        dplyr::mutate(adj.p = stats::p.adjust(p = "p.value",
            method="bonferroni",n=length(target_quantifs)))
    return(list(z.test = z.output,target_quantifs=target_quantifs,
        bg_quantifs=bg_quantifs))
}
