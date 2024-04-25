# Data
data("Beaf32_Peaks.gnr")
data("TSS_Peaks.gnr")
data("TADs_Domains.gnr")
data("HiC_Ctrl.cmx_lst")
data("HiC_HS.cmx_lst")

# Global Variables
seqlengths.num <- c('2L'=23513712, '2R'=25286936)
chromSizes  <- data.frame(
    seqnames   = names(seqlengths.num ), 
    seqlengths = seqlengths.num
)

# Test BalanceHiC
bal_hic_tester <- BalanceHiC(HiC_Ctrl.cmx_lst,
    interactionType="all", method="ICE")
testthat::expect(
    ok = all(
        checkmate::checkList(
            bal_hic_tester,types = "ContactMatrix",
            null.ok = F),
        checkmate::checkTRUE(
            x = "cis" %in% as.character(
                attr(bal_hic_tester,"matricesKind")$type),
            na.ok = FALSE),
        checkmate::checkTRUE(
            x = "trans" %in% as.character(
                attr(bal_hic_tester,"matricesKind")$type),
            na.ok = FALSE),
        checkmate::checkNumeric(
            attributes(bal_hic_tester[[1]])$metadata$normalizer,
            null.ok = FALSE, any.missing = FALSE)),
    failure_message = "Error on BalanceHiC ICE normalization!")

bal_hic_tester <- BalanceHiC(
    HiC_Ctrl.cmx_lst, interactionType="all",
    method="VC")
testthat::expect(
    ok = all(
        checkmate::checkList(
            bal_hic_tester,types = "ContactMatrix",
            null.ok = F),
        checkmate::checkTRUE(
            x = "cis" %in% as.character(
                attr(bal_hic_tester,"matricesKind")$type),
            na.ok = FALSE),
        checkmate::checkTRUE(
            x = "trans" %in% as.character(
                attr(bal_hic_tester,"matricesKind")$type),
            na.ok = FALSE),
        checkmate::checkNumeric(
            attributes(bal_hic_tester[[1]])$metadata$normalizer,
            null.ok = FALSE, any.missing = FALSE)),
    failure_message = "Error on BalanceHiC VC normalization!")

bal_hic_tester <- BalanceHiC(
    HiC_Ctrl.cmx_lst, interactionType="all",
    method="VC_SQRT")

testthat::expect(
    ok = all(
        checkmate::checkList(
            bal_hic_tester,types = "ContactMatrix",
            null.ok = FALSE),
        checkmate::checkTRUE(
            x = "cis" %in% as.character(
                attr(bal_hic_tester,"matricesKind")$type),
            na.ok = FALSE),
        checkmate::checkTRUE(
            x = "trans" %in% as.character(
                attr(bal_hic_tester,"matricesKind")$type),
            na.ok = FALSE),
        checkmate::checkNumeric(
            attributes(bal_hic_tester[[1]])$metadata$normalizer,
            null.ok = FALSE, any.missing = FALSE)),
    failure_message = "Error on BalanceHiC VC normalization!")

HiC_Ctrl.cmx_lst <- BalanceHiC(HiC_Ctrl.cmx_lst, interactionType="cis")
HiC_HS.cmx_lst   <- BalanceHiC(HiC_HS.cmx_lst,   interactionType="cis")

# Test OverExpectedHiC
HiC_Ctrl.cmx_lst <- OverExpectedHiC(HiC_Ctrl.cmx_lst)
testthat::expect_type(
    attributes(
        HiC_Ctrl.cmx_lst[[1]])$metadata$expected,
    type = "double")
HiC_HS.cmx_lst   <- OverExpectedHiC(HiC_HS.cmx_lst, cores = 2)
testthat::expect_type(
    attributes(
        HiC_HS.cmx_lst[[1]])$metadata$expected,
    type = "double")

# Test SwitchMatrix
testthat::expect_identical(
    object = SwitchMatrix(HiC_Ctrl.cmx_lst,matrixKind = "norm")[[1]]@matrix@x,
    expected = (HiC_Ctrl.cmx_lst[[1]]@metadata$observed*
        HiC_Ctrl.cmx_lst[[1]]@metadata$normalizer)
)

testthat::expect_equal(
    object = SwitchMatrix(HiC_Ctrl.cmx_lst,matrixKind = "o/e")[[1]]@matrix@x,
    expected = ((HiC_Ctrl.cmx_lst[[1]]@metadata$observed*
        HiC_Ctrl.cmx_lst[[1]]@metadata$normalizer) /
        HiC_Ctrl.cmx_lst[[1]]@metadata$expected)
)

# Test IndexFeatures
Beaf32_Index.gnr <- IndexFeatures(
    gRangeList        = list(Beaf=Beaf32_Peaks.gnr), 
    genomicConstraint        = TADs_Domains.gnr,
    chromSizes         = chromSizes,
    binSize           = 100000,
    method            = "max",
    metadataColName = "score",
    cores             = 2,
    verbose           = TRUE
)
testthat::expect_match(
    Beaf32_Index.gnr$bin,
    regexp = "[2L|2R]:[1-9]",
    all = TRUE)
testthat::expect_match(
    Beaf32_Index.gnr$constraint,
    regexp = "Tad_[1-9]")

TSS_Index.gnr <- IndexFeatures(
    gRangeList        = list(TSS=TSS_Peaks.gnr), 
    genomicConstraint        = TADs_Domains.gnr,
    chromSizes         = chromSizes,
    binSize           = 100000,
    method            = "max",
    metadataColName = "score",
    cores             = 2,
    verbose           = TRUE
)

# Test SearchPairs
Beaf_TSS.gni <- SearchPairs(
    indexAnchor = Beaf32_Index.gnr,
    indexBait   = TSS_Index.gnr,
    minDist     = NULL, 
    maxDist     = NULL,
    cores       = 1,
    verbose     = TRUE
)

testthat::expect(
    ok = all(
        checkmate::checkClass(
            x = Beaf_TSS.gni,
            classes = "GInteractions",
            null.ok = FALSE
        ),
        checkmate::checkTRUE(
            x = all(
                Beaf_TSS.gni$anchor.Beaf.bln & Beaf_TSS.gni$bait.TSS.bln
            ),
            na.ok = FALSE
        )
    ), 
    failure_message = "SearchPairs didn't produce
        correct GInteractions object"
)
testthat::expect(ok = 
    checkmate::checkClass(
        SearchPairs(
            indexAnchor = Beaf32_Index.gnr,
            minDist     = "1", 
            maxDist     = "1MB",
            cores       = 1,
            verbose     = TRUE
        ),
    classes = "GInteractions",null.ok = FALSE),
failure_message = "SearchPairs potential call to GenomicSystem failing!")

testthat::expect(ok = 
    checkmate::checkClass(
        SearchPairs(
            indexAnchor = Beaf32_Index.gnr,
            indexBait   = TSS_Index.gnr,
            minDist     = 1, 
            maxDist     = 1000000,
            cores       = 2,
            verbose     = TRUE
        ),
    classes = "GInteractions",null.ok = FALSE),
failure_message = "SearchPairs parallel computing failing!")



# Test ExtractSubmatrix
submatrixPF_Ctrl.mtx_lst <- ExtractSubmatrix(
    genomicFeature         = Beaf_TSS.gni,
    hicLst        = HiC_Ctrl.cmx_lst,
    hicResolution            = NULL,
    referencePoint = "pf",
    matriceDim     = 21,
    cores          = 1,
    verbose        = TRUE
)
testthat::expect(
    ok = .validSubmatrices(submatrixPF_Ctrl.mtx_lst),
    failure_message = "ExtractSubmatrices failed in pf mode")

submatrixPF_HS.mtx_lst <- ExtractSubmatrix(
    genomicFeature         = Beaf_TSS.gni,
    hicLst        = HiC_HS.cmx_lst,
    hicResolution            = NULL,
    referencePoint = "pf",
    matriceDim     = 21,
    cores          = 1,
    verbose        = TRUE
)
submatrixRF_Ctrl.mtx_lst <- ExtractSubmatrix(
    genomicFeature         = Beaf_TSS.gni,
    hicLst        = HiC_Ctrl.cmx_lst,
    hicResolution            = NULL,
    referencePoint = "rf",
    matriceDim     = 21,
    cores          = 1,
    verbose        = TRUE
)
testthat::expect(
    ok = .validSubmatrices(submatrixRF_Ctrl.mtx_lst),
    failure_message = "ExtractSubmatrices failed in pf mode")
testthat::expect(
    ok = .validSubmatrices(
        ExtractSubmatrix(
            genomicFeature         = TADs_Domains.gnr,
            hicLst        = HiC_Ctrl.cmx_lst,
            referencePoint = "rf",
            matriceDim     = 101,
            cores          = 2,
            verbose        = FALSE
            )
        ),
    failure_message = "ExtractSubmatrices failed: look into ResizeMatrix")


# Test FilterInteractions
targets <- list(
    anchor.Beaf.name = c("Beaf32_204"),
    bait.TSS.name    = c("FBgn0264943")
)
selectionFun = function(){
    Reduce(union, list(anchor.Beaf.name, bait.TSS.name) )
}
testthat::expect(
    ok = .validSubmatrices(
        FilterInteractions(
            matrices      = submatrixPF_Ctrl.mtx_lst,
            targets        = targets,
            selectionFun     = selectionFun)
            ),
    failure_message = "FilterInteractions with anchor and bait names fails")

testthat::expect(
    ok = checkmate::checkList(
        x = FilterInteractions(
            genomicInteractions =
                attributes(submatrixPF_Ctrl.mtx_lst)$interactions,
            targets        = targets,
            selectionFun     = NULL
            ), types = 'numeric', null.ok = FALSE,
            names = "named", len = length(targets)
        ),
    failure_message = "FilterInteractions without submatrices fails!")

targets <- list(interactions = 
    attributes(submatrixPF_Ctrl.mtx_lst)$interactions[1:2])


testthat::expect(ok =
    checkmate::checkNumeric(
        x = GetQuantif(
            matrices  = submatrixRF_Ctrl.mtx_lst,
            areaFun      = "center",
            operationFun = "mean"
        ), null.ok = FALSE,
        any.missing = FALSE,
        len = length(submatrixRF_Ctrl.mtx_lst),
        names = 'named')
    , failure_message = "GetQuantif failed for rf submatrices"
)
testthat::expect(ok =
    checkmate::checkNumeric(
        x = GetQuantif(
            matrices  = submatrixPF_Ctrl.mtx_lst,
            areaFun      = "center",
            operationFun = "mean"
        ), null.ok = FALSE,
        any.missing = FALSE,
        len = length(submatrixPF_Ctrl.mtx_lst),
        names = 'named')
    , failure_message = "GetQuantif failed for pf submatrices"
)

# Test OrientateMatrix
testthat::expect(ok = 
    .validSubmatrices(OrientateMatrix(submatrixPF_Ctrl.mtx_lst)),
    failure_message = "OrientateMatrix failed"
)
testthat::expect_message(
    object = PrepareMtxList(
        submatrixPF_Ctrl.mtx_lst,orientate = TRUE),
        regexp = "[1-9] matrices are oriented"
)
testthat::expect(ok = 
    .validSubmatrices(
        PrepareMtxList(
            matrices = submatrixPF_HS.mtx_lst,
            transFun    = "qtl",
            minDist = 1,
            maxDist = "5Mb",
            rm0      = FALSE
        )),
    failure_message = "quantilization failed in PrepareMtxList")

# Test Aggregation
testthat::expect(ok = 
    checkmate::checkClass(x = 
        Aggregation(
            matrices = submatrixRF_Ctrl.mtx_lst, 
            aggFun      = "sum"
        ),classes="matrix"),
    failure_message = "Aggregation failed")

diffAggreg.mtx <- Aggregation(
    ctrlMatrices    = submatrixPF_Ctrl.mtx_lst,
    matrices        = submatrixPF_HS.mtx_lst,
    aggFun             = "mean",
    diffFun            = "substraction",
    scaleCorrection = TRUE,
        correctionArea  =  list(
        i = c(1:30),
        j = c(72:101)
    ),
    statCompare = TRUE
)
testthat::expect(
    ok = all(
        checkmate::checkClass(
            x = diffAggreg.mtx,
            classes = "matrix",
            null.ok = FALSE
        ),
        checkmate::checkTRUE(
            x = "differentialMethod"%in%names(attributes(diffAggreg.mtx)),
            na.ok = FALSE
        ),
        checkmate::checkTRUE(
            x = "correctedFact"%in%names(attributes(diffAggreg.mtx)),
            na.ok = FALSE
        )
    ),
    failure_message = "differential aggregation with correction area failed"
)