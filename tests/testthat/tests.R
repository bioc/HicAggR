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
BalanceHiC(HiC_Ctrl.cmx_lst, interactionType=c("cis", "trans"), method="ICE")
BalanceHiC(HiC_Ctrl.cmx_lst, interactionType=c("cis", "trans"), method="VC")
BalanceHiC(HiC_Ctrl.cmx_lst, interactionType=c("cis", "trans"), method="VC_SQRT")
BalanceHiC(HiC_Ctrl.cmx_lst, interactionType="trans", method="ICE")
BalanceHiC(HiC_Ctrl.cmx_lst, interactionType="trans", method="VC")
BalanceHiC(HiC_Ctrl.cmx_lst, interactionType="trans", method="VC_SQRT")
BalanceHiC(HiC_Ctrl.cmx_lst, interactionType="cis", method="ICE")
BalanceHiC(HiC_Ctrl.cmx_lst, interactionType="cis", method="VC")
BalanceHiC(HiC_Ctrl.cmx_lst, interactionType="cis", method="VC_SQRT")
BalanceHiC(HiC_Ctrl.cmx_lst, interactionType="all", method="ICE")
BalanceHiC(HiC_Ctrl.cmx_lst, interactionType="all", method="VC")
BalanceHiC(HiC_Ctrl.cmx_lst, interactionType="all", method="VC_SQRT")

HiC_Ctrl.cmx_lst <- BalanceHiC(HiC_Ctrl.cmx_lst, interactionType="cis")
HiC_HS.cmx_lst   <- BalanceHiC(HiC_HS.cmx_lst,   interactionType="cis")

# Test OverExpectedHiC
HiC_Ctrl.cmx_lst <- OverExpectedHiC(HiC_Ctrl.cmx_lst)
HiC_HS.cmx_lst   <- OverExpectedHiC(HiC_HS.cmx_lst, cores = 2)

# Test SwitchMatrix
SwitchMatrix(HiC_Ctrl.cmx_lst, matrixKind="norm")
SwitchMatrix(HiC_Ctrl.cmx_lst, matrixKind="o/e")

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
IndexFeatures(
    gRangeList        = list(TSS_1=TSS_Peaks.gnr, TSS_2=TSS_Peaks.gnr), 
    genomicConstraint        = NULL,
    chromSizes         = chromSizes,
    binSize           = 100000,
    cores             = 1,
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
SearchPairs(
    indexAnchor = Beaf32_Index.gnr,
    minDist     = "1", 
    maxDist     = "1MB",
    cores       = 1,
    verbose     = TRUE
)
SearchPairs(
    indexAnchor = Beaf32_Index.gnr,
    indexBait   = TSS_Index.gnr,
    minDist     = 1, 
    maxDist     = 1000000,
    cores       = 2,
    verbose     = TRUE
)


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
ExtractSubmatrix(
    genomicFeature         = TADs_Domains.gnr,
    hicLst        = HiC_Ctrl.cmx_lst,
    referencePoint = "rf",
    matriceDim     = 101,
    cores          = 2,
    verbose        = FALSE
)


# Test FilterInteractions
targets <- list(
    anchor.Beaf.name = c("Beaf32_113"),
    bait.TSS.name    = c("FBgn0267378")
)
selectionFun = function(){
    Reduce(union, list(anchor.Beaf.name, bait.TSS.name) )
}
FilterInteractions(
    matrices      = submatrixPF_Ctrl.mtx_lst,
    targets        = targets,
    selectionFun     = selectionFun
)
FilterInteractions(
  genomicInteractions = attributes(submatrixPF_Ctrl.mtx_lst)$interactions,
  targets        = targets,
  selectionFun     = NULL
)
targets <- list(interactions = attributes(submatrixPF_Ctrl.mtx_lst)$interactions[1:2])
FilterInteractions(
    genomicInteractions = attributes(submatrixPF_Ctrl.mtx_lst)$interactions,
    targets        = targets,
    selectionFun     = NULL
)
targets <- list(first = InteractionSet::anchors(attributes(submatrixPF_Ctrl.mtx_lst)$interactions)[["first"]][1:2])
FilterInteractions(
    genomicInteractions = attributes(submatrixPF_Ctrl.mtx_lst)$interactions,
    targets        = targets,
    selectionFun     = NULL
)

# Test GetQuantif
GetQuantif(
    matrices  = submatrixRF_Ctrl.mtx_lst,
    areaFun      = "center",
    operationFun = "mean"
)
GetQuantif(
    matrices  = submatrixPF_Ctrl.mtx_lst,
    areaFun      = "center",
    operationFun = "mean"
)

# Test OrientateMatrix
OrientateMatrix(submatrixPF_Ctrl.mtx_lst)

# Test Aggregation
Aggregation(
    matrices = submatrixRF_Ctrl.mtx_lst, 
    aggFun      = "sum",
    transFun    = "qtl", 
    rm0      = FALSE
    )
diffAggreg.mtx <- Aggregation(
    ctrlMatrices    = submatrixPF_Ctrl.mtx_lst,
    matrices        = submatrixPF_HS.mtx_lst,
    minDist             = 1,
    maxDist             = "5Mb",
    aggFun             = "mean",
    rm0             = FALSE,
    diffFun            = "substraction",
    scaleCorrection = TRUE,
        correctionArea  =  list(
        i = c(1:30),
        j = c(72:101)
    ),
    statCompare = TRUE
)

# Test ggAPA and PlotAPA
ggAPA(
    aggregatedMtx      = diffAggreg.mtx,
    title    = "APA",
    colMin   = 0,
    colMax   = 10,
    trim = 20,
    tails   = "both",
    blurPass = 1,
    stdev   = 0.5,
    colors  = NULL
)


ggAPA(
    aggregatedMtx      = diffAggreg.mtx,
    title    = "APA",
    colMid   = 0,
    trim = NULL,
    tails   = "both",
    blurPass = 1,
    stdev   = 0.5,
    colors  = viridis(6)
)

# Complete Tests
TrimOutliers(rnorm(1000))
GaussBox(kernScale="int")
Gauss(x=1, y=1)
GenomicSystem(1000000000,2)
GRange_1.grn <- StrToGRanges(c("chr1:1-100:+","chr2:400-500:-"))
GRange_2.grn <- StrToGRanges("chr1:10-50:*")
MergeGRanges(list(GRange_1.grn,GRange_2.grn), reduceRanges=TRUE, sortRanges=TRUE)
Hue(paletteLength=1)
Hue(paletteLength=2)
PadMtx(mtx=matrix(1:25,5,5), padSize=1, val=0, side=c('top','bot','right','left') )
ReduceRun(firstRle=rle(c("A","A","B")), secondRle=rle(c("A","B","B")), reduceMethod="paste", sep="_" )
MeanScale(rnorm(500,500))
BreakVector(x=rnorm(500,500), n=50, method="density")
BreakVector(x=rnorm(500,500), n=50, center=-500)
set.seed(123)
spMtx = as(matrix(floor(runif(7*13,0,2)),7,13), "dgCMatrix")
Rise0(spMtx=spMtx, indices=c(1,3,6,10,12))
Rise0(spMtx=spMtx, coords=data.frame(i=c(1,5,3), j=c(1,2,3) ) )
Rise0(spMtx=spMtx)
GetFileExtension(file="my/path/to/my/file.txt")