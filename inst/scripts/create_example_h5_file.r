library(HicAggR)
## Get example data from package
data("HiC_Ctrl.cmx_lst")
## create h5 file and groups
rhdf5::h5createFile(file = "inst/extdata/Control_HIC_10k_2L.h5")
rhdf5::h5createGroup(file = "inst/extdata/Control_HIC_10k_2L.h5",
    group = "/intervals")
rhdf5::h5createGroup(file = "inst/extdata/Control_HIC_10k_2L.h5",
    group = "/matrix")
## Write matrix data
rhdf5::h5write(HiC_Ctrl.cmx_lst[["2L_2L"]]@matrix@x,
    file = "inst/extdata/Control_HIC_10k_2L.h5",
    name = "/matrix/data")
rhdf5::h5write(HiC_Ctrl.cmx_lst[["2L_2L"]]@matrix@i,
    file = "inst/extdata/Control_HIC_10k_2L.h5",
    name = "/matrix/indices")
rhdf5::h5write(HiC_Ctrl.cmx_lst[["2L_2L"]]@matrix@p,
    file = "inst/extdata/Control_HIC_10k_2L.h5",
    name = "/matrix/indptr")
rhdf5::h5write(dim(HiC_Ctrl.cmx_lst[["2L_2L"]]@matrix),
    file = "inst/extdata/Control_HIC_10k_2L.h5",
    name = "/matrix/shape")

## Write itervals iformation
rhdf5::h5write(as.character(
    GenomeInfoDb::seqnames(HiC_Ctrl.cmx_lst[["2L_2L"]]@regions)),
    file = "inst/extdata/Control_HIC_10k_2L.h5",
    name = "/intervals/chr_list")

rhdf5::h5write(GenomicRanges::end(HiC_Ctrl.cmx_lst[["2L_2L"]]@regions),
    file = "inst/extdata/Control_HIC_10k_2L.h5",
    name = "/intervals/end_list")
rhdf5::h5write(rep(1,times=length(HiC_Ctrl.cmx_lst[["2L_2L"]]@regions)),
    file = "inst/extdata/Control_HIC_10k_2L.h5",
    name = "/intervals/extra_list")
rhdf5::h5write(GenomicRanges::start(HiC_Ctrl.cmx_lst[["2L_2L"]]@regions)-1,
    file = "inst/extdata/Control_HIC_10k_2L.h5",
    name = "/intervals/start_list")

## write the coordinates of nan_bins
rhdf5::h5write(which(is.na(HiC_Ctrl.cmx_lst[["2L_2L"]]@matrix@x)),
    file = "inst/extdata/Control_HIC_10k_2L.h5",
    name = "/nan_bins")
rhdf5::h5closeAll()
## Check validity of data
ImportHiC("inst/extdata/Control_HIC_10k_2L.h5",
    hicResolution = 100000,
    chrom_1 = "2L")
GetInfo("inst/extdata/Control_HIC_10k_2L.h5")
