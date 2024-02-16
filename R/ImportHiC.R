#' Import Hic data
#'
#' ImportHiC
#' @description Import ..hic, .cool, .mcool or .bedpe data
#' @param file <GRanges or Pairs[GRanges] or GInteractions>:
#'  The genomic feature on which compute the extraction of HiC submatrix.
#'  Extension should be .hic, .cool, .mcool, .h5, .hdf5, .HDF5 or .bedpe"
#'  assuming .h5 and .hdf5 are only for cool (not mcool).
#' @param hicResolution <numeric>: The HiC resolution.
#' @param chromSizes <data.frame>: A data.frame where first colum correspond
#'  to the chromosomes names, and the second column correspond to the
#'  chromosomes lengths in base pairs.
#' @param chrom_1 <numeric>: The seqnames of firsts chromosmes (rows in matrix).
#' @param chrom_2 <numeric>: The seqnames of second chromosmes (col in matrix).
#'  If is NULL is equal to chrom_1 (Defalt NULL)
#' @param cores <numerical> : An integer to specify the number of cores.
#'  (Default 1)
#' @param verbose <logical>: If TRUE show the progression in console.
#'  (Default FALSE)
#' @return A matrices list.
#' @examples
#' \donttest{
#'
#' # Prepare Temp Directory
#' options(timeout = 3600)
#' temp.dir <- file.path(tempdir(), "HIC_DATA")
#' dir.create(temp.dir)
#'
#' # Download .hic file
#' Hic.url <- paste0(
#'     "https://4dn-open-data-public.s3.amazonaws.com/",
#'     "fourfront-webprod/wfoutput/",
#'     "7386f953-8da9-47b0-acb2-931cba810544/4DNFIOTPSS3L.hic"
#' )
#' HicOutput.pth <- file.path(temp.dir, "Control_HIC.hic")
#' download.file(Hic.url, HicOutput.pth, method = "auto", extra = "-k")
#'
#' # Import .hic file
#' HiC_Ctrl.cmx_lst <- ImportHiC(
#'     file = HicOutput.pth,
#'     hicResolution = 100000,
#'     chrom_1 = c("2L", "2L", "2R"),
#'     chrom_2 = c("2L", "2R", "2R")
#' )
#'
#' # Download .mcool file
#' Mcool.url <- paste0(
#'     "https://4dn-open-data-public.s3.amazonaws.com/",
#'     "fourfront-webprod/wfoutput/",
#'     "4f1479a2-4226-4163-ba99-837f2c8f4ac0/4DNFI8DRD739.mcool"
#' )
#' McoolOutput.pth <- file.path(temp.dir, "HeatShock_HIC.mcool")
#' download.file(Mcool.url, McoolOutput.pth, method = "auto", extra = "-k")
#'
#' # Import .mcool file
#' HiC_HS.cmx_lst <- ImportHiC(
#'     file = McoolOutput.pth,
#'     hicResolution = 100000,
#'     chrom_1 = c("2L", "2L", "2R"),
#'     chrom_2 = c("2L", "2R", "2R")
#' )
#' }
#'

ImportHiC <- function(
    file = NULL, hicResolution = NULL, chromSizes = NULL, chrom_1 = NULL,
    chrom_2 = NULL, verbose = FALSE, cores = 1, 
    hic_norm="NONE",hic_matrix = "observed"
) {
    # Resolution Format
    options(scipen = 999)
    if (inherits(hicResolution, "character")) {
        hicResolution <- GenomicSystem(hicResolution)
    }
    # Chromosomes Format
    if (is.null(chrom_2)) {
        chrom_2 <- chrom_1
    } else if (length(chrom_1) != length(chrom_2)) {
        stop("chrom_1 and chrom_2 must have the same length")
    }
    if ("ALL" %in% toupper(chrom_1)){
        chrom_1 <- chrom_1[-which(toupper(chrom_1) == "ALL")]
        message("ALL removed from chrom_1")
    }
    if ("ALL" %in% toupper(chrom_2)){
        chrom_2 <- chrom_2[-which(toupper(chrom_2) == "ALL")]
        message("ALL removed from chrom_2")
    }
    chrom.chr <- c(chrom_1, chrom_2) |>
        unlist() |>
        unique()
    if (grepl(pattern = "chr", chrom.chr[1], fixed = TRUE)) {
        seqlevelsStyleHiC <- "UCSC"
    } else {
        seqlevelsStyleHiC <- "ensembl"
    }
    # This throws error when building vignettes, if chromSizes
    #  object is not supplied
    # line 115-120 in HicAggR.Rmd & 168-174 in InDepth.Rmd
    # # chromSizs needs to have colnames = c("name", "length")
    # colnames(chromSizes) = c("name", "length")
    # Get SeqInfo
    # These lines make no sense to me, why check "index" presence,
    # if we are getting the chromSizes from file anyways...
    if (GetFileExtension(file) == "hic") {
        if ("index" %in% colnames(chromSizes)) {
            chromSizes <- strawr::readHicChroms(file) |>
                dplyr::select(-"index")
        } else {
            chromSizes <- strawr::readHicChroms(file)
        }
        if("index" %in% colnames(chromSizes)){
            chromSizes <- chromSizes |>
                dplyr::select(-"index")}
        colnames(chromSizes) <- c("name", "length")
    } else if (GetFileExtension(file) %in%
        c("cool", "mcool", "HDF5", "hdf5", "h5")
    ) {
        # Define HDF5groups
        chr.group <- ifelse(
            GetFileExtension(file) %in% c("cool", "HDF5", "hdf5", "h5"),
            yes = "/chroms",
            no = paste("resolutions", hicResolution, "chroms", sep = "/")
        )
        # Get SeqInfo
        chromSizes <- data.frame(rhdf5::h5read(file, name = chr.group))
    } else if (GetFileExtension(file) == "bedpe" &&
        !is.null(chromSizes)
    ) {
        colnames(chromSizes) <- c("name", "length")
        hic.gnp <- rtracklayer::import(file, format = "bedpe")
        megaHic.dtf <- data.frame(
            chrom_1 = as.vector(hic.gnp@first@seqnames),
            i = ceiling(hic.gnp@first@ranges@start/hicResolution),
            chrom_2 = as.vector(hic.gnp@second@seqnames),
            j = ceiling(hic.gnp@second@ranges@start/hicResolution),
            counts = hic.gnp@elementMetadata$score
        )
    } else {
        stop("file must be .hic, .cool, .mcool, .hdf5, .HDF5 or .bedpe")
    }
    rownames(chromSizes) <- chromSizes$name
    # vignette building throws error InDepth.Rmd (177-183), 
    # data is in UCSC chrom_1 in Ensembl
    # gives a void chromSizes table
    # # Standardize chromSizes and chrom.chr
    # chromSizes <- dplyr::filter(
    #     chromSizes,
    #     chromSizes$name %in% chrom.chr
    # )
    # Standardize seqlevelsStyle of chromSizes according to
    # chrom.chr
    # This brings a problem when GetFileExtension(file) == "hic", 
    # sometimes it adds a first chrom called "All" with out chr, 
    # So I changed rownames(chromSizes)[1] to 
    # rownames(chromSizes)[length(rownames(chromSizes))]
    if (grepl("chr", rownames(chromSizes)[
        length(rownames(chromSizes))],fixed = TRUE) &&
        seqlevelsStyleHiC == "ensembl"
    ) {
        chromSizes$name <- unlist(lapply(
            strsplit(rownames(chromSizes),"chr"),`[[`, 2
        ))
        rownames(chromSizes) <- unlist(lapply(
            strsplit(rownames(chromSizes),"chr"),`[[`, 2
        ))
    } else if (!grepl("chr", rownames(chromSizes)[
        length(rownames(chromSizes))], fixed = TRUE) &&
        seqlevelsStyleHiC == "UCSC"
    ) {
        chromSizes$name <- paste0("chr", rownames(chromSizes))
        rownames(chromSizes) <- paste0("chr", rownames(chromSizes))
    }
    chromSizes <- dplyr::mutate(
        chromSizes,
        dimension = ceiling(chromSizes$length/hicResolution)
    )
    chrom.chr <- chrom.chr[chrom.chr %in% chromSizes$name]
    chrom_1 <- chrom_1[chrom_1 %in% chromSizes$name]
    chrom_2 <- chrom_2[chrom_2 %in% chromSizes$name]
    # Create Genome as GRanges
    binnedGenome.grn <- chromSizes |>
        dplyr::pull("length") |>
        stats::setNames(chromSizes$name) |>
        GenomicRanges::tileGenome(
            tilewidth = hicResolution,
            cut.last.tile.in.chrom = TRUE
        )
    GenomeInfoDb::seqlengths(binnedGenome.grn) <- chromSizes$length |>
        stats::setNames(chromSizes$name)
    chromComb.lst <- paste(chrom_1, chrom_2, sep = "_")
    matrixSymmetric.bln <- strsplit(chromComb.lst, "_") |>
        lapply(function(name.chr) {name.chr[[1]] == name.chr[[2]]}) |>
        unlist()
    matrixType.str <- chromComb.lst
    matrixType.str[which(matrixSymmetric.bln)] <- "cis"
    matrixType.str[which(!matrixSymmetric.bln)] <- "trans"
    matrixKind <- chromComb.lst
    matrixKind[which(matrixSymmetric.bln)] <- "U"
    matrixKind[which(!matrixSymmetric.bln)] <- NA
    attributes.tbl <- dplyr::bind_cols(
        name = chromComb.lst, type = matrixType.str, kind = matrixKind,
        symmetric = matrixSymmetric.bln
    )
    # Dump file
    multicoreParam <- MakeParallelParam(
        cores = cores,
        verbose = verbose
    )
    hic.lst_cmx <- BiocParallel::bplapply(
        BPPARAM = multicoreParam, seq_along(chromComb.lst),
        function(ele.ndx) {
            # Chromosomes
            ele.lst <- unlist(strsplit(chromComb.lst[[ele.ndx]], "_"))
            chrom_1 <- ele.lst[[1]]
            chrom_2 <- ele.lst[[2]]
            # Dimension
            dims.num <- ele.lst |>
                lapply(
                    function(chrom) {
                        dplyr::filter(
                            chromSizes,
                            chromSizes$name == chrom) |>
                        dplyr::pull("dimension")
                    }
                ) |>
                unlist()
            if (GetFileExtension(file) == "hic") {
                # Read .hic file
                hic.dtf <- strawr::straw(
                    hic_norm,
                    file,
                    chrom_1,
                    chrom_2,
                    "BP",
                    hicResolution,
                    hic_matrix
                )
                hic.dtf$j <- ceiling((hic.dtf$y + 1)/hicResolution)
                hic.dtf$i <- ceiling((hic.dtf$x + 1)/hicResolution)
            } else if (GetFileExtension(file) %in%
                c("cool", "mcool", "HDF5", "hdf5", "h5")
            ) {
                # Define HDF5groups
                indexes.group <- ifelse(
                    GetFileExtension(file) %in%
                        c("cool", "HDF5", "hdf5", "h5"),
                    yes = "/indexes",
                    no = paste("resolutions",hicResolution,"indexes",sep = "/")
                )
                pixels.group <- ifelse(
                    GetFileExtension(file) %in%
                        c("cool", "HDF5", "hdf5", "h5"),
                    yes = "/pixels",
                    no = paste("resolutions", hicResolution,"pixels",sep = "/")
                )
                # Define start and end of chromosomes
                ends.ndx <- chromSizes$dimension |>
                    cumsum() |>
                    stats::setNames(chromSizes$name)
                starts.ndx <- 1 + c(0, ends.ndx[-length(ends.ndx)]) |>
                    stats::setNames(chromSizes$name)
                # Read .mcool file
                bin1.ndx <- as.vector(rhdf5::h5read(
                    file,
                    name = paste(indexes.group, "bin1_offset", sep = "/"),
                    index = list(starts.ndx[chrom_1]:ends.ndx[chrom_1])
                ))
                slice.num <- sum(
                    bin1.ndx[-1] - bin1.ndx[-length(bin1.ndx)]
                    ) - 1
                chunk.num <- seq(bin1.ndx[1] + 1, bin1.ndx[1] + 1 + slice.num)
                hic.dtf <- data.frame(
                    i = as.vector(rhdf5::h5read(
                        file,
                        name = paste(pixels.group, "bin1_id", sep = "/"),
                        index = list(chunk.num)
                    )) + 1,
                    j = as.vector(rhdf5::h5read(
                        file,
                        name = paste(pixels.group, "bin2_id", sep = "/"),
                        index = list(chunk.num)
                    )) + 1,
                    counts = as.vector(rhdf5::h5read(
                        file,
                        name = paste(pixels.group, "count", sep = "/"),
                        index = list(chunk.num)
                    ))
                )
                filter.bin2 <- hic.dtf$j %in%
                    starts.ndx[chrom_2]:ends.ndx[chrom_2]
                hic.dtf <- hic.dtf[filter.bin2, ]
                hic.dtf <- dplyr::mutate(
                    hic.dtf,
                    i = hic.dtf$i - starts.ndx[chrom_1] + 1
                )
                hic.dtf <- dplyr::mutate(
                    hic.dtf,
                    j = hic.dtf$j - starts.ndx[chrom_2] + 1
                )
            } else if (GetFileExtension(file) == "bedpe") {
                hic.dtf <- dplyr::filter(
                    megaHic.dtf,
                    # https://github.com/tidyverse/dplyr/issues/3139
                    # conflicting variable name
                    megaHic.dtf$chrom_1 == (!!chrom_1) &
                    megaHic.dtf$chrom_2 == (!!chrom_2)
                )
            }
            # Create Contact matrix
            hic.spm <- Matrix::sparseMatrix(
                i = hic.dtf$i,
                j = hic.dtf$j,
                x = hic.dtf$counts,
                dims = dims.num
            )
            row.regions <- binnedGenome.grn[which(
                as.vector(binnedGenome.grn@seqnames) == chrom_1
            )]
            col.regions <- binnedGenome.grn[which(
                as.vector(binnedGenome.grn@seqnames) == chrom_2
            )]
            hic <- InteractionSet::ContactMatrix(
                hic.spm,
                row.regions,
                col.regions
            )
            # Metadata
            hic@metadata <- dplyr::filter(
                attributes.tbl,
                attributes.tbl$name == paste(ele.lst, collapse = "_")
            ) |>
            tibble::add_column(resolution = hicResolution) |>
            as.list()
            if(hic_norm!="NONE"){
                hic@metadata <- append(hic@metadata,
                list(observed = NULL,
                normalizer = NULL,
                mtx = hic_norm))
                }
            if(hic_matrix!="observed"){
                hic@metadata <- append(hic@metadata,
                list(expected=hic_matrix))
                }
            
            return(hic)
        }
    )
    if(hic_matrix!="observed"){
        # Add attributes
        hic.lst_cmx <- hic.lst_cmx |>
            stats::setNames(chromComb.lst) |>
            AddAttr(
                attrs = list(
                    resolution = hicResolution,
                    chromSize = tibble::as_tibble(chromSizes),
                    matricesKind = attributes.tbl,
                    mtx = hic_matrix
                )
            )
    }else{
        # Add attributes
        hic.lst_cmx <- hic.lst_cmx |>
            stats::setNames(chromComb.lst) |>
            AddAttr(
                attrs = list(
                    resolution = hicResolution,
                    chromSize = tibble::as_tibble(chromSizes),
                    matricesKind = attributes.tbl,
                    mtx = "obs"
                )
            )

    }
    return(hic.lst_cmx)
}
