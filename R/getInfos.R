#' get some basic information about your hic file
#' 
#' getInfos
#'
#' @param file <character> path to your file.
#' @param printInfos <logical> print info on console? (Default: TRUE)
#' @param returnInfos <logical> return info? (Default: FALSE)
#'
#' @return list of characters if `returnInfos` is TRUE.
#' @importFrom rlang .data
#' @importFrom stringr str_split str_remove str_c str_replace str_detect
#' @importFrom rhdf5 h5ls h5read h5closeAll
#' @importFrom strawr readHicBpResolutions readHicChroms readHicNormTypes
#' @export
#' @examples
#' h5_path <- system.file("extdata",
#'     "HiC_Ctrl_2L_100Kres.h5",
#'     package = "HicAggR", mustWork = TRUE
#' )
#' getInfos(h5_path)
#'

getInfos <- function(file=NULL, printInfos = TRUE, returnInfos = FALSE) {
    options(scipen = 999)
    #
    if(file.exists(file)){
        # strawr doesn't like the tildes, this line will correct that.
        file <- normalizePath(file)
    }else{
        stop("File is not found!")
    }
    file_name <- stringr::str_split(file, "/", simplify = TRUE)[
        length(stringr::str_split(file, "/", simplify = TRUE))]
    file_path <- stringr::str_remove(file, paste0(file_name,"$"))
    if (! stringr::str_detect(file_path, "^[\\.]{0,1}/")) {
        file_path <- paste0("./",file_path)
    }
    if(GetFileExtension(file) == "hic"){
        file_res <- as.character(sort(
            strawr::readHicBpResolutions(file), decreasing = TRUE))
        
        chr_infos <- strawr::readHicChroms(file)
        file_chroms <- chr_infos[which(chr_infos$name != "ALL"),"name"]
        chrom_lengths <- chr_infos[which(chr_infos$name != "ALL"),"length"]
        file_norms <- strawr::readHicNormTypes(file)
        file_seqStyle <- ifelse(
            any(stringr::str_detect(file_chroms,"^chr")),
            "UCSC",
            "Ensembl"
        )
        file_units <- "should be <BP/FRAG>"
    } else {
        #
    file_head <- rhdf5::h5ls(file = file,
        recursive = FALSE, datasetinfo = FALSE)  # [['group']]
    rhdf5::h5closeAll()
    #
    if (dim(file_head)[1] == 1) {
        # for mcool file :
        file_res <- as.character(sort(
            as.integer(names(rhdf5::h5read(file = file,
            name = paste0(file_head[["group"]],
                file_head[["name"]]),
            recursive = FALSE,
            datasetinfo = FALSE))),
            decreasing = TRUE))
        rhdf5::h5closeAll()
        file_hNames <- names(rhdf5::h5read(file,
            name = paste0(file_head[["group"]],
                file_head[["name"]],
                "/", file_res[1]),
            recursive = FALSE,
            datasetinfo = FALSE))
        rhdf5::h5closeAll()
        if (length(intersect(file_hNames,
            c("bins", "chroms",
                "indexes", "pixels"))) !=
            4) {
            stop("invalid format")
        } else {
            file_normsLst <- rhdf5::h5read(file = file,
                name = paste0(file_head[["group"]],
                  file_head[["name"]],
                  "/", file_res[1],
                  "/bins"), recursive = FALSE)
            rhdf5::h5closeAll()
            file_norms <- stringr::str_replace(
                names(file_normsLst)[which(!names(file_normsLst) %in%
                c("chrom", "end",
                  "start"))], "weight",
                "ICE_balanced")
            file_seqStyle <- ifelse(any(stringr::str_detect(file_normsLst$chrom,
                "^chr")), "UCSC",
                "Ensembl")
            file_units <- "to be checked"
            file_chroms <- file_normsLst$chrom
            if (is.factor(file_chroms)) {
                file_chroms <- levels(file_chroms)
            }
            chrom_lengths <- as.data.frame(rhdf5::h5read(file = file,
                                           name = paste0(file_head[["group"]],
                                                         file_head[["name"]],
                                                         "/", file_res[1],
                                                         "/chroms")))$length
            rhdf5::h5closeAll()
        }
    } else if (dim(file_head)[1] ==
        4 && unique(file_head[["group"]] ==
        "/") && length(intersect(file_head[["name"]],
        c("bins", "chroms", "indexes",
            "pixels"))) == 4) {
        # for cool file :
        file_normsLst <- rhdf5::h5read(file = file,
            name = "/bins", recursive = FALSE)
        rhdf5::h5closeAll()
        file_norms <- stringr::str_replace(
            names(file_normsLst)[which(!names(file_normsLst) %in%
            c("chrom", "end", "start"))],
            "weight", "ICE_balanced")
        file_res <- as.character(file_normsLst$end[1] -
            file_normsLst$start[1])
        file_seqStyle <- ifelse(any(stringr::str_detect(file_normsLst$chrom,
            "^chr")), "UCSC", "Ensembl")
        file_units <- "to be checked"
        file_chroms <- file_normsLst$chrom
        if (is.factor(file_chroms)) {
            file_chroms <- levels(file_chroms)
        }
        chrom_lengths <- as.data.frame(rhdf5::h5read(file = file,
                                       name = "/chroms"))$length
    } else if (GetFileExtension(file) ==
        "h5") {
        # for h5 :
        file_norms <- ifelse("correction_factors" %in% file_head,
            "Correction factors are included in the h5",
            "No correction found!")
        bins <- data.frame(rhdf5::h5read(file, name = "intervals"))
        file_res <- bins$end_list[1]-bins$start_list[1]
        file_seqStyle <- ifelse(any(stringr::str_detect(bins$chr_list,
            "^chr")), "UCSC", "Ensembl")
        file_units <- "to be checked"
        bins <- bins |>
            dplyr::group_by(.data$chr_list) |> 
            dplyr::summarise("length" = max(.data$end_list))|>
            dplyr::rename("name"="chr_list")
        file_chroms <- as.character(bins$name)
        if (is.factor(file_chroms)) {
            file_chroms <- levels(file_chroms)
        }
        chrom_lengths <- as.numeric(bins$length)
    } else {
        stop("invalid format")
    }
    }
    
    if (printInfos) {
        message("file name: ",
                paste(file_name, sep = " "))
        message("file path: ", paste(file_path,
            sep = " "))
        message("resolution(s): ",
                paste(file_res, collapse = " "))
        message("normalization(s): ",
                paste(file_norms, collapse = " "))
        message("unit(s): ", paste(file_units,
                                   collapse = " "))
        message("seqStyle: ", paste(file_seqStyle,
                                    collapse = " "))
        message("chromosom(s): ",
                paste(ifelse(length(file_chroms) <=
                25, stringr::str_c(file_chroms,
                collapse = " "),
                paste0(stringr::str_c(head(file_chroms,
                  n = 4), collapse = " ; "),
                  " (...) ", stringr::str_c(tail(file_chroms,
                    n = 4), collapse = " ; "),
                  " (", length(file_chroms),
                  " total)")),
                collapse = " "))
        message("chromosom length(s): ",
                paste(ifelse(length(chrom_lengths) <=
                25, stringr::str_c(chrom_lengths,
                collapse = " "),
                paste0(stringr::str_c(head(chrom_lengths,
                  n = 4), collapse = " ; "),
                  " (...) ", stringr::str_c(tail(chrom_lengths,
                    n = 4), collapse = " ; "),
                  " (", length(chrom_lengths),
                  " total)")),
                collapse = " "))
    }
    if (returnInfos) {
        return(list(name = file_name,
            path = file_path, res = file_res,
            norms = file_norms,
            units = file_units,
            seqStyle = file_seqStyle,
            chroms = file_chroms,
            lengths = chrom_lengths))
    }
}
