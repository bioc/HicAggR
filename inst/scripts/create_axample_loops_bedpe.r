library(tidyverse)
## HiCCUPS loops from GSE63525, called on HeLa cells
readr::read_tsv(paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/",
        "GSE63nnn/GSE63525/suppl/GSE63525_HeLa_HiCCUPS_looplist.txt.gz"),
        col_names = T) %>% 
    dplyr::mutate(start1 = plyr::round_any(((x2+x1)/2)-2500,5000,f = ceiling)) %>% 
    dplyr::mutate(end1 = plyr::round_any(((x2+x1)/2)+2500,5000,f = ceiling)) %>% 
    dplyr::mutate(start2 = plyr::round_any(((y2+y1)/2)-2500,5000,f = ceiling)) %>% 
    dplyr::mutate(end2 = plyr::round_any(((y2+y1)/2)+2500,5000,f = ceiling)) %>%
    dplyr::mutate(name=".") %>% 
    dplyr::mutate(score=".") %>% 
    dplyr::mutate(strand1="*") %>% 
    dplyr::mutate(strand2="*") %>% 
    dplyr::select(c(chr1,start1,end1,
        chr2,start2,end2,
        name,score,strand1,strand2)) %>% 
    dplyr::filter(chr1%in%c(1,2))
    readr::write_tsv("inst/extdata/postprocessed_pixels_5000.bedpe",
        col_names = FALSE)