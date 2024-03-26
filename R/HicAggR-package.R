#' HicAggR
#'
#' HicAggR is a package that allows to integrate 1D genomics data with
#' 3D genomics data.
#' This package provides a set of functions useful in the analysis of 3D
#' genomic interactions. It includes the import of standard HiC data
#' formats into R and HiC normalisation procedures. The main objective of
#' this package is to facilitate the visualization and quantification of the
#' analysis of HiC contacts through aggregation.
#' The package also provides options to import externally normalized HiC data
#' to perform an in-depth analysis by investigating genome wide interactions
#' between features of interest.
#' The package can use 1D genomics data (such as annotation data or peaks
#' of features of interest in a GRanges object) to form potential couples
#' between features of interest under user specified constraints (TADs
#' coordinates or fixed minimum/max distance). Using these formed couples
#' it can extract HiC contact data for these specific couples. The submatrices
#' extracted as such can then be aggregated to summarize genome-wide
#' interactions, to perform per submatrix operation or compare between
#' conditions.
#' It also allows to identify couples with significantly signal-enriched
#' pixels at specific positions (such as the central pixel) relative to
#' background less plausible couples.
#'
#' @author Nicolas Chanard
#' @author David Depierre
#' @author Robel A Tesfaye
#' @author Naomi Schickele
#' @author Refka Askri
#' @author Pascal Martin
#' @author Stéphane Schaack
#' @author Olivier Cuvier
#'
#' @keywords internal
"_PACKAGE"