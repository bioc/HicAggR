# HicAggR 0.99.4

* removed CITATION file

# HicAggR 0.99.3

* NAMESPACE was generated with roxygen2 to define exportable functions.
* compareToBackground: to correct the skewedness of o/e values towards long distances,
    computation of z.scores is now calculated using residuals from a polynomial
    model that fits the background couples (log(counts)~distance).
* SearchPairs: added option to remove self interacting bins.
* all internal functions are now in utilities.R.
* Docs were reviewed.

# HicAggR 0.99.2

* Implemented import of corrected matrices for data in .hic and cool/mcool format.
* Implemented import of O/E matrix for data in .hic format.
* Implemented import of raw data in .h5 format.
* Added getInfos function to get info on a hic data (.hic, cool/mcool/h5 formats).
* Removed dependency to BSDA::z.test in compareToBackground.
* Removed dependency to InteractionSet and added it as package to import in NAMESPACE to remove the PackageStartUpMessages.
* Removed chatty package start up message and replaced it with nicer message.
* Some BiocCheck NOTES were also takedn into consideration: changing sapply to vapply etc.

# HicAggR 0.99.1

* Corrected with Bioconductor's reviews
* Added PrepareMtxList, importLoops, plotMultiAPA & compareToBackground functions
* Corrected bugs on seqlevels consistancy and name column for GRanges objects
* Encapsulated small and internal functions in utilities.R
* ExtractSubMatrix has option to remove duplicated submatrices
* Corrected quantilization operations in Aggregation
* Corrected over all code with suggestions from BiocCheck

# HicAggR 0.99.0

* Submitted to Bioconductor