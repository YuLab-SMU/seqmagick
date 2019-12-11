##' combine 2 fasta files into 1
##'
##' 
##' @title fa_combine
##' @param file1 fasta file 1
##' @param file2 fasta file 2
##' @param outfile output file
##' @param type one of interleaved and sequential
##' @return BStringSet
##' @export
##' @author Guangchuang Yu
##' @examples
##' fa_file <- system.file("extdata/HA.fas", package="seqmagick")
##' x <- fa_read(fa_file)
##' fa1 <- tempfile(fileext=".fa")
##' fa2 <- tempfile(fileext=".fa")
##' fa_write(x[1:5], fa1)
##' fa_write(x[6:10], fa2)
##' fa_combine(fa1, fa2)
fa_combine <- function(file1, file2, outfile=NULL, type="interleaved") {
    fa1 <- readBStringSet(file1)
    fa2 <- readBStringSet(file2)
    fa <- c(fa1, fa2)

    if (!is.null(outfile)) {
        fa_write(fa, outfile, type)
    }
    invisible(fa)
}

