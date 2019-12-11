##' convert fasta file to interleaved format
##'
##' 
##' @title fa_to_interleaved
##' @rdname fa_conversion
##' @param file fasta file
##' @param outfile output file
##' @return None
##' @export
##' @author Guangchuang Yu
##' @examples
##' fa_file <- system.file("extdata/HA.fas", package="seqmagick")
##' fa1 <- tempfile(fileext = '.fa')
##' fa2 <- tempfile(fileext = '.fa')
##' fa_to_interleaved(fa_file, fa1)
##' fa_to_sequential(fa_file, fa2)
fa_to_interleaved <- function(file, outfile) {
    fa_write(readBStringSet(file), outfile, type="interleaved")
}

##' convert fasta file to sequential format
##'
##' 
##' @title fa_to_sequential
##' @rdname fa_conversion
##' @param file fasta file
##' @param outfile output file
##' @return None
##' @export
##' @author Guangchuang Yu
fa_to_sequential <- function(file, outfile) {
    fa_write(readBStringSet(file), outfile, type="sequential")
}
