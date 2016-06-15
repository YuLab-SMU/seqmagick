##' convert fasta file to interleaved format
##'
##' 
##' @title fa_to_interleaved
##' @param file fasta file
##' @param outfile output file
##' @return NULL
##' @export
##' @author Guangchuang Yu
fa_to_interleaved <- function(file, outfile) {
    fa_write(readBStringSet(file), outfile, type="interleaved")
}

##' convert fasta file to sequential format
##'
##' 
##' @title fa_to_sequential
##' @param file fasta file
##' @param outfile output file
##' @return NULL
##' @export
##' @author Guangchuang Yu
fa_to_sequential <- function(file, outfile) {
    fa_write(readBStringSet(file), outfile, type="sequential")
}
