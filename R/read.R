##' read fasta file
##'
##' 
##' @title fa_read
##' @param file fasta file
##' @return BStringSet object
##' @importFrom Biostrings readBStringSet
##' @export
##' @author Guangchuang Yu
fa_read <- function(file) {
    readBStringSet(file)
}
