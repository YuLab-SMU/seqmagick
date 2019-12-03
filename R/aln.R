##' sequence alignment
##'
##' 
##' @title bs_aln
##' @param x XStringSet object
##' @param method alignment method
##' @param ... additional parameter
##' @return aligned sequences, XStringSet object
##' @import muscle muscle
##' @importFrom Biostrings DNAStringSet
##' @importFrom Biostrings RNAStringSet
##' @importFrom Biostrings AAStringSet
##' @export
##' @author Guangchuang Yu
bs_aln <- function(x, method = 'muscle', ...) {
    method <- match.arg(method, c('muscle'))
    if (method == 'muscle') {
        res <- muscle(x, ...)
    }

    switch(class(res),
           DNAMultipleAlignment = DNAStringSet(res),
           RNAMultipleAlignment = RNAStringSet(res),
           AAMultipleAlignment = AAStringSet(res)
           )
}
