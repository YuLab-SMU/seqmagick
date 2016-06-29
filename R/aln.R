##' sequence alignment
##'
##' 
##' @title bs_aln
##' @param x XStringSet object
##' @param method alignment method
##' @param ... additional parameter
##' @return aligned sequences, XStringSet object
##' @import muscle
##' @export
##' @author Guangchuang Yu
bs_aln <- function(x, method = 'muscle', ...) {
    method <- match.arg(method, c('muscle'))
    if (method == 'muscle') {
        res <- muscle2(x, ...)
    }
    return(res)
}
