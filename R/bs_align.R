##' alignment
##'
##' 
##' @title bs_align
##' @param x BStringSet or DNAStringSet
##' @param method alignment method
##' @param ... additional parameters
##' @return BStringSet
##' @importFrom muscle muscle
##' @importFrom Biostrings DNAStringSet
##' @importFrom Biostrings BStringSet
##' @export
##' @author Guangchuang Yu
bs_align <- function(x, method="muscle", ...) {
    method <- match.arg(method, c("muscle")) ## currently, only muscle is supported
    if (method == 'muscle') 
        y <- muscle(DNAStringSet(x), ...)
    
    BStringSet(y)
}
