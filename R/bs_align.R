##' alignment
##'
##' 
##' @title bs_align
##' @param x DNAStringSet
##' @param method alignment method
##' @param ... additional parameters
##' @return DNAStringSet
##' @importFrom muscle muscle
##' @export
##' @author Guangchuang Yu
bs_align <- function(x, method="muscle", ...) {
    method <- match.arg(method, c("muscle")) ## currently, only muscle is supported
    if (method == 'muscle') 
        y <- muscle(x, ...)
    
    DNAStringSet(y)
}
