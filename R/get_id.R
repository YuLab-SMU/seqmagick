##' get id at specific position
##'
##'
##' @title get_id
##' @param x x
##' @param sep separator to split x
##' @param position id position
##' @return id
##' @export
##' @author guangchuang yu
get_id <- function(x, sep=" ", position) {
    res <- tokenize(x, sep)
    sapply(res, function(x) x[position])
}
