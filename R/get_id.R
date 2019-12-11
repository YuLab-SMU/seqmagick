##' get id at specific position
##'
##'
##' @title get_id
##' @param x sequence description line
##' @param sep separator to split x
##' @param position id position
##' @return id
##' @export
##' @author Guangchuang Yu
##' @examples
##' fa_file <- system.file("extdata/HA.fas", package="seqmagick")
##' x <- fa_read(fa_file)
##' get_id(names(x)[1:5], sep = " ", position=1)
get_id <- function(x, sep=" ", position) {
    res <- tokenize(x, sep)
    sapply(res, function(x) x[position])
}
