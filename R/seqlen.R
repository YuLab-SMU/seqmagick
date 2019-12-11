##' sequence length
##'
##'
##' @title seqlen
##' @param fasfile fasta file
##' @return numeric vector
##' @export
##' @author Guangchuang Yu
seqlen <- function(fasfile) {
  x <- fa_read(fasfile)
  yy <- toString(x)
  yy <- gsub("-", "", yy)
  yy <- gsub(" ", "", yy)
  xx <- strsplit(yy, ",")
  sapply(xx, nchar)
}

